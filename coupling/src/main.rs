//! Acople operativo snowmelt-rs → rainflow sobre el Río Choapa en Cuncumén.
//!
//! snowmelt-rs resuelve la fase nival (acumulación, balance de energía,
//! ablación) sobre bandas de elevación y entrega el **aporte líquido**
//! (lluvia + derretimiento) por día. rainflow (GR4J) cierra el balance
//! suelo-escorrentía: lo toma como "precipitación" y produce el caudal.
//!
//! Compara, en split-sample, GR4J alimentado por:
//!   (a) precipitación cruda (sin nieve) — el baseline ingenuo, y
//!   (b) el aporte líquido de snowmelt — el modelo acoplado.
//!
//! En una cuenca nival el acoplado debe ganar holgado: la precipitación
//! cruda llega como caudal en invierno (cuando en realidad se acumula como
//! nieve), y falta en primavera (cuando derrite).
//!
//! Requiere un checkout de rainflow en `../rainflow`. Ejecutar con:
//! `cargo run -p snowmelt-cli --example couple_rainflow --features rainflow --release`

use std::path::Path;

use ndarray::Array2;
use rainflow_core::{DdsConfig, Gr4j, Objective, Optimizer, calibrate_gr4j, metrics};
use snowmelt_cli::{asc, forcing_csv};
use snowmelt_core::{AlbedoDecay, DegreeDayParams, Dem, EnergyBalanceParams, Forcing, SnowModel};

const BASE: &str = "validation/choapa-cuncumen/data";
const Z_REF: f64 = 3142.0;
const LATITUDE: f64 = -31.95;
const WARMUP: usize = 365;
// Límites de los parámetros GR4J (airGR): x1, x2, x3, x4.
const BOUNDS: [(f64, f64); 4] = [(1.0, 2000.0), (-10.0, 5.0), (1.0, 1000.0), (0.5, 10.0)];

/// Aporte líquido diario (lluvia + derretimiento) medio de cuenca, vía el
/// modelo nival de snowmelt-core sobre las bandas de elevación.
fn snowmelt_input() -> Vec<f64> {
    let dem_grid = asc::read(Path::new(&format!("{BASE}/bands_dem.asc"))).expect("bands_dem");
    let elevation = dem_grid.data.clone();
    let dem = Dem::new(dem_grid.data).expect("dem");
    let records = forcing_csv::read(Path::new(&format!("{BASE}/forcing.csv"))).expect("forcing");

    let params = DegreeDayParams {
        lapse_rate: -0.0075,
        energy_balance: Some(EnergyBalanceParams::default()),
        albedo_decay: Some(AlbedoDecay {
            tau_days: 9.0,
            albedo_min: 0.4,
            ..AlbedoDecay::default()
        }),
        ..DegreeDayParams::default()
    };
    let mut model = SnowModel::new(dem, params).expect("model");

    // Radiación potencial de cielo despejado por banda (terreno plano, el
    // pseudo-DEM tiene slopes ~0). Constante para una latitud y día del año;
    // aquí basta una aproximación estacional simple por día del año.
    records
        .iter()
        .map(|rec| {
            let doy = day_of_year(&rec.date);
            let rad = clear_sky_radiation(LATITUDE, doy, &elevation);
            let forcing = Forcing::Uniform {
                t_ref: rec.temp_c,
                z_ref: Z_REF,
                precip: rec.precip_mm,
            };
            let out = model
                .step_radiation(&forcing, Some(rad.view()), 1.0)
                .expect("step");
            let s = model.summarize(&out);
            s.mean_runoff // lluvia + derretimiento, mm/día
        })
        .collect()
}

/// Radiación de onda corta de tope de atmósfera modulada por transmitancia
/// (W m⁻², media diaria) para una latitud y día del año, terreno plano.
/// Suficiente para el término ETI sobre bandas sin aspecto.
fn clear_sky_radiation(lat_deg: f64, doy: u32, elevation: &Array2<f64>) -> Array2<f64> {
    let lat = lat_deg.to_radians();
    let decl = 0.409 * (2.0 * std::f64::consts::PI * doy as f64 / 365.0 - 1.39).sin();
    let ws = (-lat.tan() * decl.tan()).clamp(-1.0, 1.0).acos();
    let dr = 1.0 + 0.033 * (2.0 * std::f64::consts::PI * doy as f64 / 365.0).cos();
    // Radiación extraterrestre diaria (MJ m⁻² día⁻¹), FAO-56.
    let ra = (24.0 * 60.0 / std::f64::consts::PI)
        * 0.0820
        * dr
        * (ws * lat.sin() * decl.sin() + lat.cos() * decl.cos() * ws.sin());
    // Cielo despejado ~0.75·Ra; MJ/m²/día → W/m² medio diario (×1e6/86400).
    let g = (0.75 * ra).max(0.0) * 1.0e6 / 86_400.0;
    elevation.mapv(|z| if z.is_finite() { g } else { f64::NAN })
}

fn day_of_year(date: &str) -> u32 {
    let p: Vec<&str> = date.split('-').collect();
    let (m, d): (u32, u32) = (p[1].parse().unwrap(), p[2].parse().unwrap());
    const CUM: [u32; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    (CUM[m as usize - 1] + d).min(365)
}

/// Lee `balance.csv` (date, precip_mm, pet_mm, qobs_mm).
fn read_balance() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut rdr = csv::Reader::from_path(format!("{BASE}/balance.csv")).expect("balance.csv");
    let (mut p, mut pet, mut q) = (Vec::new(), Vec::new(), Vec::new());
    let num = |s: &str| -> f64 {
        let s = s.trim();
        if s.is_empty() || s.eq_ignore_ascii_case("na") {
            f64::NAN
        } else {
            s.parse().unwrap()
        }
    };
    for rec in rdr.records() {
        let rec = rec.unwrap();
        p.push(num(&rec[1]));
        pet.push(num(&rec[2]));
        q.push(num(&rec[3]));
    }
    (p, pet, q)
}

/// Calibra GR4J (DDS) con la mitad A y evalúa NSE en la mitad B (y viceversa).
fn split_sample(precip: &[f64], pet: &[f64], qobs: &[f64], label: &str) {
    let n = precip.len();
    let mid = n / 2;
    let opt = Optimizer::Dds(DdsConfig {
        max_iter: 3000,
        r: 0.2,
        seed: 42,
    });
    let folds = [("A→B", 0, mid, mid, n), ("B→A", mid, n, 0, mid)];
    print!("{label:<22}");
    for (_name, cs, ce, vs, ve) in folds {
        let cal = calibrate_gr4j(
            &precip[cs..ce],
            &pet[cs..ce],
            &qobs[cs..ce],
            WARMUP,
            Objective::Nse,
            &BOUNDS,
            &opt,
        )
        .expect("calibrate");
        let model = Gr4j::new(cal.params).expect("gr4j");
        let qsim = model.run(&precip[vs..ve], &pet[vs..ve]).expect("run");
        let val = metrics::nse(&qobs[vs + WARMUP..ve], &qsim[WARMUP..]).unwrap_or(f64::NAN);
        print!("  {val:>6.3}");
    }
    println!();
}

fn main() {
    eprintln!("Corriendo el modelo nival sobre las bandas...");
    let melt_input = snowmelt_input();
    let (precip, pet, qobs) = read_balance();
    assert_eq!(melt_input.len(), precip.len(), "series desalineadas");

    println!("\nAcople snowmelt-rs → rainflow (GR4J) — Río Choapa en Cuncumén");
    println!("Split-sample NSE de validación (DDS, 3000 evals, warm-up 365 d)\n");
    println!("{:<22}  {:>6}  {:>6}", "forzante de GR4J", "A→B", "B→A");
    println!("{}", "-".repeat(38));
    split_sample(&precip, &pet, &qobs, "precip cruda (sin nieve)");
    split_sample(&melt_input, &pet, &qobs, "aporte snowmelt (nieve)");
    println!("\nEl aporte de snowmelt mueve el agua de invierno (acumulación) a");
    println!("primavera (deshielo); GR4J cierra el balance suelo-escorrentía.");
}
