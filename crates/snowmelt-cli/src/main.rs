//! `snowmelt` — modelo grado-día distribuido sobre un DEM.
//!
//! Lee un DEM en ESRI ASCII Grid y una serie diaria de forzantes CSV,
//! simula acumulación/derretimiento de SWE y escribe la serie agregada
//! más la grilla final de SWE.

mod asc;
mod forcing_csv;
mod solar;

use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use ndarray::Array2;

use anyhow::{Context, Result};
use clap::Parser;
use snowmelt_core::{AlbedoDecay, DegreeDayParams, Dem, Forcing, SnowModel};
use surtgis_algorithms::terrain::HorizonParams;

#[derive(Parser, Debug)]
#[command(
    name = "snowmelt",
    version,
    about = "Modelo grado-día distribuido de derretimiento nival sobre un DEM"
)]
struct Cli {
    /// DEM en formato ESRI ASCII Grid (.asc); celdas NODATA quedan fuera del balance
    #[arg(long)]
    dem: PathBuf,

    /// Serie diaria de forzantes CSV: date,temp_c,precip_mm (header opcional)
    #[arg(long)]
    forcing: PathBuf,

    /// Directorio de salida (se crea si no existe)
    #[arg(long, default_value = "out")]
    out_dir: PathBuf,

    /// Factor grado-día [mm °C⁻¹ día⁻¹]
    #[arg(long, default_value_t = 4.0)]
    ddf: f64,

    /// Temperatura umbral de fusión [°C]
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    t_melt: f64,

    /// Bajo esta temperatura toda la precipitación es nieve [°C]
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    t_snow: f64,

    /// Sobre esta temperatura toda la precipitación es lluvia [°C]
    #[arg(long, default_value_t = 2.0, allow_hyphen_values = true)]
    t_rain: f64,

    /// Gradiente térmico vertical [°C m⁻¹]
    #[arg(long, default_value_t = -0.0065, allow_hyphen_values = true)]
    lapse_rate: f64,

    /// Elevación a la que se midió la temperatura forzante [m].
    /// Si se omite, se usa la elevación media del DEM.
    #[arg(long, allow_hyphen_values = true)]
    z_ref: Option<f64>,

    /// Factor de radiación de onda corta del modelo ETI
    /// [mm día⁻¹ (W m⁻²)⁻¹]. 0 = grado-día puro; típico ~0.2.
    /// Si es > 0 se calcula radiación potencial desde el DEM (requiere --latitude).
    #[arg(long, default_value_t = 0.0)]
    srf: f64,

    /// Albedo de la nieve (0-1) para el término radiativo
    #[arg(long, default_value_t = 0.6)]
    albedo: f64,

    /// Gradiente orográfico de precipitación [m⁻¹]:
    /// p(z) = p_ref·(1 + grad·(z − z_ref)), acotado a ≥ 0
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    precip_gradient: f64,

    /// Latitud del dominio [°, negativa en el hemisferio sur].
    /// Requerida cuando --srf > 0.
    #[arg(long, allow_hyphen_values = true)]
    latitude: Option<f64>,

    /// Transmitancia atmosférica de cielo despejado (0-1)
    #[arg(long, default_value_t = 0.7)]
    transmittance: f64,

    /// Factor de turbidez de Linke (2 muy claro, 3 claro, 4 brumoso).
    /// Si se especifica, reemplaza a --transmittance (modelo Kasten 1996).
    #[arg(long)]
    linke_turbidity: Option<f64>,

    /// Sombreado por horizonte topográfico (cast shadows). Precalcula
    /// ángulos de horizonte; memoria ≈ 8·direcciones·celdas bytes.
    #[arg(long, default_value_t = false)]
    horizon_shading: bool,

    /// Radio de búsqueda del horizonte [celdas]
    #[arg(long, default_value_t = 100)]
    horizon_radius: usize,

    /// Direcciones acimutales para el horizonte
    #[arg(long, default_value_t = 36)]
    horizon_directions: usize,

    /// Activa albedo dinámico con decaimiento exponencial por edad de la
    /// nieve: α(t) = α_min + (α_fresh − α_min)·exp(−t/τ). Valor = τ [días].
    #[arg(long)]
    albedo_tau: Option<f64>,

    /// Albedo de nieve fresca (modo decaimiento)
    #[arg(long, default_value_t = 0.85)]
    albedo_fresh: f64,

    /// Albedo asintótico de nieve vieja (modo decaimiento)
    #[arg(long, default_value_t = 0.4)]
    albedo_min: f64,

    /// Nevada por paso [mm] que reinicia el albedo a fresco
    #[arg(long, default_value_t = 1.0)]
    albedo_refresh: f64,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let grid = asc::read(&cli.dem)
        .with_context(|| format!("error leyendo el DEM {}", cli.dem.display()))?;
    let header = grid.header.clone();
    let elevation = grid.data.clone();
    let dem = Dem::new(grid.data).context("DEM inválido")?;

    let terrain = if cli.srf > 0.0 {
        if cli.latitude.is_none() {
            anyhow::bail!("--srf > 0 requiere --latitude para calcular la radiación potencial");
        }
        let horizon = cli.horizon_shading.then_some(HorizonParams {
            radius: cli.horizon_radius,
            directions: cli.horizon_directions,
        });
        Some(
            solar::Terrain::from_dem(&elevation, &header, horizon)
                .context("derivando slope/aspect")?,
        )
    } else {
        None
    };

    let records = forcing_csv::read(&cli.forcing)
        .with_context(|| format!("error leyendo forzantes {}", cli.forcing.display()))?;

    let z_ref = match cli.z_ref {
        Some(z) => z,
        None => {
            let z = dem.mean_elevation();
            eprintln!("z_ref no especificado: usando elevación media del DEM ({z:.1} m)");
            z
        }
    };

    let params = DegreeDayParams {
        ddf: cli.ddf,
        t_melt: cli.t_melt,
        t_snow: cli.t_snow,
        t_rain: cli.t_rain,
        lapse_rate: cli.lapse_rate,
        srf: cli.srf,
        albedo: cli.albedo,
        albedo_decay: cli.albedo_tau.map(|tau| AlbedoDecay {
            albedo_fresh: cli.albedo_fresh,
            albedo_min: cli.albedo_min,
            tau_days: tau,
            refresh_swe_mm: cli.albedo_refresh,
        }),
        precip_gradient: cli.precip_gradient,
    };
    let mut model = SnowModel::new(dem, params).context("parámetros inválidos")?;

    fs::create_dir_all(&cli.out_dir)
        .with_context(|| format!("no se pudo crear {}", cli.out_dir.display()))?;

    let mut series = String::from(
        "date,snowfall_mm,rain_mm,melt_mm,runoff_mm,swe_mm,albedo,snow_cover_fraction\n",
    );
    let mut total_melt = 0.0;
    let mut total_precip = 0.0;
    // Cache de radiación potencial por día del año (se repite entre años).
    let mut rad_cache: HashMap<u32, Array2<f64>> = HashMap::new();
    for rec in &records {
        let forcing = Forcing::Uniform {
            t_ref: rec.temp_c,
            z_ref,
            precip: rec.precip_mm,
        };
        let radiation = match &terrain {
            Some(terrain) => {
                let doy =
                    solar::day_of_year(&rec.date).with_context(|| format!("paso {}", rec.date))?;
                let rad = match rad_cache.entry(doy) {
                    std::collections::hash_map::Entry::Occupied(e) => e.into_mut(),
                    std::collections::hash_map::Entry::Vacant(e) => {
                        let rad = terrain
                            .potential_radiation(
                                &elevation,
                                doy,
                                cli.latitude.expect("validado arriba"),
                                cli.transmittance,
                                cli.linke_turbidity,
                                cli.albedo,
                            )
                            .with_context(|| format!("radiación para el día {doy}"))?;
                        e.insert(rad)
                    }
                };
                Some(&*rad)
            }
            None => None,
        };
        let out = model
            .step_radiation(&forcing, radiation.map(|r| r.view()), 1.0)
            .with_context(|| format!("fallo en el paso {}", rec.date))?;
        let s = model.summarize(&out);
        total_melt += s.mean_melt;
        total_precip += rec.precip_mm;
        let _ = writeln!(
            series,
            "{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.4},{:.4}",
            rec.date,
            s.mean_snowfall,
            s.mean_rain,
            s.mean_melt,
            s.mean_runoff,
            s.mean_swe,
            s.mean_albedo,
            s.snow_cover_fraction
        );
    }

    let series_path = cli.out_dir.join("series.csv");
    fs::write(&series_path, series)
        .with_context(|| format!("no se pudo escribir {}", series_path.display()))?;

    let swe_path = cli.out_dir.join("swe_final.asc");
    asc::write(&swe_path, &header, model.swe())
        .with_context(|| format!("no se pudo escribir {}", swe_path.display()))?;

    let final_swe = {
        let (sum, n) = model
            .swe()
            .iter()
            .filter(|v| v.is_finite())
            .fold((0.0, 0usize), |(s, n), &v| (s + v, n + 1));
        if n == 0 { f64::NAN } else { sum / n as f64 }
    };
    println!("Simulación completada: {} pasos diarios", records.len());
    println!("  celdas válidas      : {}", model.dem().valid_cells());
    println!("  precipitación total : {total_precip:.1} mm");
    println!("  derretimiento medio : {total_melt:.1} mm");
    println!("  SWE medio final     : {final_swe:.1} mm");
    println!("  serie agregada      : {}", series_path.display());
    println!("  SWE final (grilla)  : {}", swe_path.display());
    Ok(())
}
