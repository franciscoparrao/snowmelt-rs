//! `snowmelt` — modelo grado-día distribuido sobre un DEM.
//!
//! Lee un DEM en ESRI ASCII Grid y una serie diaria de forzantes CSV,
//! simula acumulación/derretimiento de SWE y escribe la serie agregada
//! más la grilla final de SWE.

mod asc;
mod forcing_csv;

use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use snowmelt_core::{DegreeDayParams, Dem, Forcing, SnowModel};

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
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let grid = asc::read(&cli.dem)
        .with_context(|| format!("error leyendo el DEM {}", cli.dem.display()))?;
    let header = grid.header.clone();
    let dem = Dem::new(grid.data).context("DEM inválido")?;

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
    };
    let mut model = SnowModel::new(dem, params).context("parámetros inválidos")?;

    fs::create_dir_all(&cli.out_dir)
        .with_context(|| format!("no se pudo crear {}", cli.out_dir.display()))?;

    let mut series =
        String::from("date,snowfall_mm,rain_mm,melt_mm,runoff_mm,swe_mm,snow_cover_fraction\n");
    let mut total_melt = 0.0;
    let mut total_precip = 0.0;
    for rec in &records {
        let forcing = Forcing::Uniform {
            t_ref: rec.temp_c,
            z_ref,
            precip: rec.precip_mm,
        };
        let out = model
            .step(&forcing)
            .with_context(|| format!("fallo en el paso {}", rec.date))?;
        let s = model.summarize(&out);
        total_melt += s.mean_melt;
        total_precip += rec.precip_mm;
        let _ = writeln!(
            series,
            "{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.4}",
            rec.date,
            s.mean_snowfall,
            s.mean_rain,
            s.mean_melt,
            s.mean_runoff,
            s.mean_swe,
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
