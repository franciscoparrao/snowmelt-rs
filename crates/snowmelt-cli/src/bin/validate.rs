//! `snowmelt-validate` — compara cobertura nival simulada contra observada.
//!
//! Recibe pares `SIM.asc:OBS.asc` (misma grilla). Una celda se considera
//! cubierta si su valor supera el umbral correspondiente; celdas nodata en
//! cualquiera de las dos grillas quedan fuera. Reporta matriz de confusión
//! y métricas por par, más el agregado.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::Parser;
use snowmelt_cli::asc;

#[derive(Parser, Debug)]
#[command(
    name = "snowmelt-validate",
    version,
    about = "Métricas de validación de cobertura nival simulada vs observada (ej. MODIS)"
)]
struct Cli {
    /// Pares simulado:observado en formato .asc (ej. cover_2019-08-01.asc:modis_2019-08-01.asc)
    #[arg(required = true, value_name = "SIM.asc:OBS.asc")]
    pairs: Vec<String>,

    /// Valor sobre el cual la grilla simulada cuenta como nieve
    #[arg(long, default_value_t = 0.5)]
    threshold_sim: f64,

    /// Valor sobre el cual la grilla observada cuenta como nieve
    #[arg(long, default_value_t = 0.5)]
    threshold_obs: f64,
}

#[derive(Debug, Default, Clone, Copy)]
struct Confusion {
    tp: u64,
    fp: u64,
    fn_: u64,
    tn: u64,
}

impl Confusion {
    fn valid(&self) -> u64 {
        self.tp + self.fp + self.fn_ + self.tn
    }
    fn accuracy(&self) -> f64 {
        (self.tp + self.tn) as f64 / self.valid() as f64
    }
    fn precision(&self) -> f64 {
        self.tp as f64 / (self.tp + self.fp) as f64
    }
    fn recall(&self) -> f64 {
        self.tp as f64 / (self.tp + self.fn_) as f64
    }
    fn f1(&self) -> f64 {
        let (p, r) = (self.precision(), self.recall());
        2.0 * p * r / (p + r)
    }
    /// Razón entre área nival simulada y observada (1 = sin sesgo).
    fn bias(&self) -> f64 {
        (self.tp + self.fp) as f64 / (self.tp + self.fn_) as f64
    }
    fn add(&mut self, other: &Confusion) {
        self.tp += other.tp;
        self.fp += other.fp;
        self.fn_ += other.fn_;
        self.tn += other.tn;
    }
}

fn compare(sim: &Path, obs: &Path, thr_sim: f64, thr_obs: f64) -> Result<Confusion> {
    let sim_grid = asc::read(sim).with_context(|| format!("leyendo {}", sim.display()))?;
    let obs_grid = asc::read(obs).with_context(|| format!("leyendo {}", obs.display()))?;
    if sim_grid.data.dim() != obs_grid.data.dim() {
        bail!(
            "grillas incompatibles: {} es {:?} y {} es {:?}",
            sim.display(),
            sim_grid.data.dim(),
            obs.display(),
            obs_grid.data.dim()
        );
    }
    let mut c = Confusion::default();
    for (&s, &o) in sim_grid.data.iter().zip(obs_grid.data.iter()) {
        if !s.is_finite() || !o.is_finite() {
            continue;
        }
        match (s >= thr_sim, o >= thr_obs) {
            (true, true) => c.tp += 1,
            (true, false) => c.fp += 1,
            (false, true) => c.fn_ += 1,
            (false, false) => c.tn += 1,
        }
    }
    if c.valid() == 0 {
        bail!("sin celdas válidas comunes entre las grillas");
    }
    Ok(c)
}

fn print_row(label: &str, c: &Confusion) {
    println!(
        "{label:<24} {:>8} {:>7} {:>7} {:>7} {:>7} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>7.3}",
        c.valid(),
        c.tp,
        c.fp,
        c.fn_,
        c.tn,
        c.accuracy(),
        c.precision(),
        c.recall(),
        c.f1(),
        c.bias()
    );
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    println!(
        "{:<24} {:>8} {:>7} {:>7} {:>7} {:>7} {:>9} {:>9} {:>9} {:>9} {:>7}",
        "par", "celdas", "TP", "FP", "FN", "TN", "accuracy", "precision", "recall", "F1", "bias"
    );
    let mut total = Confusion::default();
    for pair in &cli.pairs {
        let (sim, obs) = pair
            .split_once(':')
            .with_context(|| format!("par inválido `{pair}` (se espera SIM.asc:OBS.asc)"))?;
        let (sim, obs) = (PathBuf::from(sim), PathBuf::from(obs));
        let c = compare(&sim, &obs, cli.threshold_sim, cli.threshold_obs)?;
        let label = sim
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| pair.clone());
        print_row(&label, &c);
        total.add(&c);
    }
    if cli.pairs.len() > 1 {
        print_row("TOTAL", &total);
    }
    Ok(())
}
