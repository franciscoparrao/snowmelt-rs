//! Multi-year mass balance and equilibrium-line diagnostics.
//!
//! Integrates per-cell **accumulation** (snowfall) and **ablation** (melt +
//! sublimation) over an arbitrary run — one season, one hydrological year,
//! or many — so the model can be used for glacier-style surface mass
//! balance, not just single-season SWE. The net balance grid
//! (`accumulation − ablation`, mm w.e.) and its crossing of zero with
//! elevation (the **equilibrium line altitude**, ELA) summarise where the
//! surface gains or loses mass over the period.
//!
//! The accumulator is independent of the melt scheme: it consumes the
//! [`StepOutput`](crate::StepOutput) grids each step, so it works for
//! degree-day, ETI and energy-balance runs alike, including the sublimation
//! mass loss of the energy-balance mode.

use ndarray::{Array2, ArrayView2, Zip};

use crate::model::StepOutput;

/// Per-cell accumulation/ablation integrator over a run.
#[derive(Debug, Clone)]
pub struct MassBalance {
    accumulation: Array2<f64>,
    ablation: Array2<f64>,
}

impl MassBalance {
    /// A zeroed accumulator for a grid of shape `(rows, cols)`. Nodata cells
    /// become `NaN` on the first [`MassBalance::add`] (the step grids carry
    /// `NaN` there) and stay `NaN` thereafter.
    pub fn new(shape: (usize, usize)) -> Self {
        Self {
            accumulation: Array2::zeros(shape),
            ablation: Array2::zeros(shape),
        }
    }

    /// Adds one step's fluxes: snowfall to accumulation, melt and
    /// sublimation to ablation.
    pub fn add(&mut self, out: &StepOutput) {
        Zip::from(&mut self.accumulation)
            .and(&out.snowfall)
            .for_each(|acc, &s| *acc += s);
        Zip::from(&mut self.ablation)
            .and(&out.melt)
            .and(&out.sublimation)
            .for_each(|abl, &m, &sub| *abl += m + sub);
    }

    /// Cumulative accumulation (mm w.e.; `NaN` on nodata).
    pub fn accumulation(&self) -> ArrayView2<'_, f64> {
        self.accumulation.view()
    }

    /// Cumulative ablation = melt + sublimation (mm w.e.; `NaN` on nodata).
    pub fn ablation(&self) -> ArrayView2<'_, f64> {
        self.ablation.view()
    }

    /// Net surface mass balance per cell: `accumulation − ablation`
    /// (mm w.e.; positive = net gain). `NaN` on nodata.
    pub fn net(&self) -> Array2<f64> {
        &self.accumulation - &self.ablation
    }
}

/// Equilibrium line altitude (m): the elevation where the mean net balance
/// changes from negative (ablation-dominated, below) to positive
/// (accumulation-dominated, above).
///
/// Valid cells are binned into `n_bands` equal-elevation bands; the mean net
/// balance per band is interpolated to find the zero crossing. Returns
/// `None` if there is no sign change (the whole domain gains or loses mass),
/// if `n_bands < 2`, or if there are too few valid cells.
pub fn equilibrium_line_altitude(
    net: &ArrayView2<'_, f64>,
    elevation: &ArrayView2<'_, f64>,
    n_bands: usize,
) -> Option<f64> {
    if n_bands < 2 {
        return None;
    }
    let mut z_min = f64::INFINITY;
    let mut z_max = f64::NEG_INFINITY;
    for (&z, &b) in elevation.iter().zip(net.iter()) {
        if z.is_finite() && b.is_finite() {
            z_min = z_min.min(z);
            z_max = z_max.max(z);
        }
    }
    if !z_min.is_finite() || z_max <= z_min {
        return None;
    }
    let width = (z_max - z_min) / n_bands as f64;
    let mut sum = vec![0.0_f64; n_bands];
    let mut count = vec![0_usize; n_bands];
    for (&z, &b) in elevation.iter().zip(net.iter()) {
        if !(z.is_finite() && b.is_finite()) {
            continue;
        }
        let idx = (((z - z_min) / width) as usize).min(n_bands - 1);
        sum[idx] += b;
        count[idx] += 1;
    }
    // Non-empty band centres with their mean balance, low to high.
    let bands: Vec<(f64, f64)> = (0..n_bands)
        .filter(|&i| count[i] > 0)
        .map(|i| {
            let centre = z_min + (i as f64 + 0.5) * width;
            (centre, sum[i] / count[i] as f64)
        })
        .collect();
    // First adjacent pair crossing from negative to positive balance.
    for pair in bands.windows(2) {
        let (z0, b0) = pair[0];
        let (z1, b1) = pair[1];
        if b0 < 0.0 && b1 >= 0.0 {
            let frac = -b0 / (b1 - b0); // linear interpolation to balance = 0
            return Some(z0 + frac * (z1 - z0));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    fn step(snowfall: Array2<f64>, melt: Array2<f64>, sublimation: Array2<f64>) -> StepOutput {
        let rain = Array2::zeros(snowfall.dim());
        StepOutput {
            snowfall,
            rain,
            melt,
            sublimation,
        }
    }

    #[test]
    fn accumulates_and_nets_over_steps() {
        let mut mb = MassBalance::new((1, 2));
        mb.add(&step(
            array![[10.0, 5.0]],
            array![[2.0, 8.0]],
            array![[1.0, 1.0]],
        ));
        mb.add(&step(
            array![[4.0, 0.0]],
            array![[1.0, 3.0]],
            array![[0.0, 2.0]],
        ));
        // Cell 0: acc 14, abl (3 + 1) = 4 → net +10. Cell 1: acc 5, abl (11 + 3) = 14 → net −9.
        let net = mb.net();
        assert!(approx(net[[0, 0]], 10.0));
        assert!(approx(net[[0, 1]], -9.0));
        assert!(approx(mb.accumulation()[[0, 0]], 14.0));
        assert!(approx(mb.ablation()[[0, 1]], 14.0));
    }

    #[test]
    fn nodata_propagates() {
        let mut mb = MassBalance::new((1, 2));
        mb.add(&step(
            array![[10.0, f64::NAN]],
            array![[2.0, f64::NAN]],
            array![[0.0, f64::NAN]],
        ));
        let net = mb.net();
        assert!(approx(net[[0, 0]], 8.0));
        assert!(net[[0, 1]].is_nan());
    }

    #[test]
    fn ela_found_at_balance_crossing() {
        // Net balance increases with elevation: negative low, positive high.
        let elevation = array![[1000.0, 2000.0, 3000.0, 4000.0]];
        let net = array![[-20.0, -10.0, 10.0, 20.0]];
        let ela =
            equilibrium_line_altitude(&net.view(), &elevation.view(), 4).expect("crossing exists");
        // Crossing is between the 2000 m and 3000 m band centres, at the
        // midpoint by symmetry (−10 → +10).
        assert!(ela > 2000.0 && ela < 3000.0, "ELA {ela}");
        assert!((ela - 2500.0).abs() < 200.0, "ELA {ela}");
    }

    #[test]
    fn ela_none_when_all_positive() {
        let elevation = array![[2000.0, 3000.0, 4000.0]];
        let net = array![[5.0, 10.0, 15.0]];
        assert!(equilibrium_line_altitude(&net.view(), &elevation.view(), 3).is_none());
    }

    #[test]
    fn ela_none_for_degenerate_inputs() {
        let elevation = array![[2000.0, 3000.0]];
        let net = array![[-1.0, 1.0]];
        assert!(equilibrium_line_altitude(&net.view(), &elevation.view(), 1).is_none());
        let flat = array![[2000.0, 2000.0]];
        assert!(equilibrium_line_altitude(&net.view(), &flat.view(), 2).is_none());
    }
}
