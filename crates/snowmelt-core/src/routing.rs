//! Linear-reservoir routing of catchment water input to a hydrograph.
//!
//! The snow model produces a per-step water input (rain + melt) averaged
//! over the catchment. Converting that into a streamflow series at the
//! outlet needs a transfer function for catchment storage and travel
//! time. A single linear reservoir (Nash, one box) is the simplest
//! defensible choice: outflow proportional to storage, exponential
//! recession with time constant `k`.
//!
//! Mass is conserved exactly each step: `outflow = input + (S_old − S_new)`.
//!
//! This is intentionally minimal — it gives the snowmelt hydrograph its
//! timing and recession so it can be compared with observed discharge or
//! handed to a full rainfall–runoff model ([rainflow]) as a forcing. It is
//! **not** a rainfall–runoff balance: it neither removes evapotranspiration
//! nor adds baseflow from non-snow processes.

use crate::error::{Result, SnowmeltError};

/// A single linear reservoir: `dS/dt = I − S/k`, `Q = S/k`.
#[derive(Debug, Clone)]
pub struct LinearReservoir {
    storage_mm: f64,
    k_days: f64,
}

impl LinearReservoir {
    /// Creates a reservoir with recession time constant `k_days` and zero
    /// initial storage.
    ///
    /// # Errors
    /// Returns [`SnowmeltError::InvalidParameter`] if `k_days` is not finite
    /// and positive.
    pub fn new(k_days: f64) -> Result<Self> {
        Self::with_storage(k_days, 0.0)
    }

    /// Creates a reservoir with a prescribed initial storage (mm).
    ///
    /// # Errors
    /// As [`Self::new`], plus a non-finite or negative `storage_mm`.
    pub fn with_storage(k_days: f64, storage_mm: f64) -> Result<Self> {
        if !k_days.is_finite() || k_days <= 0.0 {
            return Err(SnowmeltError::InvalidParameter {
                name: "k_days",
                reason: format!("must be finite and > 0, got {k_days}"),
            });
        }
        if !storage_mm.is_finite() || storage_mm < 0.0 {
            return Err(SnowmeltError::InvalidParameter {
                name: "storage_mm",
                reason: format!("must be finite and >= 0, got {storage_mm}"),
            });
        }
        Ok(Self { storage_mm, k_days })
    }

    /// Current storage (mm).
    pub fn storage(&self) -> f64 {
        self.storage_mm
    }

    /// Routes one step: adds `input_mm` (water input over `dt_days`) and
    /// returns the outflow (mm) released during the step.
    ///
    /// Uses the exact solution of the linear reservoir over the step with
    /// the input applied as a constant rate, so mass is conserved:
    /// `outflow = input_mm + (S_old − S_new)`.
    ///
    /// Non-finite or non-positive `dt_days`, or a non-finite `input_mm`,
    /// yield `NaN` and leave storage `NaN` (propagating nodata).
    pub fn step(&mut self, input_mm: f64, dt_days: f64) -> f64 {
        if !input_mm.is_finite() || !dt_days.is_finite() || dt_days <= 0.0 {
            self.storage_mm = f64::NAN;
            return f64::NAN;
        }
        let decay = (-dt_days / self.k_days).exp();
        let rate = input_mm / dt_days; // mm/day inflow over the step
        let s_old = self.storage_mm;
        let s_new = s_old * decay + rate * self.k_days * (1.0 - decay);
        self.storage_mm = s_new;
        input_mm + s_old - s_new
    }

    /// Routes a whole input series (mm per step), returning the outflow
    /// series (mm per step). Convenience over repeated [`Self::step`].
    pub fn route(&mut self, input: &[f64], dt_days: f64) -> Vec<f64> {
        input.iter().map(|&i| self.step(i, dt_days)).collect()
    }
}

/// Converts a depth over a catchment (mm) to a mean discharge (m³ s⁻¹)
/// over `dt_days`, given the catchment `area_km2`.
pub fn depth_to_discharge(depth_mm: f64, area_km2: f64, dt_days: f64) -> f64 {
    // mm over km² = 1e-3 m × 1e6 m² = 1e3 m³; per (dt_days × 86400 s).
    depth_mm * area_km2 * 1.0e3 / (dt_days * 86_400.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(LinearReservoir::new(0.0).is_err());
        assert!(LinearReservoir::new(-1.0).is_err());
        assert!(LinearReservoir::new(f64::NAN).is_err());
        assert!(LinearReservoir::with_storage(5.0, -1.0).is_err());
    }

    #[test]
    fn pure_recession_is_exponential() {
        let mut r = LinearReservoir::with_storage(10.0, 100.0).unwrap();
        let out = r.step(0.0, 1.0);
        // S decays by e^{-1/10}; outflow is the released difference.
        let expected_s = 100.0 * (-0.1_f64).exp();
        assert!(approx(r.storage(), expected_s));
        assert!(approx(out, 100.0 - expected_s));
    }

    #[test]
    fn conserves_mass_over_a_series() {
        let mut r = LinearReservoir::new(8.0).unwrap();
        let input = [12.0, 4.0, 0.0, 30.0, 5.0, 0.0, 0.0, 1.0];
        let out = r.route(&input, 1.0);
        let total_in: f64 = input.iter().sum();
        let total_out: f64 = out.iter().sum();
        // input == outflow + storage still held back.
        assert!(approx(total_in, total_out + r.storage()));
        assert!(out.iter().all(|&q| q >= 0.0));
    }

    #[test]
    fn steady_state_outflow_approaches_input() {
        let mut r = LinearReservoir::new(3.0).unwrap();
        let mut out = 0.0;
        for _ in 0..200 {
            out = r.step(5.0, 1.0);
        }
        // Constant 5 mm/day input → outflow converges to 5 mm/day.
        assert!((out - 5.0).abs() < 1e-6, "out = {out}");
    }

    #[test]
    fn small_k_is_near_passthrough() {
        // With k ≪ dt the reservoir empties almost fully each step, holding
        // back only ≈ rate·k = 10·0.01 = 0.1 mm.
        let mut r = LinearReservoir::new(0.01).unwrap();
        let out = r.step(10.0, 1.0);
        assert!(out > 9.8, "out = {out}");
    }

    #[test]
    fn nan_input_propagates() {
        let mut r = LinearReservoir::new(5.0).unwrap();
        assert!(r.step(f64::NAN, 1.0).is_nan());
        assert!(r.storage().is_nan());
    }

    #[test]
    fn depth_to_discharge_units() {
        // 1 mm/day over 86.4 km² = 1 m³/s.
        assert!((depth_to_discharge(1.0, 86.4, 1.0) - 1.0).abs() < 1e-9);
    }
}
