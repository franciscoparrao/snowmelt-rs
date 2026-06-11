//! Degree-day model parameters.

use crate::error::{Result, SnowmeltError};

/// Snow albedo decay by age (Verseghy 1991-style exponential).
///
/// Albedo relaxes from `albedo_fresh` towards `albedo_min` with e-folding
/// time `tau_days`, and resets to fresh when a step's snowfall reaches
/// `refresh_swe_mm`: `α(t) = α_min + (α_fresh − α_min)·exp(−t/τ)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlbedoDecay {
    /// Albedo of fresh snow (0–1). Typical: 0.80–0.90.
    pub albedo_fresh: f64,
    /// Asymptotic albedo of old snow (0–1). Typical: 0.35–0.50.
    pub albedo_min: f64,
    /// E-folding decay time in days. Typical: 4–8.
    pub tau_days: f64,
    /// Snowfall per step (mm w.e.) that resets the surface to fresh snow.
    pub refresh_swe_mm: f64,
}

impl Default for AlbedoDecay {
    fn default() -> Self {
        Self {
            albedo_fresh: 0.85,
            albedo_min: 0.4,
            tau_days: 6.0,
            refresh_swe_mm: 1.0,
        }
    }
}

impl AlbedoDecay {
    /// Albedo for snow of age `age_days`. `NaN` for `NaN` input.
    pub fn albedo(&self, age_days: f64) -> f64 {
        self.albedo_min + (self.albedo_fresh - self.albedo_min) * (-age_days / self.tau_days).exp()
    }
}

/// Parameters of the degree-day (temperature-index) snow model.
///
/// Units: temperatures in °C, degree-day factor in mm w.e. °C⁻¹ day⁻¹,
/// lapse rate in °C m⁻¹ (negative for the usual decrease with elevation).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DegreeDayParams {
    /// Degree-day factor (mm w.e. °C⁻¹ day⁻¹). Typical snow values: 2–6.
    pub ddf: f64,
    /// Temperature above which melt occurs (°C).
    pub t_melt: f64,
    /// At or below this temperature all precipitation falls as snow (°C).
    pub t_snow: f64,
    /// At or above this temperature all precipitation falls as rain (°C).
    /// Between `t_snow` and `t_rain` the snow fraction decreases linearly.
    pub t_rain: f64,
    /// Vertical temperature lapse rate (°C m⁻¹), applied when extrapolating
    /// a reference temperature over the DEM. Typically `-0.0065`.
    pub lapse_rate: f64,
    /// Shortwave radiation factor (mm day⁻¹ (W m⁻²)⁻¹) of the enhanced
    /// temperature-index model (Pellicciotti et al. 2005). With `srf > 0`
    /// melt becomes `ddf·(T − t_melt) + srf·(1 − albedo)·G` for `T > t_melt`,
    /// and a radiation grid `G` (W m⁻², daily mean) must be supplied via
    /// [`SnowModel::step_radiation`](crate::SnowModel::step_radiation).
    /// `0` (default) disables the radiative term (pure degree-day).
    /// Typical daily value: ~0.2 (0.0094 mm h⁻¹ (W m⁻²)⁻¹ · 24).
    pub srf: f64,
    /// Snow albedo (0–1) used by the radiative melt term when
    /// [`albedo_decay`](Self::albedo_decay) is `None`. Typical: 0.4–0.8.
    pub albedo: f64,
    /// Optional age-dependent albedo. When set, the per-cell albedo decays
    /// with days since the last significant snowfall and `albedo` is
    /// ignored.
    pub albedo_decay: Option<AlbedoDecay>,
    /// Orographic precipitation gradient (m⁻¹), applied to uniform forcings:
    /// `p(z) = p_ref · (1 + precip_gradient·(z − z_ref))`, clamped to ≥ 0.
    /// `0` (default) means uniform precipitation. Typical: 0.0002–0.001.
    pub precip_gradient: f64,
}

impl Default for DegreeDayParams {
    fn default() -> Self {
        Self {
            ddf: 4.0,
            t_melt: 0.0,
            t_snow: 0.0,
            t_rain: 2.0,
            lapse_rate: -0.0065,
            srf: 0.0,
            albedo: 0.6,
            albedo_decay: None,
            precip_gradient: 0.0,
        }
    }
}

impl DegreeDayParams {
    /// Checks that all parameters are finite and mutually consistent.
    ///
    /// # Errors
    /// Returns [`SnowmeltError::InvalidParameter`] for a non-finite value,
    /// a negative `ddf`, or `t_snow > t_rain`.
    pub fn validate(&self) -> Result<()> {
        let finite = [
            ("ddf", self.ddf),
            ("t_melt", self.t_melt),
            ("t_snow", self.t_snow),
            ("t_rain", self.t_rain),
            ("lapse_rate", self.lapse_rate),
            ("srf", self.srf),
            ("albedo", self.albedo),
            ("precip_gradient", self.precip_gradient),
        ];
        for (name, value) in finite {
            if !value.is_finite() {
                return Err(SnowmeltError::InvalidParameter {
                    name,
                    reason: format!("must be finite, got {value}"),
                });
            }
        }
        if self.ddf < 0.0 {
            return Err(SnowmeltError::InvalidParameter {
                name: "ddf",
                reason: format!("must be >= 0, got {}", self.ddf),
            });
        }
        if self.srf < 0.0 {
            return Err(SnowmeltError::InvalidParameter {
                name: "srf",
                reason: format!("must be >= 0, got {}", self.srf),
            });
        }
        if !(0.0..=1.0).contains(&self.albedo) {
            return Err(SnowmeltError::InvalidParameter {
                name: "albedo",
                reason: format!("must be in [0, 1], got {}", self.albedo),
            });
        }
        if let Some(decay) = &self.albedo_decay {
            for (name, value) in [
                ("albedo_fresh", decay.albedo_fresh),
                ("albedo_min", decay.albedo_min),
            ] {
                if !(0.0..=1.0).contains(&value) {
                    return Err(SnowmeltError::InvalidParameter {
                        name,
                        reason: format!("must be in [0, 1], got {value}"),
                    });
                }
            }
            if decay.albedo_min > decay.albedo_fresh {
                return Err(SnowmeltError::InvalidParameter {
                    name: "albedo_min",
                    reason: format!(
                        "albedo_min ({}) must be <= albedo_fresh ({})",
                        decay.albedo_min, decay.albedo_fresh
                    ),
                });
            }
            if !decay.tau_days.is_finite() || decay.tau_days <= 0.0 {
                return Err(SnowmeltError::InvalidParameter {
                    name: "tau_days",
                    reason: format!("must be finite and > 0, got {}", decay.tau_days),
                });
            }
            if !decay.refresh_swe_mm.is_finite() || decay.refresh_swe_mm < 0.0 {
                return Err(SnowmeltError::InvalidParameter {
                    name: "refresh_swe_mm",
                    reason: format!("must be finite and >= 0, got {}", decay.refresh_swe_mm),
                });
            }
        }
        if self.t_snow > self.t_rain {
            return Err(SnowmeltError::InvalidParameter {
                name: "t_snow",
                reason: format!(
                    "t_snow ({}) must be <= t_rain ({})",
                    self.t_snow, self.t_rain
                ),
            });
        }
        Ok(())
    }

    /// Fraction of precipitation falling as snow at air temperature `t_c`.
    ///
    /// 1 at or below [`t_snow`](Self::t_snow), 0 at or above
    /// [`t_rain`](Self::t_rain), linear in between. Returns `NaN` for `NaN`
    /// input.
    pub fn snow_fraction(&self, t_c: f64) -> f64 {
        if t_c <= self.t_snow {
            1.0
        } else if t_c >= self.t_rain {
            0.0
        } else {
            (self.t_rain - t_c) / (self.t_rain - self.t_snow)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_params_are_valid() {
        DegreeDayParams::default().validate().unwrap();
    }

    #[test]
    fn snow_fraction_edges_and_midpoint() {
        let p = DegreeDayParams::default(); // t_snow = 0, t_rain = 2
        assert_eq!(p.snow_fraction(-5.0), 1.0);
        assert_eq!(p.snow_fraction(0.0), 1.0);
        assert_eq!(p.snow_fraction(1.0), 0.5);
        assert_eq!(p.snow_fraction(2.0), 0.0);
        assert_eq!(p.snow_fraction(10.0), 0.0);
        assert!(p.snow_fraction(f64::NAN).is_nan());
    }

    #[test]
    fn snow_fraction_is_a_step_when_thresholds_coincide() {
        let p = DegreeDayParams {
            t_snow: 1.0,
            t_rain: 1.0,
            ..DegreeDayParams::default()
        };
        p.validate().unwrap();
        assert_eq!(p.snow_fraction(1.0), 1.0);
        assert_eq!(p.snow_fraction(1.0001), 0.0);
    }

    #[test]
    fn rejects_negative_ddf_and_crossed_thresholds() {
        let p = DegreeDayParams {
            ddf: -1.0,
            ..DegreeDayParams::default()
        };
        assert!(p.validate().is_err());

        let p = DegreeDayParams {
            t_snow: 3.0,
            t_rain: 1.0,
            ..DegreeDayParams::default()
        };
        assert!(p.validate().is_err());
    }

    #[test]
    fn rejects_non_finite_values() {
        let p = DegreeDayParams {
            lapse_rate: f64::NAN,
            ..DegreeDayParams::default()
        };
        assert!(p.validate().is_err());
    }
}
