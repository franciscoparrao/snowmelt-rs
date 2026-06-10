//! Degree-day model parameters.

use crate::error::{Result, SnowmeltError};

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
}

impl Default for DegreeDayParams {
    fn default() -> Self {
        Self {
            ddf: 4.0,
            t_melt: 0.0,
            t_snow: 0.0,
            t_rain: 2.0,
            lapse_rate: -0.0065,
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
