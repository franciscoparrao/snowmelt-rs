//! Single-layer snowpack energy balance.
//!
//! Net energy available to the pack, per cell and step (all W m⁻²):
//!
//! ```text
//! Q = (1 − α)·G  +  LW_in − LW_out  +  Q_H  +  Q_E  +  Q_G
//! ```
//!
//! - `(1 − α)·G`: net shortwave from the potential-radiation forcing and
//!   the model's (constant or age-decayed) albedo.
//! - `LW_in = ε_a·σ·T_a⁴` with clear-sky atmospheric emissivity after
//!   Brutsaert (1975); `LW_out = ε_s·σ·T_s⁴` with the snow surface at
//!   `T_s = min(T_a, 0 °C)`.
//! - `Q_H`, `Q_E`: bulk-aerodynamic sensible and latent fluxes with a
//!   single exchange coefficient, parameterised wind speed and relative
//!   humidity, and air density from the elevation-derived pressure.
//! - `Q_G`: constant ground heat flux.
//!
//! Negative balances build **cold content** (J m⁻², capped by the pack's
//! heat capacity over [`EnergyBalanceParams::t_cold_max`] kelvin); positive
//! balances first pay it off, then melt at `L_f = 334 kJ kg⁻¹`.
//!
//! Simplifications (documented for v0.4): no rain-on-snow heat, no mass
//! loss by sublimation, clear-sky longwave (consistent with the clear-sky
//! shortwave forcing).

use crate::error::{Result, SnowmeltError};

/// Stefan–Boltzmann constant (W m⁻² K⁻⁴).
const SIGMA: f64 = 5.670374419e-8;
/// Specific heat of air (J kg⁻¹ K⁻¹).
const CP_AIR: f64 = 1005.0;
/// Specific heat of ice (J kg⁻¹ K⁻¹).
const C_ICE: f64 = 2100.0;
/// Latent heat of fusion (J kg⁻¹).
const L_FUSION: f64 = 334_000.0;
/// Latent heat of sublimation (J kg⁻¹).
const L_SUBLIMATION: f64 = 2.834e6;
/// Gas constant of dry air (J kg⁻¹ K⁻¹).
const R_DRY: f64 = 287.05;
/// Seconds per day.
const DAY_S: f64 = 86_400.0;

/// Parameters of the energy-balance melt mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyBalanceParams {
    /// Wind speed (m s⁻¹). Typical: 1–4.
    pub wind_speed: f64,
    /// Relative humidity (0–1). Typical: 0.4–0.8.
    pub rel_humidity: f64,
    /// Snow surface emissivity (0–1]. Typical: 0.97–0.99.
    pub snow_emissivity: f64,
    /// Bulk exchange coefficient for sensible and latent heat (–).
    /// Typical: 0.001–0.003.
    pub exchange_coeff: f64,
    /// Ground heat flux into the pack (W m⁻²). Typical: 0–2.
    pub ground_heat: f64,
    /// Maximum pack cooling below 0 °C (K) used to cap the cold content
    /// at `c_ice·SWE·t_cold_max`. Typical: 5–15.
    pub t_cold_max: f64,
}

impl Default for EnergyBalanceParams {
    fn default() -> Self {
        Self {
            wind_speed: 2.0,
            rel_humidity: 0.6,
            snow_emissivity: 0.98,
            exchange_coeff: 0.0015,
            ground_heat: 1.0,
            t_cold_max: 10.0,
        }
    }
}

impl EnergyBalanceParams {
    /// Checks that all parameters are finite and within their domains.
    ///
    /// # Errors
    /// Returns [`SnowmeltError::InvalidParameter`] on the first violation.
    pub fn validate(&self) -> Result<()> {
        let checks: [(&'static str, f64, bool); 6] = [
            ("wind_speed", self.wind_speed, self.wind_speed >= 0.0),
            (
                "rel_humidity",
                self.rel_humidity,
                (0.0..=1.0).contains(&self.rel_humidity),
            ),
            (
                "snow_emissivity",
                self.snow_emissivity,
                self.snow_emissivity > 0.0 && self.snow_emissivity <= 1.0,
            ),
            (
                "exchange_coeff",
                self.exchange_coeff,
                self.exchange_coeff > 0.0,
            ),
            ("ground_heat", self.ground_heat, true),
            ("t_cold_max", self.t_cold_max, self.t_cold_max >= 0.0),
        ];
        for (name, value, in_domain) in checks {
            if !value.is_finite() || !in_domain {
                return Err(SnowmeltError::InvalidParameter {
                    name,
                    reason: format!("out of domain: {value}"),
                });
            }
        }
        Ok(())
    }
}

/// Saturation vapour pressure over ice/water (hPa), Magnus form.
fn e_sat_hpa(t_c: f64) -> f64 {
    6.112 * (17.62 * t_c / (243.12 + t_c)).exp()
}

/// Air pressure (Pa) at elevation `z` (m), standard atmosphere.
pub(crate) fn air_pressure_pa(z_m: f64) -> f64 {
    101_325.0 * (1.0 - 2.25577e-5 * z_m).powf(5.2559)
}

/// Net energy flux into the snowpack (W m⁻²) for one cell.
///
/// `sw_in` is the incoming shortwave (W m⁻², daily mean), `albedo` the
/// current surface albedo, `pressure_pa` the cell's air pressure.
pub(crate) fn net_energy(
    p: &EnergyBalanceParams,
    t_air_c: f64,
    sw_in: f64,
    albedo: f64,
    pressure_pa: f64,
) -> f64 {
    let t_air_k = t_air_c + 273.15;
    // Snow surface cannot exceed the melting point.
    let t_surf_c = t_air_c.min(0.0);
    let t_surf_k = t_surf_c + 273.15;

    // Shortwave.
    let q_sw = (1.0 - albedo) * sw_in;

    // Longwave: Brutsaert (1975) clear-sky emissivity.
    let e_air_hpa = p.rel_humidity * e_sat_hpa(t_air_c);
    let emiss_air = 1.24 * (e_air_hpa / t_air_k).powf(1.0 / 7.0);
    let q_lw = emiss_air * SIGMA * t_air_k.powi(4) - p.snow_emissivity * SIGMA * t_surf_k.powi(4);

    // Turbulent fluxes (bulk aerodynamic).
    let rho_air = pressure_pa / (R_DRY * t_air_k);
    let q_h = rho_air * CP_AIR * p.exchange_coeff * p.wind_speed * (t_air_c - t_surf_c);
    let q_air = 0.622 * (e_air_hpa * 100.0) / pressure_pa;
    let q_surf = 0.622 * (e_sat_hpa(t_surf_c) * 100.0) / pressure_pa;
    let q_e = rho_air * L_SUBLIMATION * p.exchange_coeff * p.wind_speed * (q_air - q_surf);

    q_sw + q_lw + q_h + q_e + p.ground_heat
}

/// Applies one step's energy to a cell's cold content and returns the
/// potential melt (mm w.e.; the caller caps it by the available SWE).
///
/// Negative energy cools the pack (cold content grows, capped by
/// `c_ice·swe·t_cold_max`); positive energy first cancels the cold
/// content, the remainder melts.
pub(crate) fn apply_energy(
    p: &EnergyBalanceParams,
    q_wm2: f64,
    dt_days: f64,
    swe_mm: f64,
    cold_content_jm2: &mut f64,
) -> f64 {
    let energy = q_wm2 * dt_days * DAY_S; // J/m²
    if swe_mm <= 0.0 {
        *cold_content_jm2 = 0.0;
        return 0.0;
    }
    if energy < 0.0 {
        let cap = C_ICE * swe_mm * p.t_cold_max; // swe en kg/m² == mm
        *cold_content_jm2 = (*cold_content_jm2 - energy).min(cap);
        return 0.0;
    }
    let after_cold = energy - *cold_content_jm2;
    if after_cold <= 0.0 {
        *cold_content_jm2 -= energy;
        return 0.0;
    }
    *cold_content_jm2 = 0.0;
    after_cold / L_FUSION // kg/m² == mm
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_params_validate() {
        EnergyBalanceParams::default().validate().unwrap();
    }

    #[test]
    fn rejects_out_of_domain_params() {
        for bad in [
            EnergyBalanceParams {
                rel_humidity: 1.5,
                ..Default::default()
            },
            EnergyBalanceParams {
                exchange_coeff: 0.0,
                ..Default::default()
            },
            EnergyBalanceParams {
                wind_speed: f64::NAN,
                ..Default::default()
            },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
    }

    #[test]
    fn pressure_decreases_with_elevation() {
        assert!((air_pressure_pa(0.0) - 101_325.0).abs() < 1.0);
        let p3000 = air_pressure_pa(3000.0);
        assert!(p3000 > 65_000.0 && p3000 < 72_000.0, "{p3000}");
    }

    #[test]
    fn warm_sunny_day_yields_positive_energy() {
        let p = EnergyBalanceParams::default();
        let q = net_energy(&p, 5.0, 250.0, 0.6, air_pressure_pa(2000.0));
        // SW neto = 100 W/m²; LW pierde, turbulentos aportan algo.
        assert!(q > 30.0 && q < 250.0, "q = {q}");
    }

    #[test]
    fn cold_clear_night_yields_negative_energy() {
        let p = EnergyBalanceParams::default();
        let q = net_energy(&p, -15.0, 0.0, 0.8, air_pressure_pa(3000.0));
        assert!(q < 0.0, "q = {q}");
    }

    #[test]
    fn cold_content_buffers_melt() {
        let p = EnergyBalanceParams::default();
        let mut cc = 0.0;
        // Día frío: -50 W/m² sobre 100 mm de SWE → acumula frío, no derrite.
        let melt = apply_energy(&p, -50.0, 1.0, 100.0, &mut cc);
        assert_eq!(melt, 0.0);
        assert!(cc > 0.0);
        let cc_after_cold = cc;
        // Día levemente positivo: +20 W/m² = 1.728 MJ < cc (4.32 MJ) → aún 0.
        let melt = apply_energy(&p, 20.0, 1.0, 100.0, &mut cc);
        assert_eq!(melt, 0.0);
        assert!(cc < cc_after_cold);
        // Día muy cálido: paga el resto del frío y derrite.
        let melt = apply_energy(&p, 150.0, 1.0, 100.0, &mut cc);
        assert!(melt > 0.0);
        assert_eq!(cc, 0.0);
    }

    #[test]
    fn cold_content_is_capped_and_cleared_without_snow() {
        let p = EnergyBalanceParams::default();
        let mut cc = 0.0;
        for _ in 0..100 {
            apply_energy(&p, -200.0, 1.0, 10.0, &mut cc);
        }
        let cap = C_ICE * 10.0 * p.t_cold_max;
        assert!(cc <= cap + 1e-6, "cc {cc} > cap {cap}");
        apply_energy(&p, -200.0, 1.0, 0.0, &mut cc);
        assert_eq!(cc, 0.0);
    }

    #[test]
    fn one_wm2_day_melts_about_a_quarter_mm() {
        let p = EnergyBalanceParams::default();
        let mut cc = 0.0;
        let melt = apply_energy(&p, 1.0, 1.0, 1000.0, &mut cc);
        assert!((melt - DAY_S / L_FUSION).abs() < 1e-12);
        assert!((melt - 0.2586).abs() < 1e-3);
    }
}
