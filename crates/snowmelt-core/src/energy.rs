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
/// Gravitational acceleration (m s⁻²).
const G_ACCEL: f64 = 9.81;
/// Von Kármán constant.
const VON_KARMAN: f64 = 0.4;
/// Critical bulk Richardson number above which turbulence is fully damped.
const RI_CRIT: f64 = 0.2;

/// Explicit aerodynamic resistance for the turbulent fluxes.
///
/// Replaces the single bulk exchange coefficient with a Monin–Obukhov-style
/// neutral aerodynamic resistance from a logarithmic wind profile,
///
/// ```text
/// r_a = ln(z/z0m)·ln(z/z0h) / (k²·u)      [s m⁻¹]
/// ```
///
/// so the turbulent conductance `1/r_a` follows physically from the surface
/// roughness instead of a tuned constant. With `stability` on, the
/// conductance is scaled by a bulk-Richardson stability function (Anderson
/// 1976; Tarboton & Luce 1996): stable stratification (warm air over cold
/// snow) damps the fluxes, unstable enhances them — the dominant control on
/// snow sublimation that a fixed coefficient misses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AeroResistance {
    /// Momentum roughness length `z0m` (m). Snow: ~1e-4 to 5e-3.
    pub z0_momentum: f64,
    /// Scalar (heat/vapour) roughness length `z0h` (m). Often `z0m/10` or
    /// smaller over snow.
    pub z0_heat: f64,
    /// Measurement height of wind, temperature and humidity (m). Typical: 2.
    pub measurement_height: f64,
    /// Apply the bulk-Richardson stability correction.
    pub stability: bool,
}

impl Default for AeroResistance {
    fn default() -> Self {
        Self {
            z0_momentum: 1e-3,
            z0_heat: 1e-4,
            measurement_height: 2.0,
            stability: true,
        }
    }
}

impl AeroResistance {
    /// Checks roughness lengths are positive and below the measurement
    /// height (so the log profile is finite and positive).
    ///
    /// # Errors
    /// Returns [`SnowmeltError::InvalidParameter`] on the first violation.
    pub fn validate(&self) -> Result<()> {
        let checks: [(&'static str, f64, bool); 3] = [
            (
                "z0_momentum",
                self.z0_momentum,
                self.z0_momentum.is_finite()
                    && self.z0_momentum > 0.0
                    && self.z0_momentum < self.measurement_height,
            ),
            (
                "z0_heat",
                self.z0_heat,
                self.z0_heat.is_finite()
                    && self.z0_heat > 0.0
                    && self.z0_heat < self.measurement_height,
            ),
            (
                "measurement_height",
                self.measurement_height,
                self.measurement_height.is_finite() && self.measurement_height > 0.0,
            ),
        ];
        for (name, value, ok) in checks {
            if !ok {
                return Err(SnowmeltError::InvalidParameter {
                    name,
                    reason: format!("out of domain: {value}"),
                });
            }
        }
        Ok(())
    }

    /// Neutral turbulent conductance `1/r_a` (m s⁻¹) at wind speed `u`.
    fn neutral_conductance(&self, u: f64) -> f64 {
        let u = u.max(0.1); // floor to keep a finite resistance in calm air
        let ln_m = (self.measurement_height / self.z0_momentum).ln();
        let ln_h = (self.measurement_height / self.z0_heat).ln();
        VON_KARMAN * VON_KARMAN * u / (ln_m * ln_h)
    }
}

/// Bulk-Richardson stability factor for the turbulent conductance.
///
/// `1` at neutral; `< 1` (down to 0 at [`RI_CRIT`]) for stable
/// stratification `Ri > 0`; `> 1` for unstable `Ri < 0`.
fn stability_factor(ri: f64) -> f64 {
    if ri >= 0.0 {
        let ri = ri.min(RI_CRIT);
        let f = 1.0 - ri / RI_CRIT;
        f * f
    } else {
        (1.0 - 16.0 * ri).powf(0.75)
    }
}

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
    /// Typical: 0.001–0.003. Ignored when `aerodynamic` is `Some`.
    pub exchange_coeff: f64,
    /// Explicit aerodynamic resistance for the turbulent fluxes. When set,
    /// it replaces the `exchange_coeff` bulk path (roughness-based
    /// conductance plus optional stability correction).
    pub aerodynamic: Option<AeroResistance>,
    /// Ground heat flux into the pack (W m⁻²). Typical: 0–2.
    pub ground_heat: f64,
    /// Maximum pack cooling below 0 °C (K) used to cap the cold content
    /// at `c_ice·SWE·t_cold_max`. Typical: 5–15.
    pub t_cold_max: f64,
    /// Effective cloud fraction (0–1). Reduces shortwave by
    /// `(1 − 0.75·N³)` and raises incoming longwave emissivity by
    /// `(1 + 0.22·N²)` (Crawford & Duchon-style). 0 = clear sky.
    pub cloud_fraction: f64,
}

impl Default for EnergyBalanceParams {
    fn default() -> Self {
        Self {
            wind_speed: 2.0,
            rel_humidity: 0.6,
            snow_emissivity: 0.98,
            exchange_coeff: 0.0015,
            aerodynamic: None,
            ground_heat: 1.0,
            t_cold_max: 10.0,
            cloud_fraction: 0.0,
        }
    }
}

impl EnergyBalanceParams {
    /// Checks that all parameters are finite and within their domains.
    ///
    /// # Errors
    /// Returns [`SnowmeltError::InvalidParameter`] on the first violation.
    pub fn validate(&self) -> Result<()> {
        let checks: [(&'static str, f64, bool); 7] = [
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
            (
                "cloud_fraction",
                self.cloud_fraction,
                (0.0..=1.0).contains(&self.cloud_fraction),
            ),
        ];
        for (name, value, in_domain) in checks {
            if !value.is_finite() || !in_domain {
                return Err(SnowmeltError::InvalidParameter {
                    name,
                    reason: format!("out of domain: {value}"),
                });
            }
        }
        if let Some(aero) = &self.aerodynamic {
            aero.validate()?;
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

/// Specific heat of liquid water (J kg⁻¹ K⁻¹).
const C_WATER: f64 = 4186.0;

/// Energy fluxes into the snowpack for one cell: `(total, latent)` W m⁻².
///
/// `sw_in` is the clear-sky incoming shortwave (W m⁻², daily mean),
/// `albedo` the current surface albedo, `pressure_pa` the cell's air
/// pressure, `rain_mm` the step's liquid precipitation (rain-on-snow
/// advective heat) over `dt_days`. The latent component is returned
/// separately so the caller can account sublimation mass loss.
pub(crate) fn energy_fluxes(
    p: &EnergyBalanceParams,
    t_air_c: f64,
    sw_in: f64,
    albedo: f64,
    pressure_pa: f64,
    rain_mm: f64,
    dt_days: f64,
) -> (f64, f64) {
    let t_air_k = t_air_c + 273.15;
    // Snow surface cannot exceed the melting point.
    let t_surf_c = t_air_c.min(0.0);
    let t_surf_k = t_surf_c + 273.15;
    let n = p.cloud_fraction;

    // Shortwave, attenuated by clouds.
    let q_sw = (1.0 - albedo) * sw_in * (1.0 - 0.75 * n.powi(3));

    // Longwave: Brutsaert (1975) clear-sky emissivity, raised by clouds.
    let e_air_hpa = p.rel_humidity * e_sat_hpa(t_air_c);
    let emiss_air = (1.24 * (e_air_hpa / t_air_k).powf(1.0 / 7.0)) * (1.0 + 0.22 * n.powi(2));
    let emiss_air = emiss_air.min(1.0);
    let q_lw = emiss_air * SIGMA * t_air_k.powi(4) - p.snow_emissivity * SIGMA * t_surf_k.powi(4);

    // Turbulent fluxes: either a fixed bulk coefficient (`C_e·u`) or an
    // explicit roughness-based conductance (`1/r_a`) with optional
    // bulk-Richardson stability correction.
    let rho_air = pressure_pa / (R_DRY * t_air_k);
    let q_air = 0.622 * (e_air_hpa * 100.0) / pressure_pa;
    let q_surf = 0.622 * (e_sat_hpa(t_surf_c) * 100.0) / pressure_pa;
    let conductance = match &p.aerodynamic {
        Some(aero) => {
            let mut g = aero.neutral_conductance(p.wind_speed);
            if aero.stability {
                let u = p.wind_speed.max(0.1);
                let ri =
                    G_ACCEL * aero.measurement_height * (t_air_c - t_surf_c) / (t_air_k * u * u);
                g *= stability_factor(ri);
            }
            g
        }
        None => p.exchange_coeff * p.wind_speed,
    };
    let q_h = rho_air * CP_AIR * conductance * (t_air_c - t_surf_c);
    let q_e = rho_air * L_SUBLIMATION * conductance * (q_air - q_surf);

    // Rain-on-snow advective heat (rain at air temperature onto a 0 °C pack).
    let q_r = C_WATER * rain_mm * t_air_c.max(0.0) / (dt_days * DAY_S);

    (q_sw + q_lw + q_h + q_e + q_r + p.ground_heat, q_e)
}

/// Sublimation mass loss (mm w.e.) for a negative latent flux over the
/// step, to be capped by the available SWE by the caller.
pub(crate) fn sublimation_mm(q_e_wm2: f64, dt_days: f64) -> f64 {
    if q_e_wm2 >= 0.0 {
        0.0 // deposition gain is neglected
    } else {
        -q_e_wm2 * dt_days * DAY_S / L_SUBLIMATION
    }
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
        let (q, _) = energy_fluxes(&p, 5.0, 250.0, 0.6, air_pressure_pa(2000.0), 0.0, 1.0);
        // SW neto = 100 W/m²; LW pierde, turbulentos aportan algo.
        assert!(q > 30.0 && q < 250.0, "q = {q}");
    }

    #[test]
    fn cold_clear_night_yields_negative_energy() {
        let p = EnergyBalanceParams::default();
        let (q, q_e) = energy_fluxes(&p, -15.0, 0.0, 0.8, air_pressure_pa(3000.0), 0.0, 1.0);
        assert!(q < 0.0, "q = {q}");
        // Aire seco sobre nieve saturada a la misma T → sublimación.
        assert!(q_e < 0.0, "q_e = {q_e}");
        assert!(sublimation_mm(q_e, 1.0) > 0.0);
        assert_eq!(sublimation_mm(10.0, 1.0), 0.0);
    }

    #[test]
    fn clouds_warm_winter_and_cool_spring() {
        let clear = EnergyBalanceParams::default();
        let cloudy = EnergyBalanceParams {
            cloud_fraction: 0.5,
            ..Default::default()
        };
        let press = air_pressure_pa(3000.0);
        // Invierno: poca SW → domina la LW extra de las nubes.
        let (q_clear_w, _) = energy_fluxes(&clear, -2.0, 60.0, 0.8, press, 0.0, 1.0);
        let (q_cloud_w, _) = energy_fluxes(&cloudy, -2.0, 60.0, 0.8, press, 0.0, 1.0);
        assert!(q_cloud_w > q_clear_w);
        // Primavera: mucha SW y albedo bajo → domina la atenuación SW.
        let (q_clear_s, _) = energy_fluxes(&clear, 8.0, 380.0, 0.45, press, 0.0, 1.0);
        let (q_cloud_s, _) = energy_fluxes(&cloudy, 8.0, 380.0, 0.45, press, 0.0, 1.0);
        assert!(q_cloud_s < q_clear_s);
    }

    #[test]
    fn rain_on_snow_adds_heat() {
        let p = EnergyBalanceParams::default();
        let press = air_pressure_pa(2000.0);
        let (dry, _) = energy_fluxes(&p, 4.0, 100.0, 0.6, press, 0.0, 1.0);
        let (wet, _) = energy_fluxes(&p, 4.0, 100.0, 0.6, press, 30.0, 1.0);
        // 30 mm a 4 °C ≈ 4186·30·4/86400 ≈ 5.8 W/m².
        assert!((wet - dry - 5.8).abs() < 0.1, "ΔQ_R = {}", wet - dry);
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
    fn aero_resistance_validates_and_rejects_bad_roughness() {
        AeroResistance::default().validate().unwrap();
        for bad in [
            AeroResistance {
                z0_momentum: 0.0,
                ..Default::default()
            },
            AeroResistance {
                z0_heat: -1e-4,
                ..Default::default()
            },
            AeroResistance {
                z0_momentum: 5.0, // above measurement height
                ..Default::default()
            },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
        // Through EnergyBalanceParams too.
        let p = EnergyBalanceParams {
            aerodynamic: Some(AeroResistance {
                z0_momentum: 0.0,
                ..Default::default()
            }),
            ..Default::default()
        };
        assert!(p.validate().is_err());
    }

    #[test]
    fn stability_factor_neutral_stable_unstable() {
        assert!((stability_factor(0.0) - 1.0).abs() < 1e-12);
        assert!(stability_factor(0.1) < 1.0); // stable damps
        assert_eq!(stability_factor(RI_CRIT), 0.0); // fully damped at critical
        assert_eq!(stability_factor(1.0), 0.0); // capped beyond critical
        assert!(stability_factor(-0.5) > 1.0); // unstable enhances
    }

    #[test]
    fn rougher_surface_increases_turbulent_exchange() {
        let press = air_pressure_pa(3000.0);
        let smooth = EnergyBalanceParams {
            aerodynamic: Some(AeroResistance {
                z0_momentum: 1e-4,
                z0_heat: 1e-5,
                stability: false,
                ..Default::default()
            }),
            ..Default::default()
        };
        let rough = EnergyBalanceParams {
            aerodynamic: Some(AeroResistance {
                z0_momentum: 5e-3,
                z0_heat: 5e-4,
                stability: false,
                ..Default::default()
            }),
            ..Default::default()
        };
        // Cold dry air over snow: latent flux (sublimation) is negative;
        // a rougher surface makes it more negative (stronger exchange).
        let (_, qe_smooth) = energy_fluxes(&smooth, -8.0, 100.0, 0.8, press, 0.0, 1.0);
        let (_, qe_rough) = energy_fluxes(&rough, -8.0, 100.0, 0.8, press, 0.0, 1.0);
        assert!(qe_smooth < 0.0 && qe_rough < 0.0);
        assert!(
            qe_rough < qe_smooth,
            "rough {qe_rough} vs smooth {qe_smooth}"
        );
    }

    #[test]
    fn stability_damps_sensible_flux_on_warm_air() {
        let press = air_pressure_pa(2000.0);
        let neutral = EnergyBalanceParams {
            aerodynamic: Some(AeroResistance {
                stability: false,
                ..Default::default()
            }),
            ..Default::default()
        };
        let stable = EnergyBalanceParams {
            aerodynamic: Some(AeroResistance {
                stability: true,
                ..Default::default()
            }),
            ..Default::default()
        };
        // Warm air (8 °C) over a 0 °C snow surface → stable stratification.
        // Isolate Q_H by comparing total energy at the same low radiation;
        // the stable run must transfer less sensible heat.
        let (q_neu, _) = energy_fluxes(&neutral, 8.0, 50.0, 0.6, press, 0.0, 1.0);
        let (q_sta, _) = energy_fluxes(&stable, 8.0, 50.0, 0.6, press, 0.0, 1.0);
        assert!(q_sta < q_neu, "stable {q_sta} should be < neutral {q_neu}");
    }

    #[test]
    fn aero_and_bulk_give_comparable_magnitude() {
        let press = air_pressure_pa(2500.0);
        let bulk = EnergyBalanceParams::default();
        let aero = EnergyBalanceParams {
            aerodynamic: Some(AeroResistance {
                stability: false,
                ..Default::default()
            }),
            ..Default::default()
        };
        let (q_bulk, _) = energy_fluxes(&bulk, 5.0, 200.0, 0.6, press, 0.0, 1.0);
        let (q_aero, _) = energy_fluxes(&aero, 5.0, 200.0, 0.6, press, 0.0, 1.0);
        // Same sign and within a factor of ~3 (roughness vs tuned constant).
        assert!(q_bulk > 0.0 && q_aero > 0.0);
        assert!(
            (q_aero / q_bulk) > 0.3 && (q_aero / q_bulk) < 3.0,
            "{q_aero} vs {q_bulk}"
        );
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
