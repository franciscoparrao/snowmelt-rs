//! Topographic downscaling of meteorological forcing (MicroMet-style).
//!
//! Turns a single reference temperature and precipitation value into
//! per-cell fields at the DEM resolution, using terrain to resolve
//! sub-grid structure that coarse reanalysis (5–25 km) misses over a
//! high-relief catchment. The scheme follows the intermediate-complexity
//! approach of Liston & Elder (2006, *MicroMet*), with three terrain
//! controls beyond the plain lapse rate:
//!
//! 1. **Temperature** — lapse extrapolation plus a curvature term that
//!    cools concave terrain and warms convex terrain, a daily-mean proxy
//!    for cold-air pooling / drainage in valleys (`temp_curvature`).
//! 2. **Wind** — a terrain speed factor `W = 1 + γ_s·Ω_s + γ_c·Ω_c`, where
//!    `Ω_s` is the slope in the wind direction and `Ω_c` the curvature,
//!    both in `[-0.5, 0.5]` (Liston & Elder 2006, eqs. 16–19).
//! 3. **Precipitation** — an elevation factor (Thornton et al. 1997) times
//!    an orographic windward factor `1 + γ_w·Ω_s` that enhances slopes
//!    facing the prevailing wind and dries leeward slopes.
//!
//! The windward precipitation term and the curvature temperature term are
//! the levers a purely elevation-based reanalysis downscaling cannot reach;
//! they are what this module adds over the existing lapse/`precip_gradient`
//! forcing. Output grids feed the model through
//! [`Forcing::Distributed`](crate::Forcing::Distributed), so the snow model
//! and energy balance are untouched.

use ndarray::{Array2, ArrayView2, Zip};

use crate::error::{Result, SnowmeltError};
use crate::terrain;

/// Parameters of the topographic downscaling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DownscaleParams {
    /// Temperature lapse rate (°C m⁻¹). Typical: −0.0065 to −0.0075.
    pub lapse_rate: f64,
    /// Curvature temperature coefficient (°C). Adds `temp_curvature·Ω_c`,
    /// so convex terrain (ridges, `Ω_c > 0`) is warmer and concave terrain
    /// (valleys, `Ω_c < 0`) is colder — a daily-mean proxy for cold-air
    /// pooling. `Ω_c ∈ [−0.5, 0.5]`, so the swing is `±temp_curvature/2`.
    /// 0 disables it. Typical: 1–4.
    pub temp_curvature: f64,
    /// Precipitation–elevation factor (km⁻¹), Thornton et al. (1997):
    /// `P(z) = P_ref·(1 + f·Δz)/(1 − f·Δz)` with `Δz = (z − z_ref)` in km
    /// (the ratio is clamped against its pole). 0 disables it.
    /// Typical: 0.05–0.3.
    pub precip_elev_factor: f64,
    /// Orographic windward factor `γ_w`: precipitation is scaled by
    /// `1 + γ_w·Ω_s`, enhancing slopes facing the wind (`Ω_s > 0`) and
    /// drying leeward slopes. 0 disables it. Typical: 0.2–0.8.
    pub precip_windward: f64,
    /// Prevailing wind direction, degrees the wind blows **from** (compass,
    /// clockwise from north). Central-Chile frontal precipitation is
    /// north-westerly, ~300°.
    pub wind_dir_from_deg: f64,
    /// Weight of the slope-in-wind term in the wind speed factor (`γ_s`).
    /// Typical: 0.5.
    pub wind_slope_weight: f64,
    /// Weight of the curvature term in the wind speed factor (`γ_c`).
    /// Typical: 0.5.
    pub wind_curvature_weight: f64,
}

impl Default for DownscaleParams {
    fn default() -> Self {
        Self {
            lapse_rate: -0.0065,
            temp_curvature: 0.0,
            precip_elev_factor: 0.0,
            precip_windward: 0.0,
            wind_dir_from_deg: 300.0,
            wind_slope_weight: 0.5,
            wind_curvature_weight: 0.5,
        }
    }
}

impl DownscaleParams {
    /// Checks that every parameter is finite and within its domain.
    ///
    /// # Errors
    /// Returns [`SnowmeltError::InvalidParameter`] on the first violation.
    pub fn validate(&self) -> Result<()> {
        let checks: [(&'static str, f64, bool); 7] = [
            ("lapse_rate", self.lapse_rate, self.lapse_rate.is_finite()),
            (
                "temp_curvature",
                self.temp_curvature,
                self.temp_curvature.is_finite(),
            ),
            (
                "precip_elev_factor",
                self.precip_elev_factor,
                self.precip_elev_factor.is_finite() && self.precip_elev_factor >= 0.0,
            ),
            (
                "precip_windward",
                self.precip_windward,
                self.precip_windward.is_finite() && self.precip_windward >= 0.0,
            ),
            (
                "wind_dir_from_deg",
                self.wind_dir_from_deg,
                self.wind_dir_from_deg.is_finite(),
            ),
            (
                "wind_slope_weight",
                self.wind_slope_weight,
                self.wind_slope_weight.is_finite() && self.wind_slope_weight >= 0.0,
            ),
            (
                "wind_curvature_weight",
                self.wind_curvature_weight,
                self.wind_curvature_weight.is_finite() && self.wind_curvature_weight >= 0.0,
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
}

/// Largest `|f·Δz|` allowed in the precipitation–elevation ratio, to keep
/// it away from the `1/(1 − f·Δz)` pole over a tall catchment.
const PRECIP_RATIO_CLAMP: f64 = 0.95;

/// Topographic downscaler: precomputes terrain derivatives once, then maps
/// scalar reference forcing to per-cell grids.
#[derive(Debug)]
pub struct Downscaler {
    elevation: Array2<f64>,
    z_ref: f64,
    params: DownscaleParams,
    /// Normalised curvature `Ω_c ∈ [−0.5, 0.5]` (`NaN` on nodata).
    curvature: Array2<f64>,
    /// Slope in the wind direction `Ω_s ∈ [−0.5, 0.5]` (`NaN` on nodata).
    omega_s: Array2<f64>,
    /// Terrain wind speed factor `W` (`NaN` on nodata).
    wind_factor: Array2<f64>,
}

impl Downscaler {
    /// Builds a downscaler for `elevation` (m, `NaN` = nodata) on a grid of
    /// spacing `cellsize` (m), with the reference forcing measured at
    /// `z_ref` (m).
    ///
    /// # Errors
    /// Returns an error if `params` fail [`DownscaleParams::validate`].
    pub fn new(
        elevation: Array2<f64>,
        cellsize: f64,
        z_ref: f64,
        params: DownscaleParams,
    ) -> Result<Self> {
        params.validate()?;
        let (slope, aspect) = terrain::slope_aspect(&elevation, cellsize);
        let curvature = terrain::curvature(&elevation);

        let max_slope = slope
            .iter()
            .filter(|s| s.is_finite())
            .fold(0.0_f64, |m, &s| m.max(s));
        // Liston & Elder normalise slope by twice the domain maximum so the
        // slope-in-wind term lands in [−0.5, 0.5].
        let slope_scale = if max_slope > 0.0 {
            1.0 / (2.0 * max_slope)
        } else {
            0.0
        };
        let theta = params.wind_dir_from_deg.to_radians();
        let (gamma_s, gamma_c) = (params.wind_slope_weight, params.wind_curvature_weight);

        let mut omega_s = Array2::zeros(elevation.dim());
        let mut wind_factor = Array2::zeros(elevation.dim());
        Zip::from(&mut omega_s)
            .and(&mut wind_factor)
            .and(&slope)
            .and(&aspect)
            .and(&curvature)
            .for_each(|os, wf, &s, &a, &c| {
                if !s.is_finite() || !a.is_finite() || !c.is_finite() {
                    *os = f64::NAN;
                    *wf = f64::NAN;
                    return;
                }
                // Ω_s = ŝ·cos(θ_wind − aspect): positive on slopes facing the
                // wind (aspect ≈ wind-from direction), negative leeward.
                let s_hat = s * slope_scale;
                let omega = s_hat * (theta - a).cos();
                *os = omega;
                *wf = 1.0 + gamma_s * omega + gamma_c * c;
            });

        Ok(Self {
            elevation,
            z_ref,
            params,
            curvature,
            omega_s,
            wind_factor,
        })
    }

    /// Per-cell temperature (°C) for a reference value `t_ref` measured at
    /// `z_ref`: lapse extrapolation plus the curvature term. `NaN` on nodata.
    pub fn temperature(&self, t_ref: f64) -> Array2<f64> {
        let lapse = self.params.lapse_rate;
        let kappa = self.params.temp_curvature;
        let z_ref = self.z_ref;
        let mut out = Array2::zeros(self.elevation.dim());
        Zip::from(&mut out)
            .and(&self.elevation)
            .and(&self.curvature)
            .par_for_each(|t, &z, &c| {
                *t = if z.is_finite() && c.is_finite() {
                    t_ref + lapse * (z - z_ref) + kappa * c
                } else {
                    f64::NAN
                };
            });
        out
    }

    /// Per-cell precipitation (mm) for a reference value `p_ref`: the
    /// elevation factor times the orographic windward factor, clamped to
    /// `≥ 0`. `NaN` on nodata.
    pub fn precip(&self, p_ref: f64) -> Array2<f64> {
        let f = self.params.precip_elev_factor;
        let gamma_w = self.params.precip_windward;
        let z_ref = self.z_ref;
        let mut out = Array2::zeros(self.elevation.dim());
        Zip::from(&mut out)
            .and(&self.elevation)
            .and(&self.omega_s)
            .par_for_each(|p, &z, &os| {
                if !z.is_finite() || !os.is_finite() {
                    *p = f64::NAN;
                    return;
                }
                let dz_km = (z - z_ref) / 1000.0;
                let x = (f * dz_km).clamp(-PRECIP_RATIO_CLAMP, PRECIP_RATIO_CLAMP);
                let f_elev = (1.0 + x) / (1.0 - x);
                let f_wind = (1.0 + gamma_w * os).max(0.0);
                *p = (p_ref * f_elev * f_wind).max(0.0);
            });
        out
    }

    /// The terrain wind speed factor `W` per cell (`NaN` on nodata).
    pub fn wind_factor(&self) -> ArrayView2<'_, f64> {
        self.wind_factor.view()
    }

    /// The slope-in-wind term `Ω_s` per cell (`NaN` on nodata).
    pub fn omega_s(&self) -> ArrayView2<'_, f64> {
        self.omega_s.view()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn flat_terrain_reduces_to_uniform_lapse() {
        let z = Array2::from_elem((4, 4), 2000.0);
        let p = DownscaleParams {
            temp_curvature: 3.0,
            precip_elev_factor: 0.2,
            precip_windward: 0.5,
            ..DownscaleParams::default()
        };
        let ds = Downscaler::new(z, 200.0, 2000.0, p).unwrap();
        // Curvature 0 and Ω_s 0 on flat ground: temp = t_ref, precip = p_ref.
        assert!(ds.temperature(5.0).iter().all(|&t| approx(t, 5.0)));
        assert!(ds.precip(10.0).iter().all(|&v| approx(v, 10.0)));
    }

    #[test]
    fn lapse_only_matches_plain_extrapolation() {
        let z = Array2::from_shape_fn((5, 5), |(i, j)| 1000.0 + 100.0 * (i + j) as f64);
        let p = DownscaleParams {
            lapse_rate: -0.0075,
            ..DownscaleParams::default()
        };
        let ds = Downscaler::new(z.clone(), 200.0, 1500.0, p).unwrap();
        let t = ds.temperature(0.0);
        // Curvature is 0 on the planar interior, so temp = lapse·(z − z_ref).
        let i = 2;
        let j = 2;
        let expect = -0.0075 * (z[[i, j]] - 1500.0);
        assert!(approx(t[[i, j]], expect), "{} vs {expect}", t[[i, j]]);
    }

    #[test]
    fn curvature_term_warms_ridge_cools_valley() {
        // V-shaped valley along the central column.
        let z = Array2::from_shape_fn((5, 5), |(_, j)| 2000.0 + 50.0 * (j as f64 - 2.0).abs());
        let p = DownscaleParams {
            lapse_rate: 0.0, // isolate the curvature term
            temp_curvature: 4.0,
            ..DownscaleParams::default()
        };
        let ds = Downscaler::new(z, 200.0, 2000.0, p).unwrap();
        let t = ds.temperature(0.0);
        // Valley floor (concave) colder than its convex flank.
        assert!(t[[2, 2]] < t[[2, 1]], "{} !< {}", t[[2, 2]], t[[2, 1]]);
    }

    #[test]
    fn windward_slope_gets_more_precip_than_leeward() {
        // West-facing slope: z increases eastward, so it descends toward the
        // west (aspect 270°). Wind from the west (270°) hits it head-on.
        let z = Array2::from_shape_fn((5, 5), |(_, j)| 2000.0 + 30.0 * j as f64);
        let west = DownscaleParams {
            wind_dir_from_deg: 270.0,
            precip_windward: 0.6,
            ..DownscaleParams::default()
        };
        let ds_w = Downscaler::new(z.clone(), 200.0, 2000.0, west).unwrap();
        // Same slope, wind from the east (90°): now leeward.
        let east = DownscaleParams {
            wind_dir_from_deg: 90.0,
            precip_windward: 0.6,
            ..DownscaleParams::default()
        };
        let ds_e = Downscaler::new(z, 200.0, 2000.0, east).unwrap();
        let pw = ds_w.precip(10.0)[[2, 2]];
        let pe = ds_e.precip(10.0)[[2, 2]];
        assert!(pw > 10.0, "windward should exceed reference: {pw}");
        assert!(pe < 10.0, "leeward should fall below reference: {pe}");
    }

    #[test]
    fn precip_elevation_factor_increases_with_height() {
        let z = Array2::from_shape_fn((3, 3), |(i, _)| 2000.0 + 1000.0 * i as f64);
        let p = DownscaleParams {
            precip_elev_factor: 0.1,
            ..DownscaleParams::default()
        };
        let ds = Downscaler::new(z, 200.0, 2000.0, p).unwrap();
        let pr = ds.precip(10.0);
        // Higher row (interior) gets more precip than the reference level.
        assert!(pr[[1, 1]] > 10.0, "{}", pr[[1, 1]]);
    }

    #[test]
    fn precip_never_negative_even_with_strong_windward_drying() {
        let z = Array2::from_shape_fn((5, 5), |(_, j)| 2000.0 + 100.0 * j as f64);
        let p = DownscaleParams {
            wind_dir_from_deg: 90.0, // leeward on this west-descending slope
            precip_windward: 5.0,    // 1 + γ_w·Ω_s could go negative → clamp
            ..DownscaleParams::default()
        };
        let ds = Downscaler::new(z, 200.0, 2000.0, p).unwrap();
        assert!(ds.precip(10.0).iter().all(|&v| v.is_nan() || v >= 0.0));
    }

    #[test]
    fn wind_factor_higher_on_windward_convex_terrain() {
        let z = Array2::from_shape_fn((5, 5), |(_, j)| 2000.0 + 30.0 * j as f64);
        let p = DownscaleParams {
            wind_dir_from_deg: 270.0,
            ..DownscaleParams::default()
        };
        let ds = Downscaler::new(z, 200.0, 2000.0, p).unwrap();
        // On a planar windward slope Ω_s > 0 → W > 1.
        assert!(ds.wind_factor()[[2, 2]] > 1.0);
        assert!(ds.omega_s()[[2, 2]] > 0.0);
    }

    #[test]
    fn nodata_propagates_through_downscaling() {
        let mut z = Array2::from_elem((4, 4), 2000.0);
        z[[1, 1]] = f64::NAN;
        let ds = Downscaler::new(z, 200.0, 2000.0, DownscaleParams::default()).unwrap();
        assert!(ds.temperature(5.0)[[1, 1]].is_nan());
        assert!(ds.precip(10.0)[[1, 1]].is_nan());
    }

    #[test]
    fn rejects_invalid_params() {
        let z = Array2::from_elem((3, 3), 2000.0);
        let bad = DownscaleParams {
            precip_windward: -1.0,
            ..DownscaleParams::default()
        };
        assert!(Downscaler::new(z, 200.0, 2000.0, bad).is_err());
    }
}
