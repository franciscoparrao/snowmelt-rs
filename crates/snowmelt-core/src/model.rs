//! Snowpack state and time integration.

use ndarray::{Array2, ArrayView2, Zip};

use crate::dem::Dem;
use crate::energy;
use crate::error::{Result, SnowmeltError};
use crate::forcing::Forcing;
use crate::params::DegreeDayParams;

/// SWE (mm) above which a cell counts as snow-covered in
/// [`SnowModel::summarize`]. Filters out floating-point melt residues.
pub const SNOW_COVER_THRESHOLD_MM: f64 = 0.1;

/// Per-cell fluxes produced by one model step (all in mm w.e.).
///
/// Mass balance per cell and step: `snowfall + rain == precip` and
/// `Δswe == snowfall - melt - sublimation`. Nodata cells are `NaN`.
#[derive(Debug)]
pub struct StepOutput {
    /// Solid precipitation added to the snowpack (mm).
    pub snowfall: Array2<f64>,
    /// Liquid precipitation, passed through to runoff (mm).
    pub rain: Array2<f64>,
    /// Snowmelt released from the snowpack (mm).
    pub melt: Array2<f64>,
    /// Sublimation mass loss (mm; non-zero only in energy-balance mode).
    pub sublimation: Array2<f64>,
}

impl StepOutput {
    /// Liquid water reaching the ground during the step: `rain + melt` (mm).
    pub fn runoff(&self) -> Array2<f64> {
        &self.rain + &self.melt
    }
}

/// Spatially aggregated fluxes and state for one step (means over valid cells).
#[derive(Debug, Clone, Copy)]
pub struct StepSummary {
    /// Mean snowfall (mm).
    pub mean_snowfall: f64,
    /// Mean rain (mm).
    pub mean_rain: f64,
    /// Mean melt (mm).
    pub mean_melt: f64,
    /// Mean sublimation loss (mm; energy-balance mode).
    pub mean_sublimation: f64,
    /// Mean runoff = rain + melt (mm).
    pub mean_runoff: f64,
    /// Mean SWE after the step (mm).
    pub mean_swe: f64,
    /// Mean albedo used by the radiative term (constant if no decay).
    pub mean_albedo: f64,
    /// Fraction of valid cells with SWE > [`SNOW_COVER_THRESHOLD_MM`].
    pub snow_cover_fraction: f64,
}

/// Distributed degree-day snow model.
///
/// Holds the per-cell SWE state (mm w.e.) and advances it one forcing step
/// at a time. Cell updates run in parallel with Rayon.
#[derive(Debug)]
pub struct SnowModel {
    dem: Dem,
    params: DegreeDayParams,
    swe: Array2<f64>,
    /// Days since the last significant snowfall (only meaningful with
    /// `params.albedo_decay`; initial snow counts as fresh).
    age: Array2<f64>,
    /// Per-cell albedo used by the radiative term in the current step.
    albedo_buf: Array2<f64>,
    /// Cold content of the pack (J m⁻²; energy-balance mode only).
    cold_content: Array2<f64>,
    /// Air pressure per cell (Pa), derived once from elevation.
    pressure_buf: Array2<f64>,
    /// Scratch: `(total, latent)` energy fluxes (W m⁻²) in energy-balance mode.
    energy_buf: Array2<[f64; 2]>,
    temp_buf: Array2<f64>,
    precip_buf: Array2<f64>,
    rad_zero: Array2<f64>,
}

impl SnowModel {
    /// Creates a model with zero initial SWE on every valid DEM cell.
    ///
    /// # Errors
    /// Returns an error if `params` fail [`DegreeDayParams::validate`].
    pub fn new(dem: Dem, params: DegreeDayParams) -> Result<Self> {
        params.validate()?;
        let swe = dem
            .elevation()
            .mapv(|z| if z.is_finite() { 0.0 } else { f64::NAN });
        let age = swe.clone();
        let initial_albedo = match params.albedo_decay {
            Some(decay) => decay.albedo_fresh,
            None => params.albedo,
        };
        let albedo_buf = dem.elevation().mapv(|z| {
            if z.is_finite() {
                initial_albedo
            } else {
                f64::NAN
            }
        });
        let cold_content = swe.clone();
        let pressure_buf = dem.elevation().mapv(|z| {
            if z.is_finite() {
                energy::air_pressure_pa(z)
            } else {
                f64::NAN
            }
        });
        let energy_buf = Array2::from_elem(dem.shape(), [0.0, 0.0]);
        let temp_buf = Array2::zeros(dem.shape());
        let precip_buf = Array2::zeros(dem.shape());
        let rad_zero = Array2::zeros(dem.shape());
        Ok(Self {
            dem,
            params,
            swe,
            age,
            albedo_buf,
            cold_content,
            pressure_buf,
            energy_buf,
            temp_buf,
            precip_buf,
            rad_zero,
        })
    }

    /// Creates a model with a prescribed initial SWE grid (mm).
    ///
    /// Nodata DEM cells are forced to `NaN` regardless of `initial_swe`.
    ///
    /// # Errors
    /// Returns an error if the shape differs from the DEM, or if any valid
    /// DEM cell has a negative or non-finite initial SWE.
    pub fn with_initial_swe(
        dem: Dem,
        params: DegreeDayParams,
        initial_swe: Array2<f64>,
    ) -> Result<Self> {
        let mut model = Self::new(dem, params)?;
        if initial_swe.dim() != model.dem.shape() {
            return Err(SnowmeltError::ShapeMismatch {
                expected: model.dem.shape(),
                got: initial_swe.dim(),
            });
        }
        for (&z, &s) in model.dem.elevation().iter().zip(initial_swe.iter()) {
            if z.is_finite() && !(s.is_finite() && s >= 0.0) {
                return Err(SnowmeltError::InvalidParameter {
                    name: "initial_swe",
                    reason: format!("valid cells need finite SWE >= 0, got {s}"),
                });
            }
        }
        Zip::from(&mut model.swe)
            .and(&initial_swe)
            .for_each(|swe, &s| {
                if swe.is_finite() {
                    *swe = s;
                }
            });
        Ok(model)
    }

    /// The model's DEM.
    pub fn dem(&self) -> &Dem {
        &self.dem
    }

    /// The model's parameters.
    pub fn params(&self) -> &DegreeDayParams {
        &self.params
    }

    /// Current SWE grid (mm; `NaN` = nodata).
    pub fn swe(&self) -> ArrayView2<'_, f64> {
        self.swe.view()
    }

    /// Days since the last significant snowfall per cell. Only updated
    /// when `params.albedo_decay` is set; initial snow counts as fresh.
    pub fn snow_age(&self) -> ArrayView2<'_, f64> {
        self.age.view()
    }

    /// Per-cell albedo used by the radiative term in the last step
    /// (constant grid when `params.albedo_decay` is `None`).
    pub fn albedo(&self) -> ArrayView2<'_, f64> {
        self.albedo_buf.view()
    }

    /// Cold content of the pack (J m⁻²; only updated in energy-balance
    /// mode). Zero means an isothermal pack at 0 °C.
    pub fn cold_content(&self) -> ArrayView2<'_, f64> {
        self.cold_content.view()
    }

    /// Advances the model by one step of `dt_days` days, optionally driven
    /// by a shortwave radiation grid (W m⁻², daily mean over the step).
    ///
    /// The radiation grid feeds the enhanced temperature-index melt term
    /// (see [`DegreeDayParams::srf`]); it is required when `srf > 0` and
    /// ignored when `srf == 0`.
    ///
    /// # Errors
    /// Returns an error for a non-positive/non-finite `dt_days`, a grid
    /// (forcing or radiation) that does not match the DEM shape, or a
    /// missing radiation grid with `srf > 0`.
    pub fn step_radiation(
        &mut self,
        forcing: &Forcing,
        radiation: Option<ArrayView2<'_, f64>>,
        dt_days: f64,
    ) -> Result<StepOutput> {
        if !dt_days.is_finite() || dt_days <= 0.0 {
            return Err(SnowmeltError::InvalidParameter {
                name: "dt_days",
                reason: format!("must be finite and > 0, got {dt_days}"),
            });
        }
        let shape = self.dem.shape();
        if (self.params.srf > 0.0 || self.params.energy_balance.is_some()) && radiation.is_none() {
            return Err(SnowmeltError::InvalidParameter {
                name: "radiation",
                reason: "srf > 0 or energy-balance mode requires a radiation grid".to_string(),
            });
        }
        if let Some(rad) = &radiation
            && rad.dim() != shape
        {
            return Err(SnowmeltError::ShapeMismatch {
                expected: shape,
                got: rad.dim(),
            });
        }
        let rad: ArrayView2<'_, f64> = match &radiation {
            Some(r) => r.view(),
            None => self.rad_zero.view(),
        };
        let mut snowfall = Array2::zeros(shape);
        let mut rain = Array2::zeros(shape);
        let mut melt = Array2::zeros(shape);
        let params = self.params;

        // Resolve per-cell temperature and precipitation views.
        let (temp, precip): (ArrayView2<'_, f64>, ArrayView2<'_, f64>) = match forcing {
            Forcing::Uniform {
                t_ref,
                z_ref,
                precip,
            } => {
                let (t_ref, z_ref, precip) = (*t_ref, *z_ref, *precip);
                let lapse = params.lapse_rate;
                let grad = params.precip_gradient;
                Zip::from(&mut self.temp_buf)
                    .and(&mut self.precip_buf)
                    .and(self.dem.elevation())
                    .par_for_each(|t, p, &z| {
                        *t = t_ref + lapse * (z - z_ref);
                        let scaled = precip * (1.0 + grad * (z - z_ref));
                        *p = if scaled.is_nan() {
                            f64::NAN
                        } else {
                            scaled.max(0.0)
                        };
                    });
                (self.temp_buf.view(), self.precip_buf.view())
            }
            Forcing::Distributed { temp, precip } => {
                for grid in [temp, precip] {
                    if grid.dim() != shape {
                        return Err(SnowmeltError::ShapeMismatch {
                            expected: shape,
                            got: grid.dim(),
                        });
                    }
                }
                (temp.view(), precip.view())
            }
        };

        // Pass 1: rain–snow partition (rain is independent of melt, so it
        // can feed the rain-on-snow term of the energy balance).
        Zip::from(&mut snowfall)
            .and(temp)
            .and(precip)
            .par_for_each(|s, &t, &p| *s = p * params.snow_fraction(t));
        Zip::from(&mut rain)
            .and(precip)
            .and(&snowfall)
            .par_for_each(|r, &p, &s| *r = p - s);

        // Pass 2: snow age and albedo (only in decay mode; otherwise
        // `albedo_buf` keeps the constant set at construction).
        if let Some(decay) = params.albedo_decay {
            Zip::from(&mut self.age)
                .and(&mut self.albedo_buf)
                .and(&snowfall)
                .par_for_each(|age, albedo, &s| {
                    if s.is_finite() {
                        *age = if s >= decay.refresh_swe_mm {
                            0.0
                        } else {
                            *age + dt_days
                        };
                        *albedo = decay.albedo(*age);
                    } else {
                        *age = f64::NAN;
                        *albedo = f64::NAN;
                    }
                });
        }

        // Pass 3: accumulation and melt.
        let mut sublimation = Array2::zeros(shape);
        match params.energy_balance {
            Some(eb) => {
                // 3a: `(total, latent)` energy fluxes per cell (W/m²).
                Zip::from(&mut self.energy_buf)
                    .and(temp)
                    .and(rad)
                    .and(&self.albedo_buf)
                    .and(&self.pressure_buf)
                    .and(&rain)
                    .par_for_each(|e, &t, &g, &albedo, &press, &r| {
                        *e = if t.is_finite()
                            && g.is_finite()
                            && albedo.is_finite()
                            && press.is_finite()
                            && r.is_finite()
                        {
                            let (q, q_e) =
                                energy::energy_fluxes(&eb, t, g, albedo, press, r, dt_days);
                            [q, q_e]
                        } else {
                            [f64::NAN, f64::NAN]
                        };
                    });
                // 3b: cold content, melt, and sublimation mass loss.
                Zip::from(&mut self.swe)
                    .and(&self.energy_buf)
                    .and(&snowfall)
                    .and(&mut self.cold_content)
                    .and(&mut melt)
                    .and(&mut sublimation)
                    .par_for_each(|swe, &[q, q_e], &s, cold, m, sub| {
                        if !swe.is_finite() || !q.is_finite() || !s.is_finite() {
                            *swe = f64::NAN;
                            *cold = f64::NAN;
                            *m = f64::NAN;
                            *sub = f64::NAN;
                            return;
                        }
                        *swe += s;
                        let potential = energy::apply_energy(&eb, q, dt_days, *swe, cold);
                        *m = potential.min(*swe);
                        *swe -= *m;
                        *sub = energy::sublimation_mm(q_e, dt_days).min(*swe);
                        *swe -= *sub;
                    });
            }
            None => {
                Zip::from(&mut self.swe)
                    .and(temp)
                    .and(rad)
                    .and(&self.albedo_buf)
                    .and(&snowfall)
                    .and(&mut melt)
                    .par_for_each(|swe, &t, &g, &albedo, &s, m| {
                        *m = melt_cell(&params, dt_days, swe, t, g, albedo, s);
                    });
                // Sin sublimación en modo índice de temperatura; propagar
                // nodata para consistencia con el resto de las salidas.
                Zip::from(&mut sublimation)
                    .and(&melt)
                    .par_for_each(|sub, &m| {
                        if !m.is_finite() {
                            *sub = f64::NAN;
                        }
                    });
            }
        }

        Ok(StepOutput {
            snowfall,
            rain,
            melt,
            sublimation,
        })
    }

    /// Advances the model by one step of `dt_days` days without radiation
    /// forcing (pure degree-day; requires `srf == 0`).
    ///
    /// # Errors
    /// See [`Self::step_radiation`].
    pub fn step_days(&mut self, forcing: &Forcing, dt_days: f64) -> Result<StepOutput> {
        self.step_radiation(forcing, None, dt_days)
    }

    /// Advances the model by one daily step.
    ///
    /// # Errors
    /// See [`Self::step_days`].
    pub fn step(&mut self, forcing: &Forcing) -> Result<StepOutput> {
        self.step_days(forcing, 1.0)
    }

    /// Runs a sequence of daily forcings, returning one summary per step.
    ///
    /// # Errors
    /// Stops at the first failing step (see [`Self::step_days`]).
    pub fn run<'a, I>(&mut self, forcings: I) -> Result<Vec<StepSummary>>
    where
        I: IntoIterator<Item = &'a Forcing>,
    {
        forcings
            .into_iter()
            .map(|f| {
                let out = self.step(f)?;
                Ok(self.summarize(&out))
            })
            .collect()
    }

    /// Aggregates a step's fluxes and the current SWE over valid cells.
    pub fn summarize(&self, out: &StepOutput) -> StepSummary {
        let mean_snowfall = nan_mean(&out.snowfall);
        let mean_rain = nan_mean(&out.rain);
        let mean_melt = nan_mean(&out.melt);
        let mean_sublimation = nan_mean(&out.sublimation);
        let mean_swe = nan_mean(&self.swe);
        let mean_albedo = nan_mean(&self.albedo_buf);
        let (covered, valid) = self
            .swe
            .iter()
            .filter(|s| s.is_finite())
            .fold((0usize, 0usize), |(c, n), &s| {
                (c + usize::from(s > SNOW_COVER_THRESHOLD_MM), n + 1)
            });
        let snow_cover_fraction = if valid == 0 {
            f64::NAN
        } else {
            covered as f64 / valid as f64
        };
        StepSummary {
            mean_snowfall,
            mean_rain,
            mean_melt,
            mean_sublimation,
            mean_runoff: mean_rain + mean_melt,
            mean_swe,
            mean_albedo,
            snow_cover_fraction,
        }
    }
}

/// Accumulates snowfall into `swe` and melts one cell; returns melt in mm.
///
/// Melt follows the enhanced temperature-index formulation (Pellicciotti
/// et al. 2005): for `T > t_melt`, `ddf·(T − t_melt) + srf·(1 − albedo)·G`,
/// which reduces to classic degree-day when `srf == 0`.
#[inline]
fn melt_cell(
    params: &DegreeDayParams,
    dt_days: f64,
    swe: &mut f64,
    t_c: f64,
    rad_wm2: f64,
    albedo: f64,
    snowfall_mm: f64,
) -> f64 {
    if !swe.is_finite()
        || !t_c.is_finite()
        || !rad_wm2.is_finite()
        || !albedo.is_finite()
        || !snowfall_mm.is_finite()
    {
        *swe = f64::NAN;
        return f64::NAN;
    }
    *swe += snowfall_mm;
    let potential = if t_c > params.t_melt {
        (params.ddf * (t_c - params.t_melt) + params.srf * (1.0 - albedo) * rad_wm2) * dt_days
    } else {
        0.0
    };
    let melt = potential.min(*swe);
    *swe -= melt;
    melt
}

/// Mean over finite cells; `NaN` if there are none.
fn nan_mean(grid: &Array2<f64>) -> f64 {
    let (sum, n) = grid
        .iter()
        .filter(|v| v.is_finite())
        .fold((0.0, 0usize), |(s, n), &v| (s + v, n + 1));
    if n == 0 { f64::NAN } else { sum / n as f64 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    fn flat_dem(rows: usize, cols: usize, z: f64) -> Dem {
        Dem::new(Array2::from_elem((rows, cols), z)).unwrap()
    }

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn cold_day_accumulates_all_precip_as_snow() {
        let mut m = SnowModel::new(flat_dem(2, 2, 1000.0), DegreeDayParams::default()).unwrap();
        let out = m
            .step(&Forcing::Uniform {
                t_ref: -5.0,
                z_ref: 1000.0,
                precip: 12.0,
            })
            .unwrap();
        assert!(out.snowfall.iter().all(|&s| approx(s, 12.0)));
        assert!(out.rain.iter().all(|&r| approx(r, 0.0)));
        assert!(out.melt.iter().all(|&x| approx(x, 0.0)));
        assert!(m.swe().iter().all(|&s| approx(s, 12.0)));
    }

    #[test]
    fn warm_day_without_snowpack_is_pure_rain() {
        let mut m = SnowModel::new(flat_dem(2, 2, 500.0), DegreeDayParams::default()).unwrap();
        let out = m
            .step(&Forcing::Uniform {
                t_ref: 10.0,
                z_ref: 500.0,
                precip: 8.0,
            })
            .unwrap();
        assert!(out.rain.iter().all(|&r| approx(r, 8.0)));
        assert!(out.melt.iter().all(|&x| approx(x, 0.0)));
        assert!(m.swe().iter().all(|&s| approx(s, 0.0)));
    }

    #[test]
    fn melt_is_capped_by_available_swe() {
        let dem = flat_dem(1, 1, 1000.0);
        let initial = array![[5.0]];
        let mut m = SnowModel::with_initial_swe(dem, DegreeDayParams::default(), initial).unwrap();
        // Potential melt = 4.0 * 10 = 40 mm >> 5 mm available.
        let out = m
            .step(&Forcing::Uniform {
                t_ref: 10.0,
                z_ref: 1000.0,
                precip: 0.0,
            })
            .unwrap();
        assert!(approx(out.melt[[0, 0]], 5.0));
        assert!(approx(m.swe()[[0, 0]], 0.0));
    }

    #[test]
    fn mixed_precipitation_at_threshold_midpoint() {
        let mut m = SnowModel::new(flat_dem(1, 1, 0.0), DegreeDayParams::default()).unwrap();
        // t = 1 °C is the midpoint of [0, 2] → half snow, half rain.
        let out = m
            .step(&Forcing::Uniform {
                t_ref: 1.0,
                z_ref: 0.0,
                precip: 10.0,
            })
            .unwrap();
        assert!(approx(out.snowfall[[0, 0]], 5.0));
        assert!(approx(out.rain[[0, 0]], 5.0));
    }

    #[test]
    fn lapse_rate_keeps_high_cells_colder() {
        // Two cells, 0 m and 2000 m; forcing measured at 0 m and +3 °C.
        let dem = Dem::new(array![[0.0, 2000.0]]).unwrap();
        let mut m = SnowModel::new(dem, DegreeDayParams::default()).unwrap();
        let out = m
            .step(&Forcing::Uniform {
                t_ref: 3.0,
                z_ref: 0.0,
                precip: 10.0,
            })
            .unwrap();
        // Low cell: 3 °C → all rain. High cell: 3 - 13 = -10 °C → all snow.
        assert!(approx(out.rain[[0, 0]], 10.0));
        assert!(approx(out.snowfall[[0, 1]], 10.0));
        assert!(m.swe()[[0, 1]] > m.swe()[[0, 0]]);
    }

    #[test]
    fn nodata_cells_propagate_nan() {
        let dem = Dem::new(array![[1000.0, f64::NAN]]).unwrap();
        let mut m = SnowModel::new(dem, DegreeDayParams::default()).unwrap();
        let out = m
            .step(&Forcing::Uniform {
                t_ref: -1.0,
                z_ref: 1000.0,
                precip: 5.0,
            })
            .unwrap();
        assert!(out.snowfall[[0, 1]].is_nan());
        assert!(out.melt[[0, 1]].is_nan());
        assert!(m.swe()[[0, 1]].is_nan());
        assert!(approx(m.swe()[[0, 0]], 5.0));
        let s = m.summarize(&out);
        assert!(approx(s.mean_snowfall, 5.0));
        assert!(approx(s.snow_cover_fraction, 1.0));
    }

    #[test]
    fn mass_balance_holds_over_a_sequence() {
        let dem = Dem::new(array![[0.0, 1500.0], [3000.0, f64::NAN]]).unwrap();
        let mut m = SnowModel::new(dem, DegreeDayParams::default()).unwrap();
        let series = [
            (-6.0, 10.0),
            (-2.0, 0.0),
            (1.0, 14.0),
            (4.0, 3.0),
            (8.0, 0.0),
            (-1.0, 7.0),
            (12.0, 5.0),
        ];
        let mut total_precip = 0.0;
        let mut total_rain = Array2::zeros((2, 2));
        let mut total_melt = Array2::zeros((2, 2));
        for (t, p) in series {
            let out = m
                .step(&Forcing::Uniform {
                    t_ref: t,
                    z_ref: 0.0,
                    precip: p,
                })
                .unwrap();
            total_precip += p;
            total_rain += &out.rain;
            total_melt += &out.melt;
        }
        // Per valid cell: precip = Δswe + rain + melt (initial SWE was 0).
        for idx in [(0, 0), (0, 1), (1, 0)] {
            let balance = m.swe()[idx] + total_rain[idx] + total_melt[idx];
            assert!(
                approx(balance, total_precip),
                "cell {idx:?}: {balance} != {total_precip}"
            );
        }
    }

    #[test]
    fn distributed_forcing_rejects_wrong_shape() {
        let mut m = SnowModel::new(flat_dem(2, 2, 100.0), DegreeDayParams::default()).unwrap();
        let err = m
            .step(&Forcing::Distributed {
                temp: Array2::zeros((3, 3)),
                precip: Array2::zeros((2, 2)),
            })
            .unwrap_err();
        assert!(matches!(err, SnowmeltError::ShapeMismatch { .. }));
    }

    #[test]
    fn distributed_matches_uniform_when_built_by_lapse() {
        // El CLI construye el forzante distribuido extrapolando por lapse;
        // debe dar el mismo resultado que el forzante uniforme equivalente.
        let dem = Dem::new(array![[1000.0, 2500.0], [4000.0, 1750.0]]).unwrap();
        let params = DegreeDayParams::default();
        let (t_ref, z_ref, precip) = (2.0, 2000.0, 8.0);

        let mut uni = SnowModel::new(dem.clone(), params).unwrap();
        let out_u = uni
            .step(&Forcing::Uniform {
                t_ref,
                z_ref,
                precip,
            })
            .unwrap();

        let lapse = params.lapse_rate;
        let temp = dem.elevation().mapv(|z| t_ref + lapse * (z - z_ref));
        let precip_grid = dem.elevation().mapv(|_| precip);
        let mut dist = SnowModel::new(dem, params).unwrap();
        let out_d = dist
            .step(&Forcing::Distributed {
                temp,
                precip: precip_grid,
            })
            .unwrap();

        for idx in [(0, 0), (0, 1), (1, 0), (1, 1)] {
            assert!(approx(out_u.snowfall[idx], out_d.snowfall[idx]));
            assert!(approx(out_u.melt[idx], out_d.melt[idx]));
            assert!(approx(uni.swe()[idx], dist.swe()[idx]));
        }
    }

    #[test]
    fn distributed_forcing_drives_each_cell_independently() {
        let mut m = SnowModel::new(flat_dem(1, 2, 0.0), DegreeDayParams::default()).unwrap();
        let out = m
            .step(&Forcing::Distributed {
                temp: array![[-5.0, 5.0]],
                precip: array![[10.0, 10.0]],
            })
            .unwrap();
        assert!(approx(out.snowfall[[0, 0]], 10.0));
        assert!(approx(out.rain[[0, 1]], 10.0));
    }

    #[test]
    fn rejects_invalid_dt_and_negative_initial_swe() {
        let mut m = SnowModel::new(flat_dem(1, 1, 0.0), DegreeDayParams::default()).unwrap();
        let f = Forcing::Uniform {
            t_ref: 0.0,
            z_ref: 0.0,
            precip: 0.0,
        };
        assert!(m.step_days(&f, 0.0).is_err());
        assert!(m.step_days(&f, f64::NAN).is_err());

        let err = SnowModel::with_initial_swe(
            flat_dem(1, 1, 0.0),
            DegreeDayParams::default(),
            array![[-1.0]],
        )
        .unwrap_err();
        assert!(matches!(err, SnowmeltError::InvalidParameter { .. }));
    }

    #[test]
    fn eti_radiation_term_adds_melt() {
        let params = DegreeDayParams {
            srf: 0.2,
            albedo: 0.5,
            ..DegreeDayParams::default()
        };
        let mut m =
            SnowModel::with_initial_swe(flat_dem(1, 2, 0.0), params, array![[500.0, 500.0]])
                .unwrap();
        // Same temperature, different radiation: 0 vs 200 W/m².
        let out = m
            .step_radiation(
                &Forcing::Uniform {
                    t_ref: 5.0,
                    z_ref: 0.0,
                    precip: 0.0,
                },
                Some(array![[0.0, 200.0]].view()),
                1.0,
            )
            .unwrap();
        // Cell 0: pure degree-day = 4*5 = 20. Cell 1: + 0.2*0.5*200 = +20.
        assert!(approx(out.melt[[0, 0]], 20.0));
        assert!(approx(out.melt[[0, 1]], 40.0));
    }

    #[test]
    fn radiation_ignored_below_melt_threshold() {
        let params = DegreeDayParams {
            srf: 0.2,
            ..DegreeDayParams::default()
        };
        let mut m =
            SnowModel::with_initial_swe(flat_dem(1, 1, 0.0), params, array![[100.0]]).unwrap();
        let out = m
            .step_radiation(
                &Forcing::Uniform {
                    t_ref: -3.0,
                    z_ref: 0.0,
                    precip: 0.0,
                },
                Some(array![[300.0]].view()),
                1.0,
            )
            .unwrap();
        assert!(approx(out.melt[[0, 0]], 0.0));
    }

    #[test]
    fn srf_without_radiation_grid_is_an_error() {
        let params = DegreeDayParams {
            srf: 0.2,
            ..DegreeDayParams::default()
        };
        let mut m = SnowModel::new(flat_dem(1, 1, 0.0), params).unwrap();
        let err = m
            .step(&Forcing::Uniform {
                t_ref: 5.0,
                z_ref: 0.0,
                precip: 0.0,
            })
            .unwrap_err();
        assert!(matches!(err, SnowmeltError::InvalidParameter { .. }));
    }

    #[test]
    fn radiation_grid_shape_is_checked() {
        let mut m = SnowModel::new(flat_dem(2, 2, 0.0), DegreeDayParams::default()).unwrap();
        let err = m
            .step_radiation(
                &Forcing::Uniform {
                    t_ref: 0.0,
                    z_ref: 0.0,
                    precip: 0.0,
                },
                Some(Array2::zeros((3, 3)).view()),
                1.0,
            )
            .unwrap_err();
        assert!(matches!(err, SnowmeltError::ShapeMismatch { .. }));
    }

    #[test]
    fn precip_gradient_scales_with_elevation() {
        let params = DegreeDayParams {
            precip_gradient: 0.0005,
            ..DegreeDayParams::default()
        };
        // Cells at 0 m and 2000 m; forcing at 0 m, 10 mm, cold (all snow).
        let dem = Dem::new(array![[0.0, 2000.0]]).unwrap();
        let mut m = SnowModel::new(dem, params).unwrap();
        let out = m
            .step(&Forcing::Uniform {
                t_ref: -20.0,
                z_ref: 0.0,
                precip: 10.0,
            })
            .unwrap();
        assert!(approx(out.snowfall[[0, 0]], 10.0));
        // p(2000) = 10 * (1 + 0.0005*2000) = 20 mm.
        assert!(approx(out.snowfall[[0, 1]], 20.0));
    }

    #[test]
    fn precip_gradient_clamps_to_zero() {
        let params = DegreeDayParams {
            precip_gradient: -0.001,
            ..DegreeDayParams::default()
        };
        // At 2000 m: 10 * (1 - 0.001*2000) = -10 → clamped to 0.
        let dem = Dem::new(array![[2000.0]]).unwrap();
        let mut m = SnowModel::new(dem, params).unwrap();
        let out = m
            .step(&Forcing::Uniform {
                t_ref: -20.0,
                z_ref: 0.0,
                precip: 10.0,
            })
            .unwrap();
        assert!(approx(out.snowfall[[0, 0]], 0.0));
        assert!(approx(out.rain[[0, 0]], 0.0));
    }

    #[test]
    fn albedo_decay_increases_melt_over_dry_days() {
        use crate::params::AlbedoDecay;
        let params = DegreeDayParams {
            ddf: 0.0, // aislar el término radiativo
            srf: 0.2,
            albedo_decay: Some(AlbedoDecay::default()),
            ..DegreeDayParams::default()
        };
        let mut m =
            SnowModel::with_initial_swe(flat_dem(1, 1, 0.0), params, array![[1000.0]]).unwrap();
        let warm_dry = Forcing::Uniform {
            t_ref: 5.0,
            z_ref: 0.0,
            precip: 0.0,
        };
        let rad = array![[200.0]];
        let mut melts = Vec::new();
        for _ in 0..4 {
            let out = m.step_radiation(&warm_dry, Some(rad.view()), 1.0).unwrap();
            melts.push(out.melt[[0, 0]]);
        }
        // El albedo decae día a día → (1 − α) crece → más melt cada día.
        assert!(
            melts.windows(2).all(|w| w[1] > w[0]),
            "melt no monotónico: {melts:?}"
        );
        // Día 1: edad 1, α = 0.4 + 0.45·exp(−1/6); melt = 0.2·(1−α)·200.
        let alpha_1 = 0.4 + 0.45 * (-1.0_f64 / 6.0).exp();
        assert!((melts[0] - 0.2 * (1.0 - alpha_1) * 200.0).abs() < 1e-9);
    }

    #[test]
    fn snowfall_resets_albedo_to_fresh() {
        use crate::params::AlbedoDecay;
        let params = DegreeDayParams {
            srf: 0.2,
            albedo_decay: Some(AlbedoDecay::default()),
            ..DegreeDayParams::default()
        };
        let mut m =
            SnowModel::with_initial_swe(flat_dem(1, 1, 0.0), params, array![[1000.0]]).unwrap();
        let rad = array![[150.0]];
        // Envejecer la nieve 10 días secos y cálidos.
        for _ in 0..10 {
            let f = Forcing::Uniform {
                t_ref: 3.0,
                z_ref: 0.0,
                precip: 0.0,
            };
            m.step_radiation(&f, Some(rad.view()), 1.0).unwrap();
        }
        let aged = m.albedo()[[0, 0]];
        // Nevada fría que supera el umbral de refresco (1 mm).
        let snow_day = Forcing::Uniform {
            t_ref: -5.0,
            z_ref: 0.0,
            precip: 10.0,
        };
        m.step_radiation(&snow_day, Some(rad.view()), 1.0).unwrap();
        let fresh = m.albedo()[[0, 0]];
        assert!(aged < 0.6, "albedo envejecido: {aged}");
        assert_eq!(fresh, 0.85, "tras nevada debe volver a fresh");
        assert_eq!(m.snow_age()[[0, 0]], 0.0);
    }

    #[test]
    fn constant_albedo_behavior_unchanged_without_decay() {
        let params = DegreeDayParams {
            srf: 0.2,
            albedo: 0.5,
            ..DegreeDayParams::default()
        };
        let mut m =
            SnowModel::with_initial_swe(flat_dem(1, 1, 0.0), params, array![[500.0]]).unwrap();
        let f = Forcing::Uniform {
            t_ref: 5.0,
            z_ref: 0.0,
            precip: 0.0,
        };
        let rad = array![[200.0]];
        // Mismo melt en pasos sucesivos: 4·5 + 0.2·0.5·200 = 40.
        for _ in 0..3 {
            let out = m.step_radiation(&f, Some(rad.view()), 1.0).unwrap();
            assert!(approx(out.melt[[0, 0]], 40.0));
        }
    }

    #[test]
    fn albedo_decay_params_are_validated() {
        use crate::params::AlbedoDecay;
        let bad_tau = DegreeDayParams {
            albedo_decay: Some(AlbedoDecay {
                tau_days: 0.0,
                ..AlbedoDecay::default()
            }),
            ..DegreeDayParams::default()
        };
        assert!(SnowModel::new(flat_dem(1, 1, 0.0), bad_tau).is_err());

        let crossed = DegreeDayParams {
            albedo_decay: Some(AlbedoDecay {
                albedo_fresh: 0.4,
                albedo_min: 0.8,
                ..AlbedoDecay::default()
            }),
            ..DegreeDayParams::default()
        };
        assert!(SnowModel::new(flat_dem(1, 1, 0.0), crossed).is_err());
    }

    #[test]
    fn energy_balance_melts_on_warm_sunny_day_only() {
        use crate::energy::EnergyBalanceParams;
        let params = DegreeDayParams {
            energy_balance: Some(EnergyBalanceParams::default()),
            ..DegreeDayParams::default()
        };
        let mut m =
            SnowModel::with_initial_swe(flat_dem(1, 1, 2000.0), params, array![[500.0]]).unwrap();
        // Día frío: nada de melt, crece el cold content (el enfriamiento
        // longwave de cielo despejado es fuerte: ~6 MJ/m² en un día a -5).
        let cold_day = Forcing::Uniform {
            t_ref: -5.0,
            z_ref: 2000.0,
            precip: 0.0,
        };
        let out = m
            .step_radiation(&cold_day, Some(array![[80.0]].view()), 1.0)
            .unwrap();
        assert_eq!(out.melt[[0, 0]], 0.0);
        assert!(m.cold_content()[[0, 0]] > 0.0);

        // Día cálido y soleado (~12 MJ/m²): paga el frío y derrite.
        let warm_day = Forcing::Uniform {
            t_ref: 10.0,
            z_ref: 2000.0,
            precip: 0.0,
        };
        let out = m
            .step_radiation(&warm_day, Some(array![[350.0]].view()), 1.0)
            .unwrap();
        let melt = out.melt[[0, 0]];
        assert!(melt > 5.0 && melt < 80.0, "melt = {melt}");
        assert_eq!(m.cold_content()[[0, 0]], 0.0);
    }

    #[test]
    fn energy_balance_cold_content_delays_melt() {
        use crate::energy::EnergyBalanceParams;
        let params = DegreeDayParams {
            energy_balance: Some(EnergyBalanceParams::default()),
            ..DegreeDayParams::default()
        };
        let mut frio =
            SnowModel::with_initial_swe(flat_dem(1, 1, 2000.0), params, array![[500.0]]).unwrap();
        let mut tibio =
            SnowModel::with_initial_swe(flat_dem(1, 1, 2000.0), params, array![[500.0]]).unwrap();
        // Enfriar uno de los packs por 5 días.
        let cold = Forcing::Uniform {
            t_ref: -15.0,
            z_ref: 2000.0,
            precip: 0.0,
        };
        for _ in 0..5 {
            frio.step_radiation(&cold, Some(array![[20.0]].view()), 1.0)
                .unwrap();
        }
        // Mismo día cálido para ambos.
        let warm = Forcing::Uniform {
            t_ref: 6.0,
            z_ref: 2000.0,
            precip: 0.0,
        };
        let melt_frio = frio
            .step_radiation(&warm, Some(array![[250.0]].view()), 1.0)
            .unwrap()
            .melt[[0, 0]];
        let melt_tibio = tibio
            .step_radiation(&warm, Some(array![[250.0]].view()), 1.0)
            .unwrap()
            .melt[[0, 0]];
        assert!(
            melt_frio < melt_tibio,
            "pack frío ({melt_frio}) debe derretir menos que el isotérmico ({melt_tibio})"
        );
    }

    #[test]
    fn energy_balance_requires_radiation_and_respects_mass_balance() {
        use crate::energy::EnergyBalanceParams;
        let params = DegreeDayParams {
            energy_balance: Some(EnergyBalanceParams::default()),
            ..DegreeDayParams::default()
        };
        let mut m =
            SnowModel::with_initial_swe(flat_dem(1, 1, 1000.0), params, array![[3.0]]).unwrap();
        // Sin radiación → error.
        assert!(
            m.step(&Forcing::Uniform {
                t_ref: 5.0,
                z_ref: 1000.0,
                precip: 0.0
            })
            .is_err()
        );
        // Día tórrido: el melt no puede exceder el SWE disponible.
        let out = m
            .step_radiation(
                &Forcing::Uniform {
                    t_ref: 20.0,
                    z_ref: 1000.0,
                    precip: 0.0,
                },
                Some(array![[400.0]].view()),
                1.0,
            )
            .unwrap();
        assert!(approx(out.melt[[0, 0]], 3.0));
        assert!(approx(m.swe()[[0, 0]], 0.0));
    }

    #[test]
    fn energy_balance_sublimation_removes_mass_and_balance_closes() {
        use crate::energy::EnergyBalanceParams;
        let params = DegreeDayParams {
            energy_balance: Some(EnergyBalanceParams {
                wind_speed: 5.0,
                rel_humidity: 0.2, // aire seco → sublimación fuerte
                ..Default::default()
            }),
            ..DegreeDayParams::default()
        };
        let mut m =
            SnowModel::with_initial_swe(flat_dem(1, 1, 4000.0), params, array![[300.0]]).unwrap();
        let cold_dry = Forcing::Uniform {
            t_ref: -10.0,
            z_ref: 4000.0,
            precip: 0.0,
        };
        let mut total_subl = 0.0;
        let mut total_melt = 0.0;
        for _ in 0..10 {
            let out = m
                .step_radiation(&cold_dry, Some(array![[100.0]].view()), 1.0)
                .unwrap();
            total_subl += out.sublimation[[0, 0]];
            total_melt += out.melt[[0, 0]];
        }
        assert!(total_subl > 1.0, "sublimación = {total_subl}");
        // Balance de masa: SWE_0 = SWE_f + melt + sublimación (precip 0).
        let balance = m.swe()[[0, 0]] + total_melt + total_subl;
        assert!((balance - 300.0).abs() < 1e-9, "balance = {balance}");
    }

    #[test]
    fn energy_balance_propagates_nodata() {
        use crate::energy::EnergyBalanceParams;
        let params = DegreeDayParams {
            energy_balance: Some(EnergyBalanceParams::default()),
            ..DegreeDayParams::default()
        };
        let dem = Dem::new(array![[1000.0, f64::NAN]]).unwrap();
        let mut m = SnowModel::new(dem, params).unwrap();
        let out = m
            .step_radiation(
                &Forcing::Uniform {
                    t_ref: -3.0,
                    z_ref: 1000.0,
                    precip: 10.0,
                },
                Some(array![[100.0, 100.0]].view()),
                1.0,
            )
            .unwrap();
        assert!(out.melt[[0, 1]].is_nan());
        assert!(m.cold_content()[[0, 1]].is_nan());
        assert!(m.swe()[[0, 0]] > 0.0);
    }

    #[test]
    fn runoff_is_rain_plus_melt() {
        let mut m = SnowModel::with_initial_swe(
            flat_dem(1, 1, 0.0),
            DegreeDayParams::default(),
            array![[100.0]],
        )
        .unwrap();
        let out = m
            .step(&Forcing::Uniform {
                t_ref: 5.0,
                z_ref: 0.0,
                precip: 6.0,
            })
            .unwrap();
        let runoff = out.runoff();
        assert!(approx(runoff[[0, 0]], out.rain[[0, 0]] + out.melt[[0, 0]]));
        // 5 °C: all rain (6 mm) + melt 4*5 = 20 mm.
        assert!(approx(runoff[[0, 0]], 26.0));
    }
}
