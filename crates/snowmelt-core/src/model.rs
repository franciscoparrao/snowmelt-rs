//! Snowpack state and time integration.

use ndarray::{Array2, ArrayView2, Zip};

use crate::dem::Dem;
use crate::error::{Result, SnowmeltError};
use crate::forcing::Forcing;
use crate::params::DegreeDayParams;

/// SWE (mm) above which a cell counts as snow-covered in
/// [`SnowModel::summarize`]. Filters out floating-point melt residues.
pub const SNOW_COVER_THRESHOLD_MM: f64 = 0.1;

/// Per-cell fluxes produced by one model step (all in mm w.e.).
///
/// Mass balance per cell and step: `snowfall + rain == precip` and
/// `Δswe == snowfall - melt`. Nodata cells are `NaN`.
#[derive(Debug)]
pub struct StepOutput {
    /// Solid precipitation added to the snowpack (mm).
    pub snowfall: Array2<f64>,
    /// Liquid precipitation, passed through to runoff (mm).
    pub rain: Array2<f64>,
    /// Snowmelt released from the snowpack (mm).
    pub melt: Array2<f64>,
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
    /// Mean runoff = rain + melt (mm).
    pub mean_runoff: f64,
    /// Mean SWE after the step (mm).
    pub mean_swe: f64,
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
    temp_buf: Array2<f64>,
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
        let temp_buf = Array2::zeros(dem.shape());
        Ok(Self {
            dem,
            params,
            swe,
            temp_buf,
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

    /// Advances the model by one step of `dt_days` days.
    ///
    /// # Errors
    /// Returns an error for a non-positive/non-finite `dt_days` or a
    /// distributed forcing whose grids do not match the DEM shape.
    pub fn step_days(&mut self, forcing: &Forcing, dt_days: f64) -> Result<StepOutput> {
        if !dt_days.is_finite() || dt_days <= 0.0 {
            return Err(SnowmeltError::InvalidParameter {
                name: "dt_days",
                reason: format!("must be finite and > 0, got {dt_days}"),
            });
        }
        let shape = self.dem.shape();
        let mut snowfall = Array2::zeros(shape);
        let mut rain = Array2::zeros(shape);
        let mut melt = Array2::zeros(shape);
        let params = self.params;

        match forcing {
            Forcing::Uniform {
                t_ref,
                z_ref,
                precip,
            } => {
                let (t_ref, z_ref, precip) = (*t_ref, *z_ref, *precip);
                let lapse = params.lapse_rate;
                Zip::from(&mut self.temp_buf)
                    .and(self.dem.elevation())
                    .par_for_each(|t, &z| *t = t_ref + lapse * (z - z_ref));
                Zip::from(&mut self.swe)
                    .and(&self.temp_buf)
                    .and(&mut snowfall)
                    .and(&mut rain)
                    .and(&mut melt)
                    .par_for_each(|swe, &t, s, r, m| {
                        (*s, *r, *m) = cell_step(&params, dt_days, swe, t, precip);
                    });
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
                Zip::from(&mut self.swe)
                    .and(temp)
                    .and(precip)
                    .and(&mut snowfall)
                    .and(&mut rain)
                    .and(&mut melt)
                    .par_for_each(|swe, &t, &p, s, r, m| {
                        (*s, *r, *m) = cell_step(&params, dt_days, swe, t, p);
                    });
            }
        }

        Ok(StepOutput {
            snowfall,
            rain,
            melt,
        })
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
        let mean_swe = nan_mean(&self.swe);
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
            mean_runoff: mean_rain + mean_melt,
            mean_swe,
            snow_cover_fraction,
        }
    }
}

/// Updates one cell for one step; returns `(snowfall, rain, melt)` in mm.
#[inline]
fn cell_step(
    params: &DegreeDayParams,
    dt_days: f64,
    swe: &mut f64,
    t_c: f64,
    precip_mm: f64,
) -> (f64, f64, f64) {
    if !swe.is_finite() || !t_c.is_finite() || !precip_mm.is_finite() {
        *swe = f64::NAN;
        return (f64::NAN, f64::NAN, f64::NAN);
    }
    let snow_frac = params.snow_fraction(t_c);
    let snowfall = precip_mm * snow_frac;
    let rain = precip_mm - snowfall;
    *swe += snowfall;
    let potential = params.ddf * (t_c - params.t_melt).max(0.0) * dt_days;
    let melt = potential.min(*swe);
    *swe -= melt;
    (snowfall, rain, melt)
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
