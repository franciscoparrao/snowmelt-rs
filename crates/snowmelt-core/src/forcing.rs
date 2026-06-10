//! Meteorological forcings for one model step.

use ndarray::Array2;

/// Forcing for a single time step.
///
/// All precipitation is in mm of water equivalent per step; temperatures
/// in °C.
#[derive(Debug, Clone)]
pub enum Forcing {
    /// A single station / basin-mean value, distributed over the DEM.
    ///
    /// Temperature is extrapolated to each cell with the model's lapse rate:
    /// `t(z) = t_ref + lapse_rate * (z - z_ref)`. Precipitation is applied
    /// uniformly.
    Uniform {
        /// Air temperature at the reference elevation (°C).
        t_ref: f64,
        /// Elevation at which `t_ref` was measured (m).
        z_ref: f64,
        /// Precipitation during the step (mm).
        precip: f64,
    },
    /// Fully distributed grids, e.g. from a reanalysis or interpolation.
    ///
    /// Both grids must match the DEM shape.
    Distributed {
        /// Per-cell air temperature (°C).
        temp: Array2<f64>,
        /// Per-cell precipitation (mm).
        precip: Array2<f64>,
    },
}
