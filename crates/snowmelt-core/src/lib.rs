//! # snowmelt-core
//!
//! Distributed snowpack model over a digital elevation model (DEM).
//!
//! Implements a degree-day (temperature-index) snow model: per-cell snow
//! water equivalent (SWE) accumulation and ablation, with rain–snow
//! partitioning by air temperature and lapse-rate extrapolation of a
//! reference temperature over the DEM.
//!
//! This crate performs **no I/O**: callers provide [`ndarray`] grids and
//! per-step forcings. Cells with `NaN` elevation are treated as nodata and
//! propagate `NaN` through all state and outputs. Per-cell computations run
//! in parallel via Rayon.
//!
//! ## Example
//!
//! ```
//! use ndarray::Array2;
//! use snowmelt_core::{Dem, DegreeDayParams, Forcing, SnowModel};
//!
//! let dem = Dem::new(Array2::from_shape_fn((3, 3), |(i, j)| {
//!     1000.0 + 500.0 * (i + j) as f64
//! }))?;
//! let mut model = SnowModel::new(dem, DegreeDayParams::default())?;
//!
//! // One cold, snowy day: temperature measured at 1000 m.
//! let out = model.step(&Forcing::Uniform { t_ref: -3.0, z_ref: 1000.0, precip: 10.0 })?;
//! assert!(out.melt.iter().all(|&m| m == 0.0));
//! # Ok::<(), snowmelt_core::SnowmeltError>(())
//! ```

pub mod dem;
pub mod downscale;
pub mod energy;
pub mod error;
pub mod forcing;
pub mod model;
pub mod params;
pub mod routing;
pub mod terrain;

pub use dem::Dem;
pub use downscale::{DownscaleParams, Downscaler};
pub use energy::EnergyBalanceParams;
pub use error::{Result, SnowmeltError};
pub use forcing::Forcing;
pub use model::{SNOW_COVER_THRESHOLD_MM, SnowModel, StepOutput, StepSummary};
pub use params::{AlbedoDecay, DegreeDayParams};
pub use routing::{LinearReservoir, depth_to_discharge};
