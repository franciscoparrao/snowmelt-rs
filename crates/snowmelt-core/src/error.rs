//! Error types for the snowmelt model.

use thiserror::Error;

/// Errors produced by `snowmelt-core`.
#[derive(Debug, Error)]
pub enum SnowmeltError {
    /// A grid does not match the DEM shape.
    #[error("shape mismatch: expected {expected:?}, got {got:?}")]
    ShapeMismatch {
        /// Expected `(rows, cols)`.
        expected: (usize, usize),
        /// Actual `(rows, cols)`.
        got: (usize, usize),
    },

    /// A parameter or input value is out of its valid domain.
    #[error("invalid parameter `{name}`: {reason}")]
    InvalidParameter {
        /// Parameter name.
        name: &'static str,
        /// Why it is invalid.
        reason: String,
    },

    /// The DEM has zero rows or columns.
    #[error("empty grid: the DEM must have at least one cell")]
    EmptyGrid,
}

/// Convenience alias for results in this crate.
pub type Result<T> = std::result::Result<T, SnowmeltError>;
