//! Digital elevation model wrapper.

use ndarray::{Array2, ArrayView2};

use crate::error::{Result, SnowmeltError};

/// A digital elevation model: per-cell elevation in metres.
///
/// Cells with non-finite elevation (`NaN`, ±∞) are nodata; the model keeps
/// them as `NaN` in every state and output grid.
#[derive(Debug, Clone)]
pub struct Dem {
    elevation: Array2<f64>,
}

impl Dem {
    /// Wraps an elevation grid. Non-finite cells are normalised to `NaN`.
    ///
    /// # Errors
    /// Returns [`SnowmeltError::EmptyGrid`] if the grid has zero rows or columns.
    pub fn new(mut elevation: Array2<f64>) -> Result<Self> {
        let (rows, cols) = elevation.dim();
        if rows == 0 || cols == 0 {
            return Err(SnowmeltError::EmptyGrid);
        }
        elevation.mapv_inplace(|z| if z.is_finite() { z } else { f64::NAN });
        Ok(Self { elevation })
    }

    /// Grid shape as `(rows, cols)`.
    pub fn shape(&self) -> (usize, usize) {
        self.elevation.dim()
    }

    /// Elevation grid view (metres; `NaN` = nodata).
    pub fn elevation(&self) -> ArrayView2<'_, f64> {
        self.elevation.view()
    }

    /// Number of valid (non-nodata) cells.
    pub fn valid_cells(&self) -> usize {
        self.elevation.iter().filter(|z| z.is_finite()).count()
    }

    /// Mean elevation over valid cells (metres). `NaN` if no cell is valid.
    pub fn mean_elevation(&self) -> f64 {
        let (sum, n) = self
            .elevation
            .iter()
            .filter(|z| z.is_finite())
            .fold((0.0, 0usize), |(s, n), &z| (s + z, n + 1));
        if n == 0 { f64::NAN } else { sum / n as f64 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn rejects_empty_grid() {
        let err = Dem::new(Array2::zeros((0, 4))).unwrap_err();
        assert!(matches!(err, SnowmeltError::EmptyGrid));
    }

    #[test]
    fn normalises_non_finite_to_nan() {
        let dem = Dem::new(array![[1000.0, f64::INFINITY], [f64::NAN, 2000.0]]).unwrap();
        assert_eq!(dem.valid_cells(), 2);
        assert!(dem.elevation()[[0, 1]].is_nan());
    }

    #[test]
    fn mean_elevation_ignores_nodata() {
        let dem = Dem::new(array![[1000.0, f64::NAN], [3000.0, f64::NAN]]).unwrap();
        assert_eq!(dem.mean_elevation(), 2000.0);
    }
}
