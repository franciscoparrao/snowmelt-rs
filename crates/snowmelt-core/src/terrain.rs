//! Terrain derivatives for topographic downscaling.
//!
//! Self-contained slope, aspect and curvature from a DEM grid — this crate
//! deliberately does **not** depend on a GIS terrain library, so the
//! downscaler (and its tests) stand alone. Slope and aspect follow Horn's
//! (1981) 3×3 finite-difference operator; curvature follows the normalised
//! topographic-position measure of Liston & Elder (2006, *MicroMet*).
//!
//! Conventions:
//! - Rows increase **southward** (row 0 is the northern edge), as in an
//!   ESRI ASCII grid; columns increase **eastward**.
//! - `aspect` is the **downslope azimuth**, radians clockwise from north
//!   (0 = downhill toward north, π/2 = toward east). Flat cells get 0.
//! - `curvature` is dimensionless in `[-0.5, 0.5]`: positive on convex
//!   terrain (ridges, peaks), negative in concave terrain (valleys, basins).
//!
//! Cells on the border, or whose 3×3 neighbourhood touches a `NaN`, are
//! reported as flat (slope 0, aspect 0, curvature 0): the downscaler then
//! falls back to the unmodified forcing there.

use ndarray::Array2;

/// Slope (radians) and downslope aspect (radians clockwise from north),
/// by Horn's (1981) method. `cellsize` is the grid spacing in metres.
///
/// Returns two grids shaped like `elevation`. Border cells and cells with
/// any `NaN` neighbour are flat (slope 0, aspect 0); `NaN` elevation cells
/// stay `NaN` in the slope grid.
pub fn slope_aspect(elevation: &Array2<f64>, cellsize: f64) -> (Array2<f64>, Array2<f64>) {
    let (rows, cols) = elevation.dim();
    let mut slope = Array2::zeros((rows, cols));
    let mut aspect = Array2::zeros((rows, cols));
    for i in 0..rows {
        for j in 0..cols {
            let z0 = elevation[[i, j]];
            if !z0.is_finite() {
                slope[[i, j]] = f64::NAN;
                aspect[[i, j]] = f64::NAN;
                continue;
            }
            let Some(nb) = neighbourhood(elevation, i, j) else {
                continue; // border or NaN neighbour → flat
            };
            // Horn: dz/dx positive eastward, dz/dy positive northward.
            let dz_dx =
                ((nb.ne + 2.0 * nb.e + nb.se) - (nb.nw + 2.0 * nb.w + nb.sw)) / (8.0 * cellsize);
            let dz_dy =
                ((nb.nw + 2.0 * nb.n + nb.ne) - (nb.sw + 2.0 * nb.s + nb.se)) / (8.0 * cellsize);
            slope[[i, j]] = dz_dx.hypot(dz_dy).atan();
            if dz_dx == 0.0 && dz_dy == 0.0 {
                aspect[[i, j]] = 0.0;
            } else {
                // Downslope direction is −gradient; azimuth clockwise from
                // north = atan2(east, north) of that vector.
                let az = (-dz_dx).atan2(-dz_dy);
                aspect[[i, j]] = az.rem_euclid(std::f64::consts::TAU);
            }
        }
    }
    (slope, aspect)
}

/// Normalised topographic curvature in `[-0.5, 0.5]` (Liston & Elder 2006).
///
/// For each cell, the mean of its height above the four pairs of opposite
/// neighbours (W–E, S–N and the two diagonals), then scaled by twice the
/// domain maximum absolute value. Positive is convex (ridge), negative is
/// concave (valley). The absolute grid spacing cancels in the
/// normalisation, so only the DEM shape matters. Border/`NaN`-adjacent
/// cells are 0; `NaN` elevation cells stay `NaN`.
pub fn curvature(elevation: &Array2<f64>) -> Array2<f64> {
    let (rows, cols) = elevation.dim();
    let mut raw = Array2::zeros((rows, cols));
    let mut valid = Array2::from_elem((rows, cols), false);
    let inv_diag = 1.0 / std::f64::consts::SQRT_2;
    let mut max_abs = 0.0_f64;
    for i in 0..rows {
        for j in 0..cols {
            if !elevation[[i, j]].is_finite() {
                raw[[i, j]] = f64::NAN;
                continue;
            }
            let Some(nb) = neighbourhood(elevation, i, j) else {
                continue; // stays 0, valid stays false
            };
            let z0 = nb.z0;
            // Orthogonal pairs at unit spacing, diagonals at √2 spacing.
            let c_we = z0 - 0.5 * (nb.w + nb.e);
            let c_sn = z0 - 0.5 * (nb.s + nb.n);
            let c_d1 = (z0 - 0.5 * (nb.sw + nb.ne)) * inv_diag;
            let c_d2 = (z0 - 0.5 * (nb.se + nb.nw)) * inv_diag;
            let value = 0.25 * (c_we + c_sn + c_d1 + c_d2);
            raw[[i, j]] = value;
            valid[[i, j]] = true;
            max_abs = max_abs.max(value.abs());
        }
    }
    if max_abs == 0.0 {
        // Flat domain: everything that is valid is 0; NaN stays NaN.
        return raw.mapv(|v| if v.is_nan() { f64::NAN } else { 0.0 });
    }
    let scale = 1.0 / (2.0 * max_abs);
    let mut curv = raw;
    for (c, &ok) in curv.iter_mut().zip(valid.iter()) {
        if c.is_nan() {
            continue;
        }
        *c = if ok {
            (*c * scale).clamp(-0.5, 0.5)
        } else {
            0.0
        };
    }
    curv
}

/// The eight neighbours of a cell plus its own value, or `None` if the cell
/// is on the border or any neighbour is non-finite.
struct Neighbourhood {
    z0: f64,
    n: f64,
    s: f64,
    e: f64,
    w: f64,
    ne: f64,
    nw: f64,
    se: f64,
    sw: f64,
}

fn neighbourhood(z: &Array2<f64>, i: usize, j: usize) -> Option<Neighbourhood> {
    let (rows, cols) = z.dim();
    if i == 0 || j == 0 || i + 1 >= rows || j + 1 >= cols {
        return None;
    }
    let nb = Neighbourhood {
        z0: z[[i, j]],
        n: z[[i - 1, j]],
        s: z[[i + 1, j]],
        e: z[[i, j + 1]],
        w: z[[i, j - 1]],
        ne: z[[i - 1, j + 1]],
        nw: z[[i - 1, j - 1]],
        se: z[[i + 1, j + 1]],
        sw: z[[i + 1, j - 1]],
    };
    let all_finite = [nb.z0, nb.n, nb.s, nb.e, nb.w, nb.ne, nb.nw, nb.se, nb.sw]
        .iter()
        .all(|v| v.is_finite());
    all_finite.then_some(nb)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    const DEG: f64 = std::f64::consts::PI / 180.0;

    /// A plane tilted so elevation decreases toward the south, on a 100 m grid.
    /// Rows increase southward, so z must decrease with row index.
    fn south_facing(grade: f64, cellsize: f64) -> Array2<f64> {
        // z(i) = -grade * (row * cellsize): higher to the north (row 0).
        Array2::from_shape_fn((5, 5), |(i, _)| -grade * (i as f64) * cellsize)
    }

    #[test]
    fn flat_terrain_is_zero() {
        let z = Array2::from_elem((4, 4), 1500.0);
        let (slope, aspect) = slope_aspect(&z, 100.0);
        assert!(slope.iter().all(|&s| s == 0.0));
        assert!(aspect.iter().all(|&a| a == 0.0));
        assert!(curvature(&z).iter().all(|&c| c == 0.0));
    }

    #[test]
    fn south_facing_slope_has_south_aspect() {
        let cellsize = 100.0;
        let grade = 0.1; // 10% slope
        let z = south_facing(grade, cellsize);
        let (slope, aspect) = slope_aspect(&z, cellsize);
        let s = slope[[2, 2]];
        let a = aspect[[2, 2]];
        // slope magnitude = atan(0.1) ≈ 5.71°.
        assert!((s - grade.atan()).abs() < 1e-9, "slope {s}");
        // Downhill toward south = 180°.
        assert!((a - 180.0 * DEG).abs() < 1e-9, "aspect {}", a / DEG);
    }

    #[test]
    fn east_facing_slope_has_east_aspect() {
        let cellsize = 100.0;
        let grade = 0.2;
        // z decreases eastward (column index up): downhill toward east = 90°.
        let z = Array2::from_shape_fn((5, 5), |(_, j)| -grade * (j as f64) * cellsize);
        let (_, aspect) = slope_aspect(&z, cellsize);
        let a = aspect[[2, 2]];
        assert!((a - 90.0 * DEG).abs() < 1e-9, "aspect {}", a / DEG);
    }

    #[test]
    fn ridge_is_convex_valley_is_concave() {
        // Pyramid peak at the centre: the apex is convex (positive curvature).
        let peak = Array2::from_shape_fn((5, 5), |(i, j)| {
            -((i as f64 - 2.0).abs() + (j as f64 - 2.0).abs())
        });
        let cp = curvature(&peak);
        assert!(cp[[2, 2]] > 0.0, "peak {}", cp[[2, 2]]);

        // Inverted pyramid: the basin floor is concave (negative curvature).
        let basin = peak.mapv(|v| -v);
        let cb = curvature(&basin);
        assert!(cb[[2, 2]] < 0.0, "basin {}", cb[[2, 2]]);

        assert!(cp.iter().all(|&c| (-0.5..=0.5).contains(&c)));
        assert!(cb.iter().all(|&c| (-0.5..=0.5).contains(&c)));
    }

    #[test]
    fn nodata_propagates() {
        let mut z = Array2::from_elem((4, 4), 1000.0);
        z[[1, 1]] = f64::NAN;
        let (slope, aspect) = slope_aspect(&z, 100.0);
        assert!(slope[[1, 1]].is_nan());
        assert!(aspect[[1, 1]].is_nan());
        assert!(curvature(&z)[[1, 1]].is_nan());
        // A finite cell next to the NaN falls back to flat, not NaN.
        assert_eq!(slope[[2, 2]], 0.0);
    }
}
