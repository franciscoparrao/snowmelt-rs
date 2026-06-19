//! Continuous skill metrics for SWE (and any gridded/series comparison).
//!
//! The cover-validation tooling scores a binary snow/no-snow map; comparing
//! modelled **snow water equivalent** against a reference (e.g. the Andes
//! Snow Reanalysis) needs continuous metrics instead. This module computes
//! them over the pairs where both values are finite, so nodata (`NaN`) is
//! skipped automatically. It is generic over any iterator of
//! `(model, reference)` pairs, so it serves both a per-date grid comparison
//! and a basin-mean time series.

/// Continuous goodness-of-fit between a modelled field/series and a reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContinuousSkill {
    /// Number of co-valid pairs.
    pub n: usize,
    /// Root mean squared error (model − reference).
    pub rmse: f64,
    /// Mean bias error (model − reference).
    pub mbe: f64,
    /// Mean absolute error.
    pub mae: f64,
    /// Pearson correlation coefficient (`NaN` if either side has zero variance).
    pub correlation: f64,
    /// Kling–Gupta efficiency `1 − √((r−1)² + (α−1)² + (β−1)²)`, with
    /// `α = σ_model/σ_ref` and `β = μ_model/μ_ref` (`NaN` if the reference
    /// has zero mean or variance).
    pub kge: f64,
}

/// Computes [`ContinuousSkill`] over the `(model, reference)` pairs where
/// both are finite. Returns `None` if fewer than two such pairs exist.
pub fn continuous_skill<I>(pairs: I) -> Option<ContinuousSkill>
where
    I: IntoIterator<Item = (f64, f64)>,
{
    let mut n = 0usize;
    let mut sum_m = 0.0;
    let mut sum_r = 0.0;
    let mut sum_abs = 0.0;
    let mut sum_sq = 0.0;
    // Co-valid pairs are buffered for the second (variance/covariance) pass.
    let mut buf: Vec<(f64, f64)> = Vec::new();
    for (m, r) in pairs {
        if m.is_finite() && r.is_finite() {
            n += 1;
            let d = m - r;
            sum_m += m;
            sum_r += r;
            sum_abs += d.abs();
            sum_sq += d * d;
            buf.push((m, r));
        }
    }
    if n < 2 {
        return None;
    }
    let nf = n as f64;
    let mean_m = sum_m / nf;
    let mean_r = sum_r / nf;
    let rmse = (sum_sq / nf).sqrt();
    let mbe = mean_m - mean_r;
    let mae = sum_abs / nf;

    let mut cov = 0.0;
    let mut var_m = 0.0;
    let mut var_r = 0.0;
    for (m, r) in &buf {
        let dm = m - mean_m;
        let dr = r - mean_r;
        cov += dm * dr;
        var_m += dm * dm;
        var_r += dr * dr;
    }
    let std_m = (var_m / nf).sqrt();
    let std_r = (var_r / nf).sqrt();
    let correlation = if std_m > 0.0 && std_r > 0.0 {
        cov / (nf * std_m * std_r)
    } else {
        f64::NAN
    };
    let kge = if std_r > 0.0 && mean_r != 0.0 && correlation.is_finite() {
        let alpha = std_m / std_r;
        let beta = mean_m / mean_r;
        1.0 - ((correlation - 1.0).powi(2) + (alpha - 1.0).powi(2) + (beta - 1.0).powi(2)).sqrt()
    } else {
        f64::NAN
    };

    Some(ContinuousSkill {
        n,
        rmse,
        mbe,
        mae,
        correlation,
        kge,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn identical_series_is_perfect() {
        let v = [1.0, 2.0, 3.0, 4.0, 5.0];
        let s = continuous_skill(v.iter().map(|&x| (x, x))).unwrap();
        assert_eq!(s.n, 5);
        assert!(approx(s.rmse, 0.0));
        assert!(approx(s.mbe, 0.0));
        assert!(approx(s.mae, 0.0));
        assert!(approx(s.correlation, 1.0));
        assert!(approx(s.kge, 1.0));
    }

    #[test]
    fn constant_offset_gives_bias_and_unit_correlation() {
        let r = [10.0, 20.0, 30.0, 40.0];
        // Model = reference + 5: perfectly correlated, biased high.
        let s = continuous_skill(r.iter().map(|&x| (x + 5.0, x))).unwrap();
        assert!(approx(s.mbe, 5.0));
        assert!(approx(s.rmse, 5.0));
        assert!(approx(s.correlation, 1.0));
        // α = 1 (same spread), β = mean_m/mean_r = 30/25 = 1.2, r = 1.
        // KGE = 1 − |β−1| = 1 − 0.2 = 0.8.
        assert!(approx(s.kge, 0.8), "kge {}", s.kge);
    }

    #[test]
    fn skips_nonfinite_pairs() {
        let pairs = [
            (1.0, 1.0),
            (f64::NAN, 2.0),
            (3.0, f64::NAN),
            (4.0, 4.0),
            (6.0, 6.0),
        ];
        let s = continuous_skill(pairs).unwrap();
        assert_eq!(s.n, 3); // only (1,1), (4,4), (6,6)
        assert!(approx(s.rmse, 0.0));
    }

    #[test]
    fn too_few_pairs_returns_none() {
        assert!(continuous_skill([(1.0, 1.0)]).is_none());
        assert!(continuous_skill([(f64::NAN, 1.0), (2.0, f64::NAN)]).is_none());
    }

    #[test]
    fn zero_reference_variance_has_nan_correlation_and_kge() {
        // Reference constant → correlation and KGE undefined, but RMSE/MBE fine.
        let s = continuous_skill([(1.0, 5.0), (3.0, 5.0), (5.0, 5.0)]).unwrap();
        assert!(approx(s.mbe, -2.0)); // mean model 3 − 5
        assert!(s.correlation.is_nan());
        assert!(s.kge.is_nan());
    }
}
