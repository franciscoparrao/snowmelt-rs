//! Integration test: a synthetic accumulation–ablation season over a
//! small sloped DEM, checking physical plausibility and mass balance.

use ndarray::Array2;
use snowmelt_core::{DegreeDayParams, Dem, Forcing, SnowModel};

fn sloped_dem() -> Dem {
    // 4x4, 1000 m at the corner up to 4000 m.
    Dem::new(Array2::from_shape_fn((4, 4), |(i, j)| {
        1000.0 + 1000.0 * (i + j) as f64 / 2.0
    }))
    .unwrap()
}

#[test]
fn seasonal_cycle_accumulates_then_melts_out() {
    let mut model = SnowModel::new(sloped_dem(), DegreeDayParams::default()).unwrap();
    let z_ref = 1000.0;

    // 30 cold, snowy days.
    let winter = Forcing::Uniform {
        t_ref: -5.0,
        z_ref,
        precip: 10.0,
    };
    let mut peak_swe = 0.0_f64;
    for _ in 0..30 {
        let out = model.step(&winter).unwrap();
        peak_swe = model.summarize(&out).mean_swe;
    }
    assert!((peak_swe - 300.0).abs() < 1e-9, "peak {peak_swe}");

    // Spring: +15 °C at 1000 m. Low cells melt fast, high cells linger.
    let spring = Forcing::Uniform {
        t_ref: 15.0,
        z_ref,
        precip: 0.0,
    };
    let out = model.step(&spring).unwrap();
    let summary = model.summarize(&out);
    let swe = model.swe();
    assert!(
        swe[[3, 3]] > swe[[0, 0]],
        "high cell should retain more snow"
    );
    assert!(summary.mean_melt > 0.0);

    // Melt everything out with a warmer forcing: at +25 °C (1000 m) even the
    // 4000 m cell sits at 25 - 6.5*3 = +5.5 °C. Total melt must equal total
    // accumulation.
    let summer = Forcing::Uniform {
        t_ref: 25.0,
        z_ref,
        precip: 0.0,
    };
    let mut total_melt = summary.mean_melt;
    for _ in 0..60 {
        let out = model.step(&summer).unwrap();
        total_melt += model.summarize(&out).mean_melt;
    }
    let final_summary = {
        let out = model.step(&summer).unwrap();
        total_melt += model.summarize(&out).mean_melt;
        model.summarize(&out)
    };
    assert!(final_summary.mean_swe.abs() < 1e-6);
    assert_eq!(final_summary.snow_cover_fraction, 0.0);
    assert!((total_melt - 300.0).abs() < 1e-6, "melt {total_melt}");
}

#[test]
fn run_returns_one_summary_per_step() {
    let mut model = SnowModel::new(sloped_dem(), DegreeDayParams::default()).unwrap();
    let forcings: Vec<Forcing> = (0..5)
        .map(|d| Forcing::Uniform {
            t_ref: -2.0 + d as f64,
            z_ref: 1000.0,
            precip: 4.0,
        })
        .collect();
    let summaries = model.run(&forcings).unwrap();
    assert_eq!(summaries.len(), 5);
    assert!(summaries.iter().all(|s| s.mean_swe.is_finite()));
}
