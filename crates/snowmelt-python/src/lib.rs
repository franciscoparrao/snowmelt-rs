//! Bindings Python del modelo snowmelt-rs.
//!
//! Expone `SnowModel` con grillas numpy (float64, 2D): SWE, albedo y edad
//! de nieve son arrays con `NaN` en nodata. Compilar con maturin:
//! `maturin develop -m crates/snowmelt-python/Cargo.toml`.

use numpy::{PyArray2, PyReadonlyArray2, ToPyArray};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use snowmelt_core as core;

fn to_py_err(e: core::SnowmeltError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Parámetros del modelo (grado-día / ETI).
#[pyclass(name = "Params", from_py_object)]
#[derive(Clone)]
struct PyParams {
    inner: core::DegreeDayParams,
}

#[pymethods]
impl PyParams {
    /// Crea parámetros; los no especificados usan los defaults del core
    /// (ddf=4.0, t_melt=0, t_snow=0, t_rain=2, lapse_rate=-0.0065,
    /// srf=0, albedo=0.6, precip_gradient=0).
    #[new]
    #[pyo3(signature = (ddf=4.0, t_melt=0.0, t_snow=0.0, t_rain=2.0,
        lapse_rate=-0.0065, srf=0.0, albedo=0.6, precip_gradient=0.0,
        albedo_tau=None, albedo_fresh=0.85, albedo_min=0.4, albedo_refresh=1.0,
        energy_balance=false, wind=2.0, rh=0.6, snow_emissivity=0.98,
        exchange_coeff=0.0015, ground_heat=1.0, t_cold_max=10.0))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        ddf: f64,
        t_melt: f64,
        t_snow: f64,
        t_rain: f64,
        lapse_rate: f64,
        srf: f64,
        albedo: f64,
        precip_gradient: f64,
        albedo_tau: Option<f64>,
        albedo_fresh: f64,
        albedo_min: f64,
        albedo_refresh: f64,
        energy_balance: bool,
        wind: f64,
        rh: f64,
        snow_emissivity: f64,
        exchange_coeff: f64,
        ground_heat: f64,
        t_cold_max: f64,
    ) -> PyResult<Self> {
        let inner = core::DegreeDayParams {
            ddf,
            t_melt,
            t_snow,
            t_rain,
            lapse_rate,
            srf,
            albedo,
            albedo_decay: albedo_tau.map(|tau| core::AlbedoDecay {
                albedo_fresh,
                albedo_min,
                tau_days: tau,
                refresh_swe_mm: albedo_refresh,
            }),
            energy_balance: energy_balance.then_some(core::EnergyBalanceParams {
                wind_speed: wind,
                rel_humidity: rh,
                snow_emissivity,
                exchange_coeff,
                ground_heat,
                t_cold_max,
            }),
            precip_gradient,
        };
        inner.validate().map_err(to_py_err)?;
        Ok(Self { inner })
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

/// Modelo de nieve distribuido sobre un DEM (NaN = nodata).
#[pyclass(name = "SnowModel")]
struct PySnowModel {
    inner: core::SnowModel,
}

#[pymethods]
impl PySnowModel {
    /// Crea el modelo desde un DEM 2D float64 (NaN = nodata) y parámetros.
    /// `initial_swe` opcional (mm, misma forma que el DEM).
    #[new]
    #[pyo3(signature = (dem, params, initial_swe=None))]
    fn new(
        dem: PyReadonlyArray2<'_, f64>,
        params: PyParams,
        initial_swe: Option<PyReadonlyArray2<'_, f64>>,
    ) -> PyResult<Self> {
        let dem = core::Dem::new(dem.as_array().to_owned()).map_err(to_py_err)?;
        let inner = match initial_swe {
            Some(swe) => {
                core::SnowModel::with_initial_swe(dem, params.inner, swe.as_array().to_owned())
            }
            None => core::SnowModel::new(dem, params.inner),
        }
        .map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Avanza un paso con forzante uniforme (estación + lapse rate).
    /// `radiation`: grilla opcional de W/m² medios diarios (requerida si srf > 0).
    /// Devuelve dict con snowfall, rain, melt (mm, grillas numpy).
    #[pyo3(signature = (t_ref, z_ref, precip, radiation=None, dt_days=1.0))]
    fn step_uniform<'py>(
        &mut self,
        py: Python<'py>,
        t_ref: f64,
        z_ref: f64,
        precip: f64,
        radiation: Option<PyReadonlyArray2<'py, f64>>,
        dt_days: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let forcing = core::Forcing::Uniform {
            t_ref,
            z_ref,
            precip,
        };
        let rad = radiation.as_ref().map(|r| r.as_array());
        let out = self
            .inner
            .step_radiation(&forcing, rad, dt_days)
            .map_err(to_py_err)?;
        step_output_dict(py, &out)
    }

    /// Avanza un paso con grillas distribuidas de temperatura (°C) y
    /// precipitación (mm), y radiación opcional (W/m²).
    #[pyo3(signature = (temp, precip, radiation=None, dt_days=1.0))]
    fn step_distributed<'py>(
        &mut self,
        py: Python<'py>,
        temp: PyReadonlyArray2<'py, f64>,
        precip: PyReadonlyArray2<'py, f64>,
        radiation: Option<PyReadonlyArray2<'py, f64>>,
        dt_days: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let forcing = core::Forcing::Distributed {
            temp: temp.as_array().to_owned(),
            precip: precip.as_array().to_owned(),
        };
        let rad = radiation.as_ref().map(|r| r.as_array());
        let out = self
            .inner
            .step_radiation(&forcing, rad, dt_days)
            .map_err(to_py_err)?;
        step_output_dict(py, &out)
    }

    /// SWE actual (mm) como array numpy.
    fn swe<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.inner.swe().to_pyarray(py)
    }

    /// Albedo por celda usado en el último paso.
    fn albedo<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.inner.albedo().to_pyarray(py)
    }

    /// Edad de la nieve en días (modo albedo dinámico).
    fn snow_age<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.inner.snow_age().to_pyarray(py)
    }

    /// Cold content del pack [J/m²] (modo balance de energía).
    fn cold_content<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.inner.cold_content().to_pyarray(py)
    }
}

fn step_output_dict<'py>(py: Python<'py>, out: &core::StepOutput) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("snowfall", out.snowfall.to_pyarray(py))?;
    d.set_item("rain", out.rain.to_pyarray(py))?;
    d.set_item("melt", out.melt.to_pyarray(py))?;
    d.set_item("runoff", out.runoff().to_pyarray(py))?;
    Ok(d)
}

/// Módulo Python `snowmelt`.
#[pymodule]
fn snowmelt(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyParams>()?;
    m.add_class::<PySnowModel>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
