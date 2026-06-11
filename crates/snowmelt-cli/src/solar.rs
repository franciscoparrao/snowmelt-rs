//! Radiación solar potencial sobre el DEM, vía SurtGIS (terrain).
//!
//! Calcula slope/aspect del DEM una vez y, por cada día del año requerido,
//! la radiación de cielo despejado (beam + difusa + reflejada) integrada
//! diaria. Devuelve W m⁻² medios diarios, que es la unidad que espera el
//! término radiativo ETI de `snowmelt-core`.

use anyhow::{Context, Result, anyhow};
use ndarray::{Array2, Zip};
use surtgis_algorithms::terrain::{
    AspectOutput, SlopeParams, SlopeUnits, SolarParams, aspect, slope, solar_radiation,
};
use surtgis_core::{GeoTransform, Raster};

use crate::asc::AscHeader;

/// Slope y aspect (radianes) precalculados del DEM.
pub struct Terrain {
    slope_rad: Raster<f64>,
    aspect_rad: Raster<f64>,
}

impl Terrain {
    /// Deriva slope/aspect del DEM (Horn). Los bordes y vecinos de nodata
    /// quedan `NaN`; ver [`Terrain::potential_radiation`] para el relleno.
    pub fn from_dem(elevation: &Array2<f64>, header: &AscHeader) -> Result<Self> {
        let mut dem = Raster::from_array(elevation.clone());
        // Origen = esquina superior izquierda; alto de píxel negativo.
        let origin_y = header.yll + header.nrows as f64 * header.cellsize;
        dem.set_transform(GeoTransform::new(
            header.xll,
            origin_y,
            header.cellsize,
            -header.cellsize,
        ));
        let slope_rad = slope(
            &dem,
            SlopeParams {
                units: SlopeUnits::Radians,
                z_factor: 1.0,
            },
        )
        .map_err(|e| anyhow!("slope: {e}"))?;
        let aspect_rad = aspect(&dem, AspectOutput::Radians).map_err(|e| anyhow!("aspect: {e}"))?;
        Ok(Self {
            slope_rad,
            aspect_rad,
        })
    }

    /// Radiación potencial de cielo despejado para `day` (1–365), en W m⁻²
    /// medios diarios.
    ///
    /// Las celdas donde slope/aspect son `NaN` pero el DEM es válido
    /// (bordes del raster, vecinos de nodata) se rellenan con la radiación
    /// de terreno plano del mismo día, para no encoger el dominio.
    pub fn potential_radiation(
        &self,
        elevation: &Array2<f64>,
        day: u32,
        latitude: f64,
        transmittance: f64,
        linke_turbidity: Option<f64>,
        albedo: f64,
    ) -> Result<Array2<f64>> {
        let params = SolarParams {
            day,
            latitude,
            transmittance,
            linke_turbidity,
            albedo,
            ..SolarParams::default()
        };
        let result = solar_radiation(&self.slope_rad, &self.aspect_rad, params.clone())
            .map_err(|e| anyhow!("solar_radiation: {e}"))?;
        let flat = flat_radiation(params).context("radiación de terreno plano")?;

        // Wh/m²/día → W/m² medio diario; relleno plano donde corresponde.
        let mut rad = result.total.data().mapv(|v| v / 24.0);
        Zip::from(&mut rad).and(elevation).for_each(|r, &z| {
            if z.is_finite() && !r.is_finite() {
                *r = flat;
            }
        });
        Ok(rad)
    }
}

/// Radiación diaria (W m⁻²) de una superficie horizontal, mismo modelo.
fn flat_radiation(params: SolarParams) -> Result<f64> {
    let zeros = Raster::from_array(Array2::zeros((3, 3)));
    let result = solar_radiation(&zeros, &zeros, params).map_err(|e| anyhow!("{e}"))?;
    let center = result.total.data()[[1, 1]];
    if !center.is_finite() {
        anyhow::bail!("la radiación de terreno plano resultó no finita");
    }
    Ok(center / 24.0)
}

/// Día del año (1–365) desde una fecha `YYYY-MM-DD`. El 29 de febrero se
/// trata como día 59 (28-feb) y los días posteriores se ajustan, de modo
/// que el resultado siempre cae en 1–365 (rango que exige SurtGIS).
pub fn day_of_year(date: &str) -> Result<u32> {
    let parts: Vec<&str> = date.split('-').collect();
    if parts.len() != 3 {
        anyhow::bail!("fecha inválida (se espera YYYY-MM-DD): `{date}`");
    }
    let month: u32 = parts[1]
        .parse()
        .with_context(|| format!("mes inválido en `{date}`"))?;
    let day: u32 = parts[2]
        .parse()
        .with_context(|| format!("día inválido en `{date}`"))?;
    if !(1..=12).contains(&month) {
        anyhow::bail!("mes fuera de rango en `{date}`");
    }
    const CUM: [u32; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let max_day = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][month as usize - 1];
    if day == 0 || day > max_day {
        anyhow::bail!("día fuera de rango en `{date}`");
    }
    // 29-feb colapsa al 28-feb (día 59).
    let doy = CUM[month as usize - 1] + day.min(if month == 2 { 28 } else { day });
    Ok(doy.min(365))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn day_of_year_handles_boundaries_and_leap() {
        assert_eq!(day_of_year("2025-01-01").unwrap(), 1);
        assert_eq!(day_of_year("2025-12-31").unwrap(), 365);
        assert_eq!(day_of_year("2025-06-21").unwrap(), 172);
        assert_eq!(day_of_year("2024-02-29").unwrap(), 59);
        assert!(day_of_year("2025-13-01").is_err());
        assert!(day_of_year("2025-02-30").is_err());
        assert!(day_of_year("no-fecha").is_err());
    }

    #[test]
    fn radiation_is_positive_and_fills_borders() {
        // DEM 5x5 inclinado, una celda nodata interior.
        let mut elev = Array2::from_shape_fn((5, 5), |(i, j)| 1000.0 + 50.0 * (i + j) as f64);
        elev[[2, 2]] = f64::NAN;
        let header = AscHeader {
            ncols: 5,
            nrows: 5,
            xll: 0.0,
            yll: 0.0,
            cellsize: 100.0,
            nodata: -9999.0,
        };
        let terrain = Terrain::from_dem(&elev, &header).unwrap();
        // Solsticio de invierno austral, latitud andina central.
        let rad = terrain
            .potential_radiation(&elev, 172, -33.5, 0.7, None, 0.6)
            .unwrap();
        // Celdas válidas (incluidos los bordes) son finitas y positivas.
        for ((i, j), &v) in rad.indexed_iter() {
            if elev[[i, j]].is_finite() {
                assert!(v.is_finite() && v > 0.0, "celda ({i},{j}): {v}");
            }
        }
        // En invierno austral a -33.5°, la radiación media diaria es modesta.
        let mean: f64 = rad.iter().filter(|v| v.is_finite()).sum::<f64>()
            / rad.iter().filter(|v| v.is_finite()).count() as f64;
        assert!(mean > 10.0 && mean < 500.0, "media {mean}");
    }
}
