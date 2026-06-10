//! Lectura/escritura de grillas en formato ESRI ASCII Grid (.asc).

use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use anyhow::{Context, Result, bail};
use ndarray::{Array2, ArrayView2};

/// Cabecera de una grilla ASCII (georreferenciación mínima).
#[derive(Debug, Clone)]
pub struct AscHeader {
    pub ncols: usize,
    pub nrows: usize,
    pub xll: f64,
    pub yll: f64,
    pub cellsize: f64,
    pub nodata: f64,
}

/// Grilla leída: cabecera + datos con `NaN` en celdas nodata.
/// La fila 0 corresponde al borde norte (convención .asc).
#[derive(Debug)]
pub struct AscGrid {
    pub header: AscHeader,
    pub data: Array2<f64>,
}

/// Lee un archivo .asc. Acepta claves de cabecera en cualquier
/// capitalización y `xllcorner`/`xllcenter` indistintamente.
pub fn read(path: &Path) -> Result<AscGrid> {
    let text =
        fs::read_to_string(path).with_context(|| format!("no se pudo leer {}", path.display()))?;
    let mut tokens = text.split_whitespace().peekable();

    let mut ncols: Option<usize> = None;
    let mut nrows: Option<usize> = None;
    let mut xll: Option<f64> = None;
    let mut yll: Option<f64> = None;
    let mut cellsize: Option<f64> = None;
    let mut nodata = -9999.0;

    while let Some(&tok) = tokens.peek() {
        if !tok.chars().next().is_some_and(|c| c.is_ascii_alphabetic()) {
            break;
        }
        let key = tok.to_ascii_lowercase();
        tokens.next();
        let value = tokens
            .next()
            .with_context(|| format!("cabecera incompleta: falta el valor de `{key}`"))?;
        match key.as_str() {
            "ncols" => ncols = Some(value.parse().context("ncols inválido")?),
            "nrows" => nrows = Some(value.parse().context("nrows inválido")?),
            "xllcorner" | "xllcenter" => xll = Some(value.parse().context("xll inválido")?),
            "yllcorner" | "yllcenter" => yll = Some(value.parse().context("yll inválido")?),
            "cellsize" => cellsize = Some(value.parse().context("cellsize inválido")?),
            "nodata_value" => nodata = value.parse().context("nodata_value inválido")?,
            other => bail!("clave de cabecera desconocida: `{other}`"),
        }
    }

    let header = AscHeader {
        ncols: ncols.context("falta `ncols` en la cabecera")?,
        nrows: nrows.context("falta `nrows` en la cabecera")?,
        xll: xll.context("falta `xllcorner` en la cabecera")?,
        yll: yll.context("falta `yllcorner` en la cabecera")?,
        cellsize: cellsize.context("falta `cellsize` en la cabecera")?,
        nodata,
    };

    let mut values = Vec::with_capacity(header.nrows * header.ncols);
    for tok in tokens {
        let v: f64 = tok
            .parse()
            .with_context(|| format!("valor de celda inválido: `{tok}`"))?;
        values.push(if v == nodata { f64::NAN } else { v });
    }
    let expected = header.nrows * header.ncols;
    if values.len() != expected {
        bail!(
            "se esperaban {expected} celdas ({} x {}), se encontraron {}",
            header.nrows,
            header.ncols,
            values.len()
        );
    }
    let data = Array2::from_shape_vec((header.nrows, header.ncols), values)
        .context("no se pudo construir la grilla")?;
    Ok(AscGrid { header, data })
}

/// Escribe una grilla .asc; las celdas `NaN` se emiten como `nodata`.
pub fn write(path: &Path, header: &AscHeader, data: ArrayView2<'_, f64>) -> Result<()> {
    let (nrows, ncols) = data.dim();
    if nrows != header.nrows || ncols != header.ncols {
        bail!(
            "la grilla ({nrows} x {ncols}) no coincide con la cabecera ({} x {})",
            header.nrows,
            header.ncols
        );
    }
    let mut out = String::new();
    let _ = writeln!(out, "ncols {}", header.ncols);
    let _ = writeln!(out, "nrows {}", header.nrows);
    let _ = writeln!(out, "xllcorner {}", header.xll);
    let _ = writeln!(out, "yllcorner {}", header.yll);
    let _ = writeln!(out, "cellsize {}", header.cellsize);
    let _ = writeln!(out, "NODATA_value {}", header.nodata);
    for row in data.rows() {
        let mut first = true;
        for &v in row {
            if !first {
                out.push(' ');
            }
            first = false;
            if v.is_finite() {
                let _ = write!(out, "{v:.3}");
            } else {
                let _ = write!(out, "{}", header.nodata);
            }
        }
        out.push('\n');
    }
    fs::write(path, out).with_context(|| format!("no se pudo escribir {}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_preserves_data_and_nodata() {
        let dir = std::env::temp_dir();
        let path = dir.join("snowmelt_asc_roundtrip_test.asc");
        let header = AscHeader {
            ncols: 3,
            nrows: 2,
            xll: 100.0,
            yll: 200.0,
            cellsize: 50.0,
            nodata: -9999.0,
        };
        let data = ndarray::array![[1.0, 2.5, f64::NAN], [4.0, 5.0, 6.0]];
        write(&path, &header, data.view()).unwrap();
        let grid = read(&path).unwrap();
        assert_eq!(grid.header.ncols, 3);
        assert_eq!(grid.header.nrows, 2);
        assert!(grid.data[[0, 2]].is_nan());
        assert_eq!(grid.data[[1, 0]], 4.0);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn rejects_cell_count_mismatch() {
        let dir = std::env::temp_dir();
        let path = dir.join("snowmelt_asc_badcount_test.asc");
        fs::write(
            &path,
            "ncols 2\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 1\n1 2 3\n",
        )
        .unwrap();
        assert!(read(&path).is_err());
        std::fs::remove_file(&path).ok();
    }
}
