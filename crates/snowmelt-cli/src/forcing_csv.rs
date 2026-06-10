//! Lectura de la serie de forzantes: CSV `date,temp_c,precip_mm`.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result, bail};

/// Un paso (diario) de forzantes meteorológicos.
#[derive(Debug, Clone)]
pub struct Record {
    pub date: String,
    pub temp_c: f64,
    pub precip_mm: f64,
}

/// Lee el CSV de forzantes. La primera línea puede ser un header
/// (se detecta porque sus columnas numéricas no parsean).
pub fn read(path: &Path) -> Result<Vec<Record>> {
    let text =
        fs::read_to_string(path).with_context(|| format!("no se pudo leer {}", path.display()))?;
    let mut records = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split(',').map(str::trim).collect();
        if parts.len() != 3 {
            bail!(
                "{}:{}: se esperaban 3 columnas (date,temp_c,precip_mm), hay {}",
                path.display(),
                i + 1,
                parts.len()
            );
        }
        match (parts[1].parse::<f64>(), parts[2].parse::<f64>()) {
            (Ok(temp_c), Ok(precip_mm)) => {
                if precip_mm < 0.0 {
                    bail!(
                        "{}:{}: precipitación negativa ({precip_mm})",
                        path.display(),
                        i + 1
                    );
                }
                records.push(Record {
                    date: parts[0].to_string(),
                    temp_c,
                    precip_mm,
                });
            }
            _ if i == 0 => continue, // header
            _ => bail!(
                "{}:{}: no se pudieron parsear temp/precip: `{line}`",
                path.display(),
                i + 1
            ),
        }
    }
    if records.is_empty() {
        bail!("{}: la serie de forzantes está vacía", path.display());
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_tmp(name: &str, content: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn parses_with_and_without_header() {
        let p = write_tmp(
            "snowmelt_forcing_hdr.csv",
            "date,temp_c,precip_mm\n2025-06-01,-2.0,5.0\n2025-06-02,3.5,0.0\n",
        );
        let recs = read(&p).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].temp_c, -2.0);
        fs::remove_file(&p).ok();

        let p = write_tmp("snowmelt_forcing_nohdr.csv", "2025-06-01,-2.0,5.0\n");
        assert_eq!(read(&p).unwrap().len(), 1);
        fs::remove_file(&p).ok();
    }

    #[test]
    fn rejects_negative_precip_and_bad_rows() {
        let p = write_tmp("snowmelt_forcing_neg.csv", "2025-06-01,-2.0,-5.0\n");
        assert!(read(&p).is_err());
        fs::remove_file(&p).ok();

        let p = write_tmp("snowmelt_forcing_bad.csv", "2025-06-01,-2.0,5.0\nfoo,bar\n");
        assert!(read(&p).is_err());
        fs::remove_file(&p).ok();
    }
}
