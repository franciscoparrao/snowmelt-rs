# snowmelt-rs

Motor de balance nival distribuido sobre DEM, en Rust. Modelo grado-día
(temperature-index) con acumulación/ablación de SWE por celda y partición
lluvia-nieve por temperatura, orientado a hidrología andina.

Parte de la familia de motores Rust: SurtGIS, Hydroflux, Smelt, Anvil,
Cantus, Criterium.

## Estructura

| Crate | Descripción |
|---|---|
| `snowmelt-core` | Modelo sin I/O: estado SWE (`ndarray`), grado-día, lapse rate, partición lluvia-nieve. Paralelo por celda con Rayon. |
| `snowmelt-cli` | Binario `snowmelt`: lee DEM (ESRI ASCII Grid) + forzantes CSV, escribe serie agregada y grilla final de SWE. |

## Uso rápido

```bash
cargo run -q -p snowmelt-cli --release -- \
  --dem examples/dem.asc \
  --forcing examples/forcing.csv \
  --out-dir out \
  --z-ref 2500
```

Salidas en `out/`:

- `series.csv` — por paso: `date,snowfall_mm,rain_mm,melt_mm,runoff_mm,swe_mm,snow_cover_fraction` (medias sobre celdas válidas).
- `swe_final.asc` — grilla final de SWE (mm), con NODATA preservado.

### Forzantes

CSV diario `date,temp_c,precip_mm` (header opcional). La temperatura se
extrapola al DEM con el gradiente vertical: `t(z) = t_ref + lapse·(z − z_ref)`.
La precipitación se aplica uniforme (gradiente orográfico: pendiente para v0.2).

### Parámetros (defaults)

| Flag | Default | Significado |
|---|---|---|
| `--ddf` | 4.0 | Factor grado-día [mm °C⁻¹ día⁻¹] |
| `--t-melt` | 0.0 | Umbral de fusión [°C] |
| `--t-snow` | 0.0 | Bajo esto, 100 % nieve [°C] |
| `--t-rain` | 2.0 | Sobre esto, 100 % lluvia [°C] (lineal entremedio) |
| `--lapse-rate` | −0.0065 | Gradiente térmico [°C m⁻¹] |
| `--z-ref` | media del DEM | Elevación de la estación de forzantes [m] |

## API (snowmelt-core)

```rust
use ndarray::Array2;
use snowmelt_core::{Dem, DegreeDayParams, Forcing, SnowModel};

let dem = Dem::new(elevation_grid)?;            // NaN = nodata
let mut model = SnowModel::new(dem, DegreeDayParams::default())?;
let out = model.step(&Forcing::Uniform { t_ref: -3.0, z_ref: 1000.0, precip: 10.0 })?;
// out.snowfall / out.rain / out.melt / out.runoff(); model.swe()
```

También acepta `Forcing::Distributed { temp, precip }` con grillas completas
(reanálisis o interpolación) y paso no diario vía `step_days`.

## Invariantes

- Balance de masa por celda y paso: `snowfall + rain == precip`, `Δswe == snowfall − melt`.
- Celdas nodata (`NaN` en el DEM) propagan `NaN` en estado y flujos; los agregados las excluyen.
- El derretimiento está acotado por el SWE disponible.

## Tests

```bash
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
```

## Roadmap (v0.2)

- Balance de energía; radiación solar desde SurtGIS (terrain).
- Gradiente orográfico de precipitación.
- Línea de nieves desde percepción remota; validación contra MODIS y caudales DGA.
- Bindings Python (PyO3) e interfaz de aporte hacia rainflow/Hydroflux.

## Licencia

MIT OR Apache-2.0
