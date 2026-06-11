# snowmelt-rs

Motor de balance nival distribuido sobre DEM, en Rust. Modelos grado-día
y ETI (enhanced temperature-index, Pellicciotti et al. 2005) con
acumulación/ablación de SWE por celda, partición lluvia-nieve por
temperatura, gradiente orográfico de precipitación y radiación potencial
derivada del terreno (SurtGIS), orientado a hidrología andina.

Parte de la familia de motores Rust: SurtGIS, Hydroflux, Smelt, Anvil,
Cantus, Criterium.

## Estructura

| Crate | Descripción |
|---|---|
| `snowmelt-core` | Modelo sin I/O: estado SWE (`ndarray`), grado-día/ETI, lapse rate, partición lluvia-nieve, gradiente orográfico. Paralelo por celda con Rayon. |
| `snowmelt-cli` | Binario `snowmelt`: lee DEM (ESRI ASCII Grid) + forzantes CSV, calcula radiación potencial (SurtGIS terrain) y escribe serie agregada + grilla final de SWE. |

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

Para el modo ETI con radiación (recomendado en cuencas andinas):

```bash
cargo run -q -p snowmelt-cli --release -- \
  --dem examples/dem.asc --forcing examples/forcing.csv \
  --out-dir out --z-ref 2500 \
  --srf 0.2 --albedo 0.6 --latitude -33.5
```

### Forzantes

CSV diario `date,temp_c,precip_mm` (header opcional). La temperatura se
extrapola al DEM con el gradiente vertical: `t(z) = t_ref + lapse·(z − z_ref)`.
La precipitación puede llevar gradiente orográfico:
`p(z) = p_ref·(1 + grad·(z − z_ref))`, acotada a ≥ 0 (`--precip-gradient`).

### Derretimiento ETI

Con `--srf > 0` el derretimiento usa el índice de temperatura mejorado
(Pellicciotti et al. 2005): para `T > t_melt`,

```
M = ddf·(T − t_melt) + srf·(1 − albedo)·G   [mm/día]
```

donde `G` es la radiación de cielo despejado (W m⁻² medios diarios)
calculada desde el DEM vía SurtGIS: slope/aspect (Horn) + geometría solar
por día del año, con transmitancia simple o turbidez de Linke
(`--linke-turbidity`). Los bordes del raster (sin vecindario 3×3) se
rellenan con radiación de terreno plano. Sin sombreado por horizonte en
v0.2 (`solar_radiation_shadowed` de SurtGIS queda para v0.3).

### Parámetros (defaults)

| Flag | Default | Significado |
|---|---|---|
| `--ddf` | 4.0 | Factor grado-día [mm °C⁻¹ día⁻¹] |
| `--t-melt` | 0.0 | Umbral de fusión [°C] |
| `--t-snow` | 0.0 | Bajo esto, 100 % nieve [°C] |
| `--t-rain` | 2.0 | Sobre esto, 100 % lluvia [°C] (lineal entremedio) |
| `--lapse-rate` | −0.0065 | Gradiente térmico [°C m⁻¹] |
| `--z-ref` | media del DEM | Elevación de la estación de forzantes [m] |
| `--srf` | 0.0 | Factor de radiación ETI [mm día⁻¹ (W m⁻²)⁻¹]; típico ~0.2 |
| `--albedo` | 0.6 | Albedo de la nieve (término radiativo) |
| `--precip-gradient` | 0.0 | Gradiente orográfico de precipitación [m⁻¹] |
| `--latitude` | — | Latitud [°]; requerida si `--srf > 0` |
| `--transmittance` | 0.7 | Transmitancia atmosférica de cielo despejado |
| `--linke-turbidity` | — | Turbidez de Linke (reemplaza a transmittance) |

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
(reanálisis o interpolación), paso no diario vía `step_days`, y radiación
para el término ETI vía `step_radiation(&forcing, Some(rad.view()), dt)`.

## Invariantes

- Balance de masa por celda y paso: `snowfall + rain == precip`, `Δswe == snowfall − melt`.
- Celdas nodata (`NaN` en el DEM) propagan `NaN` en estado y flujos; los agregados las excluyen.
- El derretimiento está acotado por el SWE disponible.

## Tests

```bash
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
```

## Roadmap (v0.3)

- Sombreado por horizonte (`solar_radiation_shadowed` de SurtGIS) y albedo
  con decaimiento por edad de la nieve.
- Balance de energía completo.
- Línea de nieves desde percepción remota; validación contra MODIS y caudales DGA.
- Bindings Python (PyO3) e interfaz de aporte hacia rainflow/Hydroflux.

## Licencia

MIT OR Apache-2.0
