# snowmelt-rs

Motor de balance nival distribuido sobre DEM, en Rust. Tres modos de
derretimiento — grado-día, ETI (Pellicciotti et al. 2005) y balance de
energía completo con cold content — con acumulación/ablación de SWE por
celda, partición lluvia-nieve por temperatura, albedo dinámico por edad,
gradiente orográfico de precipitación y radiación potencial derivada del
terreno (SurtGIS, con sombreado por horizonte), orientado a hidrología
andina.

Parte de la familia de motores Rust: SurtGIS, Hydroflux, Smelt, Anvil,
Cantus, Criterium.

## Estructura

| Crate | Descripción |
|---|---|
| `snowmelt-core` | Modelo sin I/O: estado SWE (`ndarray`), grado-día/ETI, albedo dinámico por edad, lapse rate, partición lluvia-nieve, gradiente orográfico. Paralelo por celda con Rayon. |
| `snowmelt-cli` | Binario `snowmelt`: lee DEM (ESRI ASCII Grid) + forzantes CSV, calcula radiación potencial (SurtGIS terrain, con sombreado por horizonte opcional) y escribe serie agregada + grilla final de SWE. |
| `snowmelt-python` | Bindings PyO3/numpy (`import snowmelt`); compilar con `maturin develop -m crates/snowmelt-python/Cargo.toml`. |

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
rellenan con radiación de terreno plano. Con `--horizon-shading` se
precalculan ángulos de horizonte y la radiación directa respeta las
sombras del terreno circundante (memoria ≈ 8·direcciones·celdas bytes).

El albedo puede ser constante (`--albedo`) o dinámico con `--albedo-tau τ`:
`α(t) = α_min + (α_fresh − α_min)·exp(−t/τ)`, reiniciando a fresco cuando
la nevada del paso supera `--albedo-refresh` mm.

### Balance de energía

Con `--energy-balance` el derretimiento se calcula desde el flujo neto de
energía (W m⁻²) en vez del índice de temperatura:

```
Q = (1 − α)·G + LW_in − LW_out + Q_H + Q_E + Q_G
```

con onda larga de Brutsaert (1975) y superficie a `min(T_a, 0 °C)`,
flujos turbulentos bulk (`--wind`, `--rh`, `--exchange-coeff`), presión
del aire desde la elevación de cada celda y calor de suelo
(`--ground-heat`). El pack acumula **cold content** en días de balance
negativo (cap `c_ice·SWE·t_cold_max`) y la energía positiva lo paga
antes de derretir (L_f = 334 kJ/kg). Simplificaciones v0.4: sin calor de
lluvia sobre nieve, sin pérdida de masa por sublimación, cielo despejado.

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
| `--horizon-shading` | off | Sombreado por horizonte topográfico |
| `--horizon-radius` | 100 | Radio de búsqueda del horizonte [celdas] |
| `--horizon-directions` | 36 | Direcciones acimutales del horizonte |
| `--albedo-tau` | — | τ del decaimiento de albedo [días]; activa el modo dinámico |
| `--albedo-fresh` | 0.85 | Albedo de nieve fresca (modo dinámico) |
| `--albedo-min` | 0.4 | Albedo asintótico de nieve vieja (modo dinámico) |
| `--albedo-refresh` | 1.0 | Nevada [mm] que reinicia el albedo a fresco |
| `--energy-balance` | off | Derretimiento por balance de energía (ignora ddf/srf) |
| `--wind` | 2.0 | Viento [m/s] (modo EB) |
| `--rh` | 0.6 | Humedad relativa 0–1 (modo EB) |
| `--snow-emissivity` | 0.98 | Emisividad de la nieve (modo EB) |
| `--exchange-coeff` | 0.0015 | Coef. de intercambio turbulento (modo EB) |
| `--ground-heat` | 1.0 | Calor de suelo [W/m²] (modo EB) |
| `--t-cold-max` | 10.0 | Enfriamiento máximo del pack [K] (cold content) |

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

## Python

```python
import numpy as np, snowmelt

dem = np.loadtxt("dem.txt")           # 2D float64, NaN = nodata
p = snowmelt.Params(srf=0.2, albedo_tau=6.0)
m = snowmelt.SnowModel(dem, p)
out = m.step_uniform(t_ref=-3.0, z_ref=2500.0, precip=10.0,
                     radiation=np.full(dem.shape, 150.0))
out["melt"], m.swe(), m.albedo()
```

Compilación: `pip install maturin && maturin develop --release -m crates/snowmelt-python/Cargo.toml`.

## Roadmap (v0.5)

- Validación contra MODIS de cobertura nival y caudales DGA en cuenca andina real.
- Interfaz de aporte de deshielo hacia rainflow/Hydroflux.
- Calor de lluvia sobre nieve y sublimación con pérdida de masa.

## Licencia

MIT OR Apache-2.0
