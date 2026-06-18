# snowmelt-rs

Motor de balance nival distribuido sobre DEM, en Rust. Tres modos de
derretimiento — grado-día, ETI (Pellicciotti et al. 2005) y balance de
energía completo (cold content, nubosidad, rain-on-snow, sublimación) —
con acumulación/ablación de SWE por celda, partición lluvia-nieve por
temperatura, albedo dinámico por edad, gradiente orográfico de
precipitación, radiación potencial derivada del terreno (SurtGIS, con
sombreado por horizonte) y ruteo por reservorio lineal hacia un
hidrograma. Orientado a hidrología andina. Validado contra MODIS
(cobertura nival, F1 0.83) y caudales CAMELS-CL (deshielo, NSE 0.66).

Parte de la familia de motores Rust: SurtGIS, Hydroflux, Smelt, Anvil,
Cantus, Criterium.

## Estructura

| Crate | Descripción |
|---|---|
| `snowmelt-core` | Modelo sin I/O: estado SWE (`ndarray`), grado-día/ETI, albedo dinámico por edad, lapse rate, partición lluvia-nieve, gradiente orográfico. Paralelo por celda con Rayon. |
| `snowmelt-cli` | Binario `snowmelt`: lee DEM (ESRI ASCII Grid) + forzantes CSV, calcula radiación potencial (SurtGIS terrain, con sombreado por horizonte opcional) y escribe serie agregada + grilla final de SWE. |
| `snowmelt-python` | Bindings PyO3/numpy (`import snowmelt`); compilar con `maturin develop -m crates/snowmelt-python/Cargo.toml`. |

El binario `snowmelt-validate` calcula métricas de cobertura (confusión,
F1, bias) entre grillas `.asc`.

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

donde `Q_R` (calor de lluvia sobre nieve) entra cuando hay precipitación
líquida. Onda larga de Brutsaert (1975) y superficie a `min(T_a, 0 °C)`,
flujos turbulentos bulk (`--wind`, `--rh`, `--exchange-coeff`), presión
del aire desde la elevación de cada celda y calor de suelo
(`--ground-heat`). La **nubosidad efectiva** (`--cloud-fraction N`)
atenúa la onda corta `(1 − 0.75·N³)` y refuerza la onda larga entrante
`(1 + 0.22·N²)`. El pack acumula **cold content** en días de balance
negativo (cap `c_ice·SWE·t_cold_max`) y la energía positiva lo paga
antes de derretir (L_f = 334 kJ/kg). El flujo latente negativo retira
masa por **sublimación** (L_s = 2.834 MJ/kg), reportada aparte; balance
de masa por celda: `Δswe = nieve − derretimiento − sublimación`.

### Downscaling topográfico (forzante sub-km)

Con `--downscale` el forzante distribuido se genera **desde el DEM** (estilo
MicroMet, Liston & Elder 2006) en lugar de leer grillas externas: a partir
del valor escalar del CSV se construyen campos de temperatura y
precipitación a la resolución del DEM con tres controles de terreno además
del lapse rate:

- **Temperatura**: `T(z) = T_ref + lapse·(z − z_ref) + κ·Ω_c`, donde `Ω_c`
  es la curvatura normalizada (`[−0.5, 0.5]`) y `κ` (`--temp-curvature`)
  enfría los valles cóncavos y templa las cumbres convexas — proxy diario
  del cold-air pooling.
- **Viento**: factor de terreno `W = 1 + γ_s·Ω_s + γ_c·Ω_c`, con `Ω_s` la
  pendiente en la dirección del viento (`--wind-dir`, grados desde donde
  sopla) y `Ω_c` la curvatura.
- **Precipitación**: factor de elevación de Thornton (1997)
  `(1 + f·Δz)/(1 − f·Δz)` (`--precip-elev-factor`) por un realce orográfico
  a barlovento `1 + γ_w·Ω_s` (`--precip-windward`) que aumenta las laderas
  que enfrentan el viento y seca el sotavento.

```bash
snowmelt --dem data/dem.asc --forcing data/forcing.csv --out-dir out \
  --z-ref 3117 --energy-balance --latitude -33.675 --lapse-rate -0.0075 \
  --downscale --temp-curvature 5 --wind-dir 300
```

Los derivados de terreno (slope/aspect por Horn, curvatura Liston-Elder) se
calculan en `snowmelt-core` sin dependencia GIS. Es excluyente con
`--temp-grids`/`--precip-grids` (ambos definen el forzante distribuido).
Sobre el Maipo alto la mejora es marginal y solo viene del término de
curvatura (ver [estudio de sensibilidad](validation/maipo-alto/FORCING_SENSITIVITY.md)).

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
| `--cloud-fraction` | 0.0 | Fracción efectiva de nubes 0–1 (modo EB) |
| `--downscale` | off | Downscaling topográfico del forzante desde el DEM |
| `--temp-curvature` | 0.0 | Coef. temperatura por curvatura [°C] (cold-air pooling) |
| `--precip-elev-factor` | 0.0 | Factor precipitación-elevación [km⁻¹] (Thornton) |
| `--precip-windward` | 0.0 | Realce orográfico a barlovento γ_w |
| `--wind-dir` | 300.0 | Dirección del viento dominante [° desde donde sopla] |
| `--wind-slope-weight` | 0.5 | Peso pendiente-en-viento γ_s |
| `--wind-curvature-weight` | 0.5 | Peso curvatura del factor de viento γ_c |

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

## Validación

Validación contra MODIS MOD10A1 en la cuenca alta del Maipo (temporada
2019, balance de energía + albedo dinámico, lapse calibrado −7.5 °C/km):
**accuracy 85.9%, F1 0.834, bias 1.11** (julio F1 0.91). Pipeline
reproducible, grid search de calibración y un
[estudio de sensibilidad de forzantes](validation/maipo-alto/FORCING_SENSITIVITY.md)
(uniforme vs precipitación/temperatura distribuida) en
[`validation/maipo-alto/`](validation/maipo-alto/README.md). El binario
`snowmelt-validate` calcula las métricas:

```bash
snowmelt-validate out/cover_2019-07-15.asc:data/modis_2019-07-15.asc ...
```

### Caudal (CAMELS-CL)

Validación de la componente hidrológica contra el Río Choapa en Cuncumén
(CAMELS-CL 4703002, cuenca nival, 38 años) con la cuenca discretizada en
bandas de elevación de igual área y ruteo por reservorio lineal: el
**deshielo reproduce la firma estacional del caudal** (corr. diaria 0.81,
NSE de forma 0.66, corr. del ciclo anual 0.88). Detalle e interfaz hacia
rainflow en [`validation/choapa-cuncumen/`](validation/choapa-cuncumen/README.md).

El forzante puede ser **distribuido por grillas** (`--precip-grids DIR`,
`--temp-grids DIR`): un `.asc` por fecha en la malla del DEM, para usar
precipitación/temperatura observada en lugar del valor de estación con
lapse rate. Ver la validación para el pipeline CR2MET. El aporte de cuenca
se rutea a un hidrograma con `--route-k DÍAS` (reservorio lineal).

El **acople con rainflow** (`coupling/`) cierra el balance hídrico: snowmelt
entrega el aporte líquido y GR4J produce el caudal. En el Río Choapa en
Cuncumén rescata a GR4J de inútil (NSE < 0 con precipitación cruda) a útil
(NSE +0.22). Ver [`coupling/README.md`](coupling/README.md).

## Roadmap (v0.12)

- Forzante sinóptico WRF real (cuando haya salidas <1 km descargables):
  el downscaling topográfico (v0.11) mostró que el detalle de terreno solo
  aporta marginalmente; el cuello de botella es la representación sinóptica
  del forzante, que WRF sí mejoraría.
- Sublimación con resistencia aerodinámica explícita y balance multi-año.
- Publicación: crates.io (core/cli) y wheel PyPI (bindings).

## Licencia

MIT OR Apache-2.0
