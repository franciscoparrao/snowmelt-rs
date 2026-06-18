# Changelog

## 0.11.0 — 2026-06-18

### Agregado
- **Downscaling topográfico de forzantes** (`snowmelt-core::downscale`,
  estilo MicroMet, Liston & Elder 2006): a partir de un valor escalar de
  temperatura y precipitación produce campos a la resolución del DEM con
  tres controles de terreno más allá del lapse rate — (1) temperatura con
  término de **curvatura** (cold-air pooling: enfría valles, templa
  cumbres), (2) **viento** por terreno `W = 1 + γ_s·Ω_s + γ_c·Ω_c`, y
  (3) precipitación con factor de elevación (Thornton 1997) y **realce
  orográfico a barlovento** `1 + γ_w·Ω_s`.
- `snowmelt-core::terrain`: derivados de terreno autocontenidos (sin
  dependencia GIS) — slope/aspect por Horn (1981) y curvatura normalizada
  Liston-Elder. 5 tests.
- `downscale`: `Downscaler` + `DownscaleParams`. 9 tests.
- CLI: `--downscale` con `--temp-curvature`, `--precip-elev-factor`,
  `--precip-windward`, `--wind-dir`, `--wind-slope-weight`,
  `--wind-curvature-weight`. Genera el forzante distribuido sin grillas
  externas (excluyente con `--temp-grids`/`--precip-grids`).

### Validación
- Re-validación MODIS en el Maipo alto 2019 (sub-km, 200 m, vs lapse
  −7.5 °C/km): el downscaling da una mejora **marginal y solo por la
  curvatura-temperatura** (F1 0.834 → **0.836**, accuracy 85.85 → 86.02%);
  la precipitación orográfica a barlovento y el factor de elevación **no
  ayudan**, consistente con el estudio nulo previo. El sesgo estructural de
  septiembre es invariante al detalle de terreno (responde solo al lapse
  rate). Confirma la hipótesis del estudio de forzantes: el cuello de
  botella es la representación sinóptica del forzante (que aportaría WRF),
  no el detalle topográfico que el DEM ya codifica. Detalle en
  `validation/maipo-alto/FORCING_SENSITIVITY.md`.

## 0.10.0 — 2026-06-15

### Agregado
- **Acople operativo snowmelt-rs → rainflow** (`coupling/`, crate excluido
  del workspace, opt-in): snowmelt resuelve la fase nival y entrega el aporte
  líquido (lluvia+derretimiento); rainflow (GR4J) cierra el balance
  suelo-escorrentía. Encadena ambos motores in-process (Rust→Rust).
- Pipeline: `build_catchment.py` ahora exporta también `balance.csv`
  (`precip_mm, pet_mm, qobs_mm`) como insumo del modelo lluvia-escorrentía.

### Resultado
- En el Río Choapa en Cuncumén (CAMELS-CL 4703002), GR4J con precipitación
  cruda es inútil (val NSE −0.38 / −0.16); alimentado por el aporte de
  snowmelt sube a **+0.22 / +0.23**. Confirma que la fase nival es decisiva
  y que la interfaz funciona. Los parámetros nivales son físicos y fijos
  (calibrados contra MODIS, no contra el hidrograma): más parsimonioso y
  transferible que el HBV+snow totalmente calibrado de rainflow.

## 0.9.0 — 2026-06-15

### Agregado
- **Ruteo por reservorio lineal** (`snowmelt-core::routing`): `LinearReservoir`
  (Nash, una caja) transforma el aporte de cuenca (lluvia+derretimiento) en
  un hidrograma con retardo/recesión, conservando masa exactamente; helper
  `depth_to_discharge` (mm → m³/s). 7 tests.
- CLI: `--route-k DÍAS` rutea el aporte medio de cuenca y agrega la columna
  `routed_mm` a `series.csv`.
- **Validación de caudal** (`validation/choapa-cuncumen/`) contra CAMELS-CL
  4703002 (Río Choapa en Cuncumén, 1132 km², nival, 38 años): cuenca
  representada por 12 bandas de elevación de igual área desde la curva
  hipsométrica DEM-derivada; pipeline reproducible (`build_catchment.py`,
  `validate_flow.py`).
- **Interfaz hacia rainflow**: `series.csv` documentado como forzante de
  aporte para el modelo lluvia-escorrentía.

### Validación
- El **deshielo ruteado** (k=90 d) reproduce la firma estacional del caudal
  observado: **corr. diaria 0.81, NSE de forma 0.66, corr. ciclo anual 0.88**
  (peak simulado mes 10 vs observado 11). Régimen nival correcto (SWE peak
  116 mm en agosto, ablación sep–oct).
- El aporte total (rain+melt) es peor (NSE 0.47): la lluvia invernal requiere
  el balance lluvia-escorrentía de rainflow — motivación de la interfaz.

## 0.8.0 — 2026-06-15

### Agregado
- **Pipeline de temperatura distribuida ERA5** (`fetch_data.py tempgrids`):
  malla de puntos ERA5 (Open-Meteo multi-celda), **lapse rate empírico**
  derivado por regresión T-vs-z diaria, reducción a nivel del mar,
  interpolación horizontal al DEM y re-extrapolación topográfica →
  `temp_FECHA.asc` para `--temp-grids`.
- **Estudio de sensibilidad de forzantes** (`validation/maipo-alto/FORCING_SENSITIVITY.md`):
  síntesis cuantitativa v0.6–v0.8 de 7 configuraciones de forzante.

### Hallazgos de validación
- **El lapse rate es el control dominante** del sesgo de cobertura;
  calibrarlo a −7.5 °C/km (dentro del rango empírico ERA5 −4…−8.4) da la
  mejor configuración: **F1 0.834, accuracy 85.9%, bias 1.11**.
- Temperatura distribuida ERA5 empeora (F1 0.79): el campo multi-celda es
  ~1–2.5 °C más frío que el punto central, revelando que el forzante
  uniforme acertaba por compensación de errores.
- Ningún forzante distribuido de reanálisis (CR2MET 0.05°, ERA5) supera al
  uniforme calibrado en esta cuenca de ~37 km; el límite es la resolución
  de los productos y la incertidumbre de MODIS NDSI en la zona de
  transición, no el modelo.

### Sin cambios de motor
- v0.8 es análisis + pipeline de datos; el soporte `--temp-grids` ya
  existía desde v0.7.

## 0.7.0 — 2026-06-14

### Agregado
- **Forzante distribuido por grillas diarias** en la CLI: `--precip-grids
  DIR` y `--temp-grids DIR` leen `precip_FECHA.asc` / `temp_FECHA.asc` (malla
  del DEM) y alimentan `Forcing::Distributed`; el lado faltante se construye
  por lapse rate (temperatura) o gradiente orográfico (precipitación). El
  CSV de forzantes aporta las fechas.
- **Pipeline de precipitación distribuida CR2MET** (`fetch_data.py grids`):
  regrillado bilinear de CR2MET 0.05° al DEM con **downscaling orográfico de
  subgrilla** (realce por anomalía de elevación, método tipo Liston & Elder).
- Test de equivalencia `Forcing::Uniform` ↔ `Forcing::Distributed`
  construido por lapse (valida la rama distribuida de la CLI).

### Validación
- Comparación de tres fuentes de precipitación en el Maipo alto. **Hallazgo
  principal**: CR2MET 0.05° distribuida ≈ uniforme (resolución de ~5 km
  insuficiente para una cuenca de ~37 km); el downscaling orográfico con γ
  físico no resuelve el sesgo de septiembre. El gradiente lineal mejora el
  bias agregado (1.15 → 1.07) pero por sobre-corrección. La infraestructura
  queda lista para fuentes de mayor resolución / temperatura distribuida.

## 0.6.0 — 2026-06-13

### Agregado
- **Física del balance de energía completa**:
  - *Nubosidad efectiva* (`EnergyBalanceParams.cloud_fraction`, flag
    `--cloud-fraction`): atenúa la onda corta `(1 − 0.75·N³)` y aumenta la
    emisividad atmosférica de onda larga `(1 + 0.22·N²)`.
  - *Calor de lluvia sobre nieve*: la lluvia a temperatura del aire aporta
    `c_w·P_liq·max(T,0)` al balance (advección).
  - *Sublimación con pérdida de masa*: el flujo latente negativo retira SWE
    (`L_s = 2.834 MJ/kg`); nueva salida `StepOutput.sublimation`,
    `StepSummary.mean_sublimation`, columna `sublimation_mm` en `series.csv`
    y clave `sublimation` en Python. Balance de masa:
    `Δswe = snowfall − melt − sublimation`.
- **Calibración contra MODIS** (`validation/maipo-alto/calibrate.py`): grid
  search sobre τ de albedo, α_min, nubosidad y cold content, evaluado con
  `snowmelt-validate`. En el Maipo alto sube el F1 agregado **0.815 →
  0.832** (accuracy 85.4%, recall 0.82 → 0.90) con `--albedo-tau 9
  --albedo-min 0.4`. Documenta el sesgo estructural de septiembre como
  limitación del forzante de temperatura de punto único.

### Cambiado
- `net_energy` → `energy_fluxes` (devuelve `(total, latente)` y recibe
  lluvia + `dt`); `StepOutput` y `StepSummary` ganan campos de sublimación
  (breaking para construcción por literal).

## 0.5.0 — 2026-06-11

### Agregado
- **Validación MODIS en cuenca andina real** (`validation/maipo-alto/`):
  cuenca alta del Maipo, temporada 2019, DEM Copernicus GLO-30 + ERA5
  (Open-Meteo) + CR2MET pr + MOD10A1 v6.1. Con defaults (sin calibrar):
  **accuracy 85.5%, F1 0.83, bias 1.07** sobre 5 fechas despejadas;
  julio F1 0.92. Script reproducible `fetch_data.py` + README con
  métricas y caveats.
- **`snowmelt-validate`**: binario de métricas de cobertura
  (confusión, accuracy, precision, recall, F1, bias) entre grillas .asc
  simuladas y observadas, por par y agregado.
- **Snapshots por fecha**: `--snapshot-dates d1,d2,...` escribe
  `swe_FECHA.asc` y `cover_FECHA.asc` (umbral `--cover-threshold`,
  default 10 mm SWE en la validación).
- `snowmelt-cli` expone lib interna (`asc`, `forcing_csv`, `solar`)
  compartida entre ambos binarios.

## 0.4.0 — 2026-06-11

### Agregado
- **Balance de energía completo** (`EnergyBalanceParams`, opcional en
  `DegreeDayParams.energy_balance`): onda corta neta `(1−α)·G`, onda
  larga con emisividad atmosférica de Brutsaert (1975) y superficie de
  nieve a `min(T_a, 0 °C)`, flujos turbulentos bulk (sensible y latente,
  con viento/HR parametrizados y densidad del aire desde presión por
  elevación), calor de suelo constante, y **cold content** por celda
  (J/m², cap `c_ice·SWE·t_cold_max`): la energía negativa enfría el pack
  y la positiva paga el déficit antes de derretir (L_f = 334 kJ/kg).
  Reusa la radiación de SurtGIS (con sombreado opcional) y el albedo
  constante o dinámico. Accessor `cold_content()`.
- CLI: `--energy-balance` + `--wind`, `--rh`, `--snow-emissivity`,
  `--exchange-coeff`, `--ground-heat`, `--t-cold-max`.
- Python: kwargs `energy_balance=True`, `wind`, `rh`, etc. y
  `SnowModel.cold_content()`.

### Simplificaciones documentadas
- Sin calor de lluvia sobre nieve ni pérdida de masa por sublimación;
  onda larga de cielo despejado (consistente con la SW de cielo
  despejado).

## 0.3.0 — 2026-06-11

### Agregado
- **Albedo dinámico por edad de la nieve** (`AlbedoDecay`, opcional en
  `DegreeDayParams.albedo_decay`): decaimiento exponencial
  `α(t) = α_min + (α_fresh − α_min)·exp(−t/τ)` con reinicio a fresco
  cuando la nevada del paso supera `refresh_swe_mm`. Accessors
  `SnowModel::albedo()` y `snow_age()`; `StepSummary.mean_albedo` y
  columna `albedo` en `series.csv`. Flags: `--albedo-tau`,
  `--albedo-fresh`, `--albedo-min`, `--albedo-refresh`.
- **Sombreado por horizonte topográfico** (cast shadows): precálculo de
  ángulos de horizonte (SurtGIS `horizon_angles`) y uso de
  `solar_radiation_shadowed`. Flags: `--horizon-shading`,
  `--horizon-radius`, `--horizon-directions`.
- **Bindings Python (PyO3 + numpy, abi3-py39)**: crate `snowmelt-python`
  con `Params`, `SnowModel` (`step_uniform`, `step_distributed`, `swe()`,
  `albedo()`, `snow_age()`); grillas numpy float64 con NaN = nodata.
  Compilar con `maturin develop -m crates/snowmelt-python/Cargo.toml`.

### Cambiado
- El paso del modelo se reorganizó en 4 pases paralelos (partición →
  edad/albedo → acumulación+melt → lluvia); sin cambio de resultados.
- `StepSummary` agrega `mean_albedo` (breaking para construcción por
  literal).

## 0.2.0 — 2026-06-10

### Agregado
- **Modelo ETI (enhanced temperature-index, Pellicciotti et al. 2005)**:
  nuevo parámetro `srf` (factor de radiación de onda corta) y `albedo`.
  Con `srf > 0`, el derretimiento es `ddf·(T − t_melt) + srf·(1 − albedo)·G`
  para `T > t_melt`. Con `srf = 0` (default) se mantiene el grado-día puro.
- **Radiación potencial desde el DEM (SurtGIS terrain)**: el CLI deriva
  slope/aspect (Horn) y calcula radiación de cielo despejado por día del
  año (`surtgis-algorithms`), con caché por día y relleno de bordes con
  radiación de terreno plano. Flags: `--srf`, `--albedo`, `--latitude`,
  `--transmittance`, `--linke-turbidity`.
- **Gradiente orográfico de precipitación**: `precip_gradient` (m⁻¹) en
  forzantes uniformes, `p(z) = p_ref·(1 + grad·(z − z_ref))` acotado a ≥ 0.
  Flag `--precip-gradient`.
- `SnowModel::step_radiation(forcing, Option<radiation>, dt_days)`: punto
  de entrada general con grilla de radiación (W m⁻² medios diarios).

### Cambiado
- `DegreeDayParams` agrega los campos `srf`, `albedo`, `precip_gradient`
  (breaking para construcción por literal; `Default` cubre el caso v0.1).
- `step_days`/`step` delegan en `step_radiation` (sin cambio de
  comportamiento con `srf = 0`).

## 0.1.0 — 2026-06-10

- MVP: modelo grado-día distribuido sobre DEM (`snowmelt-core`) + CLI
  (`snowmelt`) con DEM ESRI ASCII Grid y forzantes CSV diarias.
- Acumulación/ablación de SWE, partición lluvia-nieve lineal, lapse rate,
  forzantes uniformes o distribuidas, nodata como NaN, Rayon por celda.
