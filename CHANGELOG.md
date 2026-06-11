# Changelog

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
