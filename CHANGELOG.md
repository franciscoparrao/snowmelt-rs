# Changelog

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
