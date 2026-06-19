# snowmelt-rs — Modelo de derretimiento nival / glaciar (Rust)

> **Estado:** v0.12 implementado (2026-06-18): **resistencia aerodinámica + balance
> multi-año**. snowmelt-core::energy gana AeroResistance (conductancia 1/r_a desde
> rugosidad z0m/z0h + estabilidad Richardson bulk, reemplaza coef. bulk fijo); nuevo
> módulo balance (MassBalance acumulación−ablación por celda + ELA). CLI --aero-resistance,
> --mass-balance; PyParams sincronizado. Validación sublimación anclada a literatura
> (sin datos DGA): 57–75% de ablación en Maipo seco-ventoso, dentro del rango andino
> publicado; ELA ~4300 m. v0.11: **downscaling topográfico de forzantes** (downscale +
> terrain, MicroMet/Liston-Elder): T con curvatura (cold-air pooling), viento por terreno,
> precip elevación+barlovento. Vía honesta a sub-km SIN WRF. Re-validación MODIS Maipo
> (200 m): mejora marginal solo por curvatura-T (F1 0.834→0.836); cuello de botella es el
> forzante sinóptico, no el terreno.
> v0.10: **acople operativo snowmelt→rainflow** (coupling/, crate excluido, opt-in):
> GR4J con precip cruda inútil (NSE −0.38/−0.16) → con aporte snowmelt +0.22/+0.23 en
> Choapa-Cuncumén. Cierra las 3 conexiones del ecosistema (SurtGIS, MODIS, rainflow).
> v0.9: ruteo reservorio lineal + validación caudal CAMELS-CL (deshielo NSE forma 0.66,
> corr ciclo 0.88). v0.8: temp distribuida ERA5 + estudio forzantes (uniforme+lapse −7.5
> óptimo). v0.7: forzante distribuido. v0.6: física EB + calibración MODIS (F1 0.83).
> v0.5: validación MODIS. v0.4: EB. v0.3: albedo/horizonte/PyO3. v0.2: ETI+SurtGIS.
> v0.1: grado-día+CLI. Pendiente v0.13: WRF real, firn/hielo explícito, SWE DGA,
> publicación crates.io/PyPI. Creado 2026-06-10.
> Familia de motores Rust del autor: SurtGIS, Hydroflux, Smelt, Anvil, Cantus, Criterium.
> Doc madre: `~/proyectos/ideas-motores-rust.md` (idea G3).

## Qué es
Motor de balance nival sobre DEM: modelos grado-día y de balance de energía
para estimar SWE y aporte de deshielo, orientado a los Andes.

## El gap que llena
Eje criosférico ausente en tu familia. Aporte nival crítico para hidrología
andina (tus 15 cuencas). El campo es scripts/modelos sueltos (DHSVM, SnowModel).

## Alcance MVP (v0.1)
- [x] Modelo grado-día (degree-day) distribuido sobre DEM.
- [x] Acumulación/ablación de SWE; partición lluvia-nieve por temperatura.
- [x] Forzantes: temperatura/precipitación (series o grillas).
- [ ] (v0.2) Balance de energía; radiación solar (de SurtGIS terrain); línea de nieves desde RS.

## Arquitectura tentativa
- `snowmelt-core`: estados de SWE por celda, integración temporal.
- Targets: native (Rayon) + Python (PyO3) + CLI.
- Reusa radiación solar y DEM de SurtGIS.

## Validación / paridad numérica
Cross-check contra productos MODIS de cobertura nival y caudales DGA andinos.

## Venue objetivo
**Journal of Hydrology** o **The Cryosphere**.

## Conexiones con tu ecosistema
- **rainflow / Hydroflux**: aporte de deshielo como forzante.
- **SurtGIS**: radiación solar, DEM, derivados de terreno.
- **datacube-rs**: cobertura nival temporal desde Sentinel/MODIS.

## Próximos pasos al retomar
1. Implementar grado-día distribuido + partición lluvia-nieve.
2. Validar cobertura simulada contra MODIS en una cuenca andina.
3. Definir interfaz de aporte hacia rainflow.
