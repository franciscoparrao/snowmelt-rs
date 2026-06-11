# snowmelt-rs — Modelo de derretimiento nival / glaciar (Rust)

> **Estado:** v0.3 implementado (2026-06-11): albedo dinámico por edad (AlbedoDecay),
> sombreado por horizonte (SurtGIS horizon_angles + solar_radiation_shadowed),
> bindings Python (snowmelt-python, PyO3+numpy abi3). v0.2 (2026-06-10): ETI
> (Pellicciotti 2005) con radiación de SurtGIS, gradiente orográfico. v0.1: grado-día
> distribuido + CLI. Pendiente v0.4: balance de energía completo, validación
> MODIS/DGA, interfaz rainflow. Creado 2026-06-10.
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
