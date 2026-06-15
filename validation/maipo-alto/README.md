# Validación: cuenca alta del Maipo (El Yeso), temporada 2019

Primera validación de snowmelt-rs contra cobertura nival MODIS en una
cuenca andina real, con el modelo **sin calibrar** (parámetros default).

## Dominio y datos

| Ítem | Fuente | Detalle |
|---|---|---|
| Dominio | — | 70.35–69.95°W, 33.85–33.50°S (~37×39 km), UTM 19S, 200 m (197×188 celdas, 913–5304 m, media 2731 m) |
| DEM | Copernicus GLO-30 (Planetary Computer STAC) | 2 tiles, bilinear a 200 m |
| Temperatura | ERA5 vía Open-Meteo archive | Punto central del box, diaria, celda a 3117 m (= z_ref) |
| Precipitación | CR2MET v2.5 (local) | Media diaria del box, 0.05° |
| Observación | MOD10A1 v6.1 NDSI_Snow_Cover (Planetary Computer) | Binarizada NDSI ≥ 40 = nieve; nubes/fill = NODATA |
| Periodo | — | 2019-04-01 a 2019-12-31 (275 días) |

Fechas MODIS elegidas por baja nubosidad (junio quedó fuera: el mejor día
tuvo 57% de celdas válidas): 07-15 (93%), 08-08 (81%), 09-08 (97%),
10-15 (96%), 11-10 (99%).

## Configuración del modelo

Balance de energía + albedo dinámico, todo en defaults:

```bash
snowmelt --dem data/dem.asc --forcing data/forcing.csv \
  --out-dir out --z-ref 3117 --energy-balance --latitude -33.675 \
  --albedo-tau 6 --cover-threshold 10 \
  --snapshot-dates 2019-07-15,2019-08-08,2019-09-08,2019-10-15,2019-11-10
```

Celda cubierta si SWE ≥ 10 mm (comparable a detección NDSI ≥ 0.4).

## Resultados (v0.6.0)

Balance de energía + albedo dinámico. Dos configuraciones:

**Sin calibrar** (defaults; `--albedo-tau 6`):

```
par                        accuracy precision    recall        F1    bias
cover_2019-07-15             0.9028    0.9162    0.9202    0.9182   1.004
cover_2019-08-08             0.8882    0.8526    0.9145    0.8825   1.073
cover_2019-09-08             0.7518    0.6479    0.9999    0.7863   1.543
cover_2019-10-15             0.7806    0.9859    0.5181    0.6793   0.525
cover_2019-11-10             0.9294    0.6469    0.3824    0.4806   0.591
TOTAL                        0.8493    0.8084    0.8213    0.8148   1.016
```

**Calibrado** (`--albedo-tau 9 --albedo-min 0.4`, vía `calibrate.py`):

```
par                        accuracy precision    recall        F1    bias
cover_2019-07-15             0.9024    0.9001    0.9396    0.9194   1.044
cover_2019-08-08             0.8762    0.8189    0.9376    0.8742   1.145
cover_2019-09-08             0.7518    0.6479    0.9999    0.7863   1.543
cover_2019-10-15             0.8568    0.9625    0.7083    0.8161   0.736
cover_2019-11-10             0.8866    0.4181    0.8353    0.5573   1.998
TOTAL                        0.8538    0.7765    0.8956    0.8318   1.153
```

**Lectura**: pleno invierno excelente (F1 0.92 en julio). La calibración
sube el F1 agregado 0.815 → 0.832 (accuracy 85.4%, recall 0.82 → 0.90):
un albedo de nieve vieja que decae más lento (τ 6 → 9 días) frena la
ablación de primavera y rescata octubre (F1 0.68 → 0.82, bias 0.53 →
0.74). El grid search (4×3×4 combos sobre τ, α_min, nubosidad) confirma
que la nubosidad efectiva no mejora el F1 agregado: atenúa la onda corta
de forma uniforme y el óptimo queda en cielo despejado.

**Sesgo estructural de septiembre**: el bias 1.54 / recall 1.0 de
septiembre es **invariante a todos los parámetros del grid** — el modelo
mantiene nieve en elevaciones bajas donde MODIS ya no la ve.

## Forzante de precipitación distribuida (v0.7)

Hipótesis: el exceso de nieve de septiembre en los valles viene de aplicar
la precipitación **uniforme** del box a todas las elevaciones. Se probaron
tres fuentes de precipitación (todas con EB, `τ=9`, `α_min=0.4`):

| Precipitación | F1 agreg. | bias agreg. | sep F1 / bias |
|---|---|---|---|
| Uniforme (media del box) | 0.832 | 1.15 | 0.79 / 1.54 |
| CR2MET 0.05° distribuida (bilinear) | 0.831 | 1.17 | 0.79 / 1.55 |
| CR2MET + downscaling orográfico (γ=8e-4) | 0.819 | 1.19 | 0.79 / 1.53 |
| Gradiente orográfico lineal (`--precip-gradient 1e-3`) | **0.834** | **1.07** | 0.86 / 1.26 |

**Hallazgos**:

1. **CR2MET distribuida ≈ uniforme.** A 0.05° (~5 km) sobre una cuenca de
   ~37 km, CR2MET resuelve pocas celdas y su gradiente vertical efectivo es
   débil; regrillarla no cambia la cobertura.
2. **El downscaling orográfico de subgrilla** (realce por la anomalía de
   elevación respecto a lo que "ve" el píxel CR2MET, método tipo Liston &
   Elder) con γ físico (8e-4 m⁻¹) tampoco mueve septiembre: las anomalías de
   subgrilla son moderadas y el realce introduce ruido en las cumbres
   (noviembre empeora).
3. **El gradiente lineal global** sí mejora el sesgo agregado (bias 1.15 →
   1.07) y septiembre (F1 0.79 → 0.86), pero **por sobre-corrección**:
   referido a `z_ref = 3117 m`, anula la precipitación bajo ~1500 m. Funciona
   contra MODIS pero no es físicamente defendible como acumulación.

**Conclusión**: para esta cuenca pequeña la precipitación distribuida
disponible (CR2MET 0.05°) no resuelve el sesgo de elevaciones bajas; el
cuello de botella es la resolución del producto, no el modelo. La
infraestructura de forzantes distribuidas (`--precip-grids`,
`--temp-grids`) queda lista para fuentes de mayor resolución y para
temperatura distribuida real (ERA5 multi-celda), pendiente para una
próxima iteración.

## Caveats

- Modelo sin calibrar; un solo punto de temperatura (lapse rate fijo
  −6.5 °C/km); precipitación media del box sin gradiente orográfico.
- LW y SW de cielo despejado (los días de tormenta reciben radiación
  sobreestimada).
- MODIS binarizado a 500 m remuestreado a 200 m (nearest); NDSI ≥ 40 como
  verdad de terreno tiene sus propios errores en terreno escarpado.

## Reproducir

```bash
python3 fetch_data.py all      # DEM + forzantes + MODIS + grillas precip (≈3 min)
# precipitación distribuida: agrega --precip-grids data/grids al comando snowmelt
# corrida calibrada:
snowmelt --dem data/dem.asc --forcing data/forcing.csv --out-dir out \
  --z-ref 3117 --energy-balance --latitude -33.675 \
  --albedo-tau 9 --albedo-min 0.4 --cover-threshold 10 \
  --snapshot-dates 2019-07-15,2019-08-08,2019-09-08,2019-10-15,2019-11-10
snowmelt-validate out/cover_2019-07-15.asc:data/modis_2019-07-15.asc ...
# o el grid search completo:
python3 calibrate.py --top 8
```
