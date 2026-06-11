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

## Resultados (v0.5.0, sin calibración)

```
par                        celdas      TP      FP      FN      TN  accuracy precision    recall        F1    bias
cover_2019-07-15            34435   19022    1910    1386   12117    0.9043    0.9088    0.9321    0.9203   1.026
cover_2019-08-08            29944   12719    2443    1023   13759    0.8843    0.8389    0.9256    0.8801   1.103
cover_2019-09-08            35793   16347    9118       0   10328    0.7453    0.6419    1.0000    0.7819   1.558
cover_2019-10-15            35602    9791     227    6175   19409    0.8202    0.9773    0.6132    0.7536   0.627
cover_2019-11-10            36661    1801    1363    1332   32165    0.9265    0.5692    0.5748    0.5720   1.010
TOTAL                      172435   59680   15061    9916   87778    0.8552    0.7985    0.8575    0.8270   1.074
```

**Lectura**: pleno invierno excelente (F1 0.92 en julio); inicio de
primavera sobreestima extensión (bias 1.56, recall 1.0 en septiembre:
derrite lento en elevaciones medias) y a mediados de octubre el sesgo se
invierte (bias 0.63: el pack simulado se agota más rápido que el
observado, consistente con albedo decayendo a α_min y LW de cielo
despejado). El cruce sugiere calibrar `t_cold_max`, `albedo_min`/τ y la
nubosidad efectiva.

## Caveats

- Modelo sin calibrar; un solo punto de temperatura (lapse rate fijo
  −6.5 °C/km); precipitación media del box sin gradiente orográfico.
- LW y SW de cielo despejado (los días de tormenta reciben radiación
  sobreestimada).
- MODIS binarizado a 500 m remuestreado a 200 m (nearest); NDSI ≥ 40 como
  verdad de terreno tiene sus propios errores en terreno escarpado.

## Reproducir

```bash
python3 fetch_data.py all     # DEM + forzantes + MODIS (≈2 min)
# luego el comando snowmelt de arriba, y:
snowmelt-validate out/cover_F.asc:data/modis_F.asc ...
```
