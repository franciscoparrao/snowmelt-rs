# Estudio de sensibilidad de forzantes — Maipo alto 2019

Síntesis del eje de validación (v0.6–v0.8): ¿qué forzante meteorológico
maximiza la habilidad del modelo de cobertura nival contra MODIS MOD10A1,
y dónde está el límite? Todas las corridas usan balance de energía,
albedo dinámico (`τ=9 d`, `α_min=0.4`), `z_ref=3117 m`, sobre 5 fechas
despejadas (jul–nov 2019). Métricas agregadas (172 435 celdas-fecha):

| # | Forzante | F1 | accuracy | bias | sep F1 / bias |
|---|---|---|---|---|---|
| 1 | Uniforme, lapse −6.5 °C/km (base) | 0.832 | 0.854 | 1.15 | 0.79 / 1.54 |
| 2 | **Uniforme, lapse −7.5 °C/km (operacional)** | **0.834** | **0.859** | **1.11** | 0.81 / 1.46 |
| 3 | Uniforme, lapse −8.5 °C/km | 0.833 | 0.857 | 1.08 | 0.83 / 1.39 |
| 4 | CR2MET 0.05° precip distribuida | 0.831 | 0.852 | 1.17 | 0.79 / 1.55 |
| 5 | CR2MET + downscaling orográfico (γ=8e-4) | 0.819 | 0.840 | 1.19 | 0.79 / 1.53 |
| 6 | Gradiente lineal precip (1e-3 m⁻¹) | 0.834 | — | 1.07 | 0.86 / 1.26 |
| 7 | ERA5 temperatura distribuida (lapse empírico) | 0.793 | 0.797 | 1.43 | 0.78 / 1.55 |

## Hallazgos

**1. El lapse rate es el control dominante** (filas 1–3). Endurecerlo de
−6.5 a −7.5 °C/km sube el calor en elevaciones medias, eleva la línea de
nieve simulada y es la única intervención que mejora F1 *y* bias *y*
septiembre a la vez. El óptimo (−7.5) coincide con el rango del lapse
empírico derivado de la nube de puntos ERA5 (media −6.6, rango diario
−4.0…−8.4 °C/km): el clima árido andino se aparta del −6.5 estándar.

**2. La precipitación distribuida de reanálisis no aporta** (filas 4–5).
CR2MET a 0.05° (~5 km) sobre una cuenca de ~37 km da esencialmente el
mismo campo que el promedio uniforme; el downscaling orográfico de
subgrilla con un γ físico no resuelve el contraste vertical y agrega
ruido en cumbres. Solo un gradiente lineal global (fila 6) mueve el
sesgo, y lo hace por **sobre-corrección** (anula la precipitación bajo
~1500 m), no por representar mejor la acumulación.

**3. La temperatura distribuida ERA5 empeora** (fila 7). El campo
multi-celda, más representativo, es ~1–2.5 °C más frío que el punto
central (sesgo creciente con la elevación), porque la celda ERA5 que
alimentaba el forzante uniforme era casualmente cálida. Es decir: **el
forzante de punto único acertaba por compensación de errores**. El
forzante distribuido lo expone, sobreestimando nieve (recall 0.96,
precision 0.67).

**4. El sesgo de septiembre es robusto.** Ninguna fuente distribuida lo
resuelve; solo responde al lapse rate (un parámetro escalar). Combinado
con que su recall es ~1.0 (el modelo nunca pierde nieve real, solo agrega
de más en la zona de transición), la causa más probable es una mezcla de
(a) incertidumbre de MODIS NDSI en laderas de transición/escarpadas y
(b) la resolución de los forzantes (5–25 km) frente al rango vertical de
4400 m de la cuenca.

## Conclusión operacional

Para esta cuenca y escala, el **forzante uniforme con lapse rate
calibrado (−7.5 °C/km)** es la configuración recomendada: máxima
habilidad (F1 0.834, accuracy 85.9%) con el menor número de supuestos.
Los forzantes distribuidos de reanálisis disponibles no superan esta
línea base; lo harían productos de temperatura/precipitación de
resolución < 1 km (p. ej. WRF dinámicamente downscaleado), pendiente
para trabajo futuro. La infraestructura de forzantes distribuidos
(`--precip-grids`, `--temp-grids`) y el pipeline de lapse empírico quedan
listos para cuando esos datos estén disponibles.
