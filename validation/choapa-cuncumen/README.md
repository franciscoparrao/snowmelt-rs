# Validación de caudal: Río Choapa en Cuncumén (CAMELS-CL 4703002)

Primera validación de la **componente hidrológica** de snowmelt-rs: del
aporte de deshielo al caudal observado, sobre una cuenca nival andina con
38 años de registro.

## Cuenca y datos

| Ítem | Valor |
|---|---|
| Gauge | CAMELS-CL 4703002, Río Choapa en Cuncumén (Norte Chico) |
| Área | 1132 km²; régimen nival, swe_ratio 0.52 |
| Elevación | 1158–5054 m (media 3142 m) |
| Forzante | CAMELS-CL: `tmean`, `p` (CR2MET), media de cuenca, diario |
| Observación | `qobs` (mm/día), 1980–2016, 13 392 días, 99% cobertura |
| Geometría | 12 bandas de elevación de igual área desde la curva hipsométrica DEM-derivada (Copernicus GLO-30) |

El modelo corre sobre un pseudo-DEM 1×12 (una celda por banda; `cellsize`
grande → slopes ≈ 0, radiación horizontal sin aspecto, lapse rate
distribuye la temperatura). Balance de energía, albedo dinámico (`τ=9`,
`α_min=0.4`), lapse −7.5 °C/km (de la calibración del Maipo).

## Método

snowmelt produce el aporte vertical (lluvia + derretimiento); un
**reservorio lineal** (`snowmelt-core::routing`, `k=90 d`) le da el
retardo/recesión de la cuenca. Como snowmelt **no** es un balance
lluvia-escorrentía (no resta ET ni añade flujo base), se valida la
**componente de deshielo** contra la firma estacional del caudal, y se
contrasta con el aporte total.

## Resultados (1980–2016, k=90 d)

| señal | corr. diaria | NSE (forma) | corr. ciclo anual |
|---|---|---|---|
| **deshielo** | **0.813** | **0.659** | **0.879** |
| lluvia + deshielo | 0.692 | 0.465 | 0.603 |

Régimen nival simulado (correcto): SWE acumula abril–agosto (peak 116 mm),
derrite septiembre–octubre. El **deshielo ruteado reproduce la firma
estacional del caudal** — ciclo anual r=0.88, peak simulado mes 10 vs
observado mes 11 (un mes de adelanto, típico de un reservorio único frente
al deshielo escalonado por elevación).

**El aporte total es peor** (NSE 0.47): incluye lluvia invernal que va
directa al caudal, cuando en la cuenca real esa agua pasa por el suelo
(ET, recarga, flujo base). Esa brecha es exactamente lo que aporta un
modelo lluvia-escorrentía — y la razón de la interfaz hacia **rainflow**.

## Interfaz hacia rainflow

`out/series.csv` (columnas `melt_mm`, `rain_mm`, `runoff_mm`) es el
forzante de aporte para [rainflow](https://github.com/franciscoparrao/rainflow):
snowmelt resuelve la fase nival (acumulación, ablación, SWE) y entrega el
agua líquida disponible; rainflow (GR4J/HBV) cierra el balance
suelo-escorrentía y el ruteo. Esto reemplaza la rutina de nieve interna de
rainflow por un modelo de energía distribuido por elevación.

Comparación con los benchmarks HBV de CAMELS-CL en esta cuenca (README de
rainflow): HBV+snow lumped alcanza NSE val 0.31–0.62; el deshielo de
snowmelt explica el 88% de la varianza del ciclo anual del caudal por sí
solo, antes de cualquier balance — coherente con el swe_ratio 0.52
reportado y consistente como insumo nival para el acople.

## Reproducir

```bash
python3 build_catchment.py 12        # bandas + forcing + qobs desde CAMELS-CL
snowmelt --dem data/bands_dem.asc --forcing data/forcing.csv --out-dir out \
  --z-ref 3142 --energy-balance --latitude -31.95 \
  --albedo-tau 9 --albedo-min 0.4 --lapse-rate -0.0075 --route-k 90
python3 validate_flow.py --k 90
```

## Caveats

- No es balance lluvia-escorrentía: el caudal absoluto requiere rainflow.
- Forzante de cuenca media (tmean/p) distribuido solo por elevación; sin
  variabilidad horizontal ni gradiente orográfico de precipitación.
- Reservorio lineal de un parámetro (k) sin calibración formal multi-objetivo.
- CR2MET subcaptura la precipitación de alta cordillera (sesgo conocido).
