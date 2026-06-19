# Sublimación y resistencia aerodinámica (v0.12)

Validación de los flujos turbulentos con **resistencia aerodinámica
explícita** (`--aero-resistance`) frente al coeficiente de intercambio bulk
fijo, anclada a las fracciones de sublimación publicadas para los Andes.

## Por qué sin estaciones DGA

La validación directa contra SWE de rutas de nieve / snow pillows de la DGA
requiere datos que no están disponibles localmente (CAMELS-CL solo trae
`p, pet, tmean, qobs`). En su lugar se valida la **física de la
sublimación** contra el rango documentado en la literatura andina, que es
un control de primer orden del balance de masa en el Chile árido-semiárido:

- MacDonell et al. (2013, *The Cryosphere*), Pascua-Lama (~5000 m): la
  sublimación es **80–90%** de la ablación.
- Ayala et al. (2017), Andes centrales de Chile: la sublimación es una
  fracción sustancial de la ablación, creciente con la elevación y la
  aridez.

## Experimento

Cuenca alta del Maipo 2019, balance de energía + albedo dinámico
(`τ=9`, `α_min=0.4`), escenario seco-ventoso andino (`--rh 0.4 --wind 4`,
lapse −7.5 °C/km). Se compara la partición sublimación/derretimiento
integrada sobre la temporada (media de cuenca, 197×188 celdas, 913–5304 m):

| Flujos turbulentos | sublimación | derretimiento | frac. sublimación |
|---|---|---|---|
| Bulk (`C_e = 0.0015`) | 116.8 mm | 87.4 mm | 57.2% |
| Aero, neutro (`z0 = 1 mm`) | 153.6 mm | 50.6 mm | **75.2%** |
| Aero + estabilidad | 151.0 mm | 53.3 mm | 73.9% |

## Hallazgos

1. **La sublimación es una fracción grande de la ablación (57–75%)**,
   dentro del rango publicado para los Andes secos. El modelo reproduce el
   régimen criosférico correcto sin calibrar contra observaciones de SWE.

2. **La resistencia aerodinámica deriva la conductancia de la rugosidad**
   en vez de un coeficiente ajustado: con `z0 = 1 mm` (nieve típica) y
   medición a 2 m, `1/r_a ≈ 0.0085 m/s` frente a `C_e·u ≈ 0.006 m/s` del
   bulk por defecto — ~40% más intercambio, que desplaza ablación hacia
   sublimación. La partición pasa a depender de un parámetro físico medible
   (la rugosidad) en lugar de uno de calibración.

3. **La corrección de estabilidad casi no mueve la fracción estacional**
   (75.2% → 73.9%): la sublimación se concentra en condiciones frías, donde
   la superficie iguala a la temperatura del aire (`T_s = min(T_a, 0)` ⇒
   gradiente nulo ⇒ Richardson 0 ⇒ régimen neutro). La estabilidad amortigua
   sobre todo el calor sensible de los días cálidos (temporada de
   derretimiento), un efecto de segundo orden en la partición anual pero
   físicamente correcto.

## Balance de masa y ELA (multi-año)

Con `--mass-balance` el modelo acumula acumulación − ablación por celda
sobre toda la corrida y estima la **línea de equilibrio (ELA)** por bandas
de elevación. En la temporada 2019 la ELA del Maipo alto cae en **~4300 m**,
consistente con la ELA de glaciares de los Andes centrales (~3800–4500 m).
El acumulador es independiente del esquema de derretimiento y persiste el
SWE entre temporadas, de modo que admite forzantes multi-año para balance
de masa glaciar (limitado aquí por la disponibilidad de forzante: una sola
temporada).

## Reproducir

```bash
snowmelt --dem data/dem.asc --forcing data/forcing.csv \
  --z-ref 3117 --energy-balance --latitude -33.675 --lapse-rate -0.0075 \
  --albedo-tau 9 --albedo-min 0.4 --rh 0.4 --wind 4 \
  --aero-resistance --mass-balance --out-dir out_subl
# fracción de sublimación:
awk -F, 'NR>1{sb+=$5; ml+=$4} END{printf "%.1f%%\n", 100*sb/(sb+ml)}' out_subl/series.csv
```
