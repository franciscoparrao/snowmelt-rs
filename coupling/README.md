# Acople snowmelt-rs → rainflow

Acople operativo de los dos motores del ecosistema: **snowmelt-rs** resuelve
la fase nival (acumulación, balance de energía, ablación) y entrega el
**aporte líquido** diario (lluvia + derretimiento); **rainflow** (GR4J)
cierra el balance suelo-escorrentía (evapotranspiración, almacenamiento,
ruteo) y produce el caudal.

Este crate está **excluido del workspace** de snowmelt-rs: depende de un
checkout de rainflow adyacente (`../rainflow`) y es opt-in, de modo que el
build estándar de snowmelt-rs no lo necesita.

## Qué hace

Sobre el Río Choapa en Cuncumén (CAMELS-CL 4703002), compara en
split-sample (DDS, 3000 evaluaciones, warm-up 365 d) GR4J alimentado por:

- **precipitación cruda** (sin nieve) — el baseline ingenuo, y
- **el aporte líquido de snowmelt** (lluvia + derretimiento) — el acoplado.

## Resultado (NSE de validación)

| forzante de GR4J | A→B | B→A |
|---|---|---|
| precipitación cruda (sin nieve) | −0.378 | −0.162 |
| **aporte snowmelt (nieve)** | **0.223** | **0.228** |

GR4J con precipitación cruda es **inútil** en esta cuenca nival (NSE < 0):
trata la precipitación de invierno como caudal inmediato cuando en realidad
se acumula como nieve y sale en primavera. El aporte de snowmelt mueve esa
agua a la temporada de deshielo y rescata el modelo a NSE positivo —
consistente con el benchmark de rainflow (GR4J sin rutina de nieve:
val NSE 0.04 → −0.34 en esta misma cuenca).

**Contexto honesto**: el HBV+snow de rainflow, con su rutina de nieve y
**todos** sus parámetros calibrados al caudal, alcanza val NSE 0.34–0.62
(lumped) a 0.63–0.76 (bandas + lapse ajustado). El acople snowmelt→GR4J
usa parámetros nivales **físicos y fijos** (calibrados de forma
independiente contra MODIS, no contra el hidrograma), así que tiene menos
grados de libertad para ajustar el caudal; a cambio, la fase nival no se
sobre-ajusta al hidrograma y es transferible. La señal nival de snowmelt
explica por sí sola el 88% de la varianza del ciclo anual del caudal
(ver `validation/choapa-cuncumen/`).

## Ejecutar

Requiere `rainflow` en `../rainflow` y los datos de la cuenca generados:

```bash
# desde el raíz de snowmelt-rs:
python3 validation/choapa-cuncumen/build_catchment.py 12
cargo run --manifest-path coupling/Cargo.toml --release
```
