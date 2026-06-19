# Validación de SWE contra Andes-SR (Cortés & Margulis)

Validación del **SWE** de snowmelt-rs (no solo cobertura) contra la
**Andes Snow Reanalysis** (Andes-SR), la reanálisis de SWE de referencia
para los Andes centrales: SWE diario, ~180 m, **años de agua 1985–2015**,
físicamente basada y restringida por cobertura nival Landsat
(Cortés & Margulis, GRL 2017; Saavedra et al., Front. Earth Sci. 2020).

Es la pieza que cierra el eje de validación que quedó bloqueado por falta de
SWE observado de la DGA: Andes-SR cubre Maipo y Choapa y da el campo de SWE
continuo contra el cual medir RMSE, sesgo, correlación y KGE.

## Estado

| Pieza | Estado |
|---|---|
| Métricas de SWE continuo (`snowmelt-core::metrics`) | ✅ listo, testeado |
| `snowmelt-validate --mode swe` (RMSE/MBE/MAE/corr/KGE) | ✅ listo |
| Ingesta HDF5 Andes-SR → `.asc` (`fetch_andessr.py`) | ⏳ pendiente de la estructura del HDF5 |
| Corrida del modelo en WY ≤2015 + comparación | ⏳ pendiente de los datos |

## Restricciones (importantes)

1. **Andes-SR termina en el año de agua 2015.** La validación corre el
   modelo en un año ≤2015 (p. ej. WY2010), no en 2019 como la de MODIS. El
   forzante (ERA5 vía Open-Meteo archive, CR2MET v2.5) está disponible para
   esos años.
2. **El dataset está en un Box de UCLA** (web app con sign-up); no se puede
   descargar de forma programática desde el pipeline. La descarga es manual.

## Qué descargar (tú)

1. Entrar a <https://ucla.box.com/v/ANDES-SWE-REANALYSIS> (sign-up gratis si
   hace falta) y abrir el README del dataset.
2. Descargar el **tile que cubre el box del Maipo alto** para **un año de
   agua ≤ 2015**:
   - Box objetivo (igual que `validation/maipo-alto/`):
     **70.35–69.95 °W, 33.85–33.50 °S**.
   - Dejar el/los `.h5` en `validation/swe-andes-sr/data/`.
3. **Compartir la estructura del archivo** para escribir la ingesta exacta:

   ```bash
   h5ls -r validation/swe-andes-sr/data/<archivo>.h5
   # o, si tienes Python:
   python3 -c "import h5py,sys; h5py.File(sys.argv[1]).visititems(lambda n,o: print(n, getattr(o,'shape',''), getattr(o,'dtype','')))" validation/swe-andes-sr/data/<archivo>.h5
   ```

   Pegar esa salida en el chat. Con los nombres reales de los datasets
   (SWE, lat/lon o el grid, DEM) finalizo `fetch_andessr.py`.

## Pipeline planificado (una vez con los datos)

1. `fetch_andessr.py`: lee el HDF5, recorta al box del Maipo, **reproyecta y
   regrilla** el SWE diario a la malla del DEM del modelo (UTM 19S, 200 m;
   ver `validation/maipo-alto/data/dem.asc`) y escribe `swe_FECHA.asc`.
2. Forzante del año de agua elegido (ERA5 + CR2MET) con el `fetch_data.py`
   del Maipo, parametrizado al año.
3. Correr snowmelt-rs (EB + albedo dinámico + lapse −7.5, idealmente con
   `--aero-resistance`) con `--snapshot-dates` en fechas clave del año
   (acumulación máxima, deshielo) y `--mass-balance`.
4. Comparar SWE:

   ```bash
   snowmelt-validate --mode swe \
     out/swe_2010-09-01.asc:data/andessr_2010-09-01.asc ...
   ```

   Métricas: RMSE, sesgo (MBE), MAE, correlación y KGE sobre celdas
   co-válidas, por fecha y agregado; además la serie de SWE medio de cuenca
   modelado vs Andes-SR (acumulación, SWE pico, fecha de melt-out).

## Referencias

- Cortés, G. & Margulis, S. (2017), *Impacts of El Niño and La Niña on
  interannual snow accumulation in the Andes*, GRL.
- Saavedra, F. et al. (2020), Front. Earth Sci. 8:328.
- Datos: Margulis Research Group, <https://margulis-group.github.io/data/>.
