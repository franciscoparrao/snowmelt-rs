#!/usr/bin/env python3
"""Construye el caso de validación de caudal para el Río Choapa en Cuncumén
(CAMELS-CL gauge 4703002, 1132 km², nival, swe_ratio 0.52, 1158–5054 m).

Genera en data/:
  - bands_dem.asc   Pseudo-DEM 1×N de bandas de elevación de igual área
                    leídas de la curva hipsométrica DEM-derivada de
                    CAMELS-CL. cellsize grande (100 km) → slopes ≈ 0 entre
                    bandas, de modo que la radiación queda horizontal
                    (sin aspecto) y el lapse rate distribuye la temperatura.
  - forcing.csv     date,temp_c,precip_mm desde el forzante de cuenca
                    CAMELS-CL (tmean, p).
  - qobs.csv        date,qobs_mm para la validación (mm/día; NA = gap).

Uso: python3 build_catchment.py [n_bands]
"""

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
CAMELS = os.path.expanduser("~/proyectos/rainflow/data/camels-cl")
GAUGE = "4703002"
CELLSIZE = 100_000.0  # m, grande a propósito (slopes ≈ 0)


def equal_area_bands(hyps_path, n):
    """n elevaciones de bandas de igual área desde la curva hipsométrica
    (area_fraction, elevation_m), evaluadas en los puntos medios de cada
    1/n de área."""
    af, ev = [], []
    for row in csv.DictReader(open(hyps_path)):
        af.append(float(row["area_fraction"]))
        ev.append(float(row["elevation_m"]))
    af, ev = np.array(af), np.array(ev)
    mids = (np.arange(n) + 0.5) / n
    return np.interp(mids, af, ev)


def write_bands_dem(path, elevations):
    n = len(elevations)
    with open(path, "w") as f:
        f.write(f"ncols {n}\nnrows 1\n")
        f.write(f"xllcorner 0\nyllcorner 0\ncellsize {CELLSIZE}\n")
        f.write("NODATA_value -9999\n")
        f.write(" ".join(f"{e:.2f}" for e in elevations) + "\n")


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    os.makedirs(DATA, exist_ok=True)

    elevations = equal_area_bands(os.path.join(CAMELS, f"{GAUGE}_hypsometry.csv"), n)
    z_ref = float(np.mean(elevations))
    write_bands_dem(os.path.join(DATA, "bands_dem.asc"), elevations)
    print(f"bands_dem.asc: {n} bandas de igual área "
          f"{elevations.min():.0f}–{elevations.max():.0f} m, media {z_ref:.0f} m")

    rows = list(csv.DictReader(open(os.path.join(CAMELS, f"{GAUGE}.csv"))))
    with open(os.path.join(DATA, "forcing.csv"), "w") as ff, \
         open(os.path.join(DATA, "qobs.csv"), "w") as fq:
        ff.write("date,temp_c,precip_mm\n")
        fq.write("date,qobs_mm\n")
        n_q = 0
        for r in rows:
            t = r.get("tmean", "")
            p = r.get("p", "0")
            if t in ("", "NA"):
                continue
            ff.write(f"{r['date']},{float(t):.2f},{max(float(p), 0.0):.2f}\n")
            q = r.get("qobs", "NA")
            fq.write(f"{r['date']},{q if q not in ('', 'NA') else 'NA'}\n")
            if q not in ("", "NA"):
                n_q += 1
        print(f"forcing.csv + qobs.csv: {len(rows)} días, {n_q} con qobs")
    print(f"\nz_ref sugerido: {z_ref:.0f} m | latitud Cuncumén ≈ -31.95")


if __name__ == "__main__":
    main()
