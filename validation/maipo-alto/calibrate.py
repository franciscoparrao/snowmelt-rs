#!/usr/bin/env python3
"""Calibración por grid search de snowmelt-rs contra cobertura MODIS.

Recorre combinaciones de parámetros del balance de energía / albedo,
corre el modelo, evalúa con snowmelt-validate y reporta las mejores por
F1 agregado. Pensado para diagnosticar los sesgos estacionales de la
validación del Maipo alto (sept sobreestima, oct/nov subestiman).

Uso: python3 calibrate.py [--top N] [--quick]
"""

import argparse
import itertools
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SNOWMELT = os.path.join(ROOT, "target", "release", "snowmelt")
VALIDATE = os.path.join(ROOT, "target", "release", "snowmelt-validate")
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "out_cal")

DATES = ["2019-07-15", "2019-08-08", "2019-09-08", "2019-10-15", "2019-11-10"]
BASE = [
    "--dem", os.path.join(DATA, "dem.asc"),
    "--forcing", os.path.join(DATA, "forcing.csv"),
    "--z-ref", "3117",
    "--energy-balance", "--latitude", "-33.675",
    "--cover-threshold", "10",
    "--snapshot-dates", ",".join(DATES),
]

# Espacio de búsqueda. Física esperada:
#  - cloud_fraction: ↓SW dominante en primavera → frena melt tardío (oct/nov).
#  - albedo_tau / albedo_min: albedo alto más tiempo → ablación más gradual.
#  - t_cold_max: más cold content retrasa el inicio del melt.
GRID = {
    "albedo-tau": [4.0, 6.0, 9.0, 14.0],
    "albedo-min": [0.4, 0.5, 0.6],
    "cloud-fraction": [0.0, 0.2, 0.4, 0.5],
    "t-cold-max": [10.0],
}
GRID_QUICK = {
    "albedo-tau": [6.0, 12.0],
    "albedo-min": [0.4, 0.55],
    "cloud-fraction": [0.0, 0.3],
    "t-cold-max": [10.0],
}

ROW_RE = re.compile(r"^TOTAL\s+(.+)$", re.M)
PAIR_RE = re.compile(r"^(cover_\S+)\s+(.+)$", re.M)


def run_combo(combo):
    args = [SNOWMELT, *BASE, "--out-dir", OUT]
    for k, v in combo.items():
        args += [f"--{k}", str(v)]
    subprocess.run(args, check=True, capture_output=True)
    pairs = [
        f"{OUT}/cover_{d}.asc:{DATA}/modis_{d}.asc" for d in DATES
    ]
    res = subprocess.run([VALIDATE, *pairs], check=True, capture_output=True, text=True)
    return res.stdout


def parse_metrics(stdout):
    # columnas: celdas TP FP FN TN accuracy precision recall F1 bias
    m = ROW_RE.search(stdout)
    cols = m.group(1).split()
    f1, bias = float(cols[-2]), float(cols[-1])
    acc = float(cols[-5])
    per_date = {}
    for name, rest in PAIR_RE.findall(stdout):
        c = rest.split()
        date = name.replace("cover_", "")
        per_date[date] = (float(c[-2]), float(c[-1]))  # F1, bias
    return acc, f1, bias, per_date


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    grid = GRID_QUICK if args.quick else GRID
    keys = list(grid)
    combos = [dict(zip(keys, vals)) for vals in itertools.product(*grid.values())]
    print(f"Evaluando {len(combos)} combinaciones...", file=sys.stderr)

    results = []
    for i, combo in enumerate(combos, 1):
        try:
            acc, f1, bias, per_date = parse_metrics(run_combo(combo))
        except subprocess.CalledProcessError as e:
            print(f"  [{i}] error: {e.stderr.decode()[:200]}", file=sys.stderr)
            continue
        results.append((f1, acc, bias, combo, per_date))
        print(f"  [{i}/{len(combos)}] F1={f1:.4f} acc={acc:.4f} bias={bias:.3f} {combo}",
              file=sys.stderr)

    results.sort(reverse=True, key=lambda r: r[0])
    print("\n=== TOP por F1 agregado ===")
    for f1, acc, bias, combo, per_date in results[: args.top]:
        params = " ".join(f"{k}={v}" for k, v in combo.items())
        print(f"F1={f1:.4f} acc={acc:.4f} bias={bias:.3f} | {params}")
        season = " ".join(f"{d[5:]}:F1={v[0]:.2f},b={v[1]:.2f}" for d, v in per_date.items())
        print(f"    {season}")


if __name__ == "__main__":
    main()
