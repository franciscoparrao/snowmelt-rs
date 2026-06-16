#!/usr/bin/env python3
"""Valida la señal de deshielo de snowmelt-rs contra el caudal observado
CAMELS-CL (mm/día) en el Río Choapa en Cuncumén (gauge 4703002).

Honestidad metodológica: snowmelt produce aporte vertical (lluvia +
derretimiento), NO un balance lluvia-escorrentía (no resta
evapotranspiración ni añade flujo base / almacenamiento de suelo). Por eso:

  * La señal de DESHIELO (melt), ruteada por un reservorio lineal, se valida
    contra la firma estacional del caudal — es la componente que un modelo
    nival debe acertar.
  * El aporte TOTAL (rain+melt) se reporta como contraste: es peor porque la
    lluvia invernal va directa al caudal sin el balance que aporta rainflow.
    Ese contraste es, justamente, la motivación de la interfaz hacia
    rainflow (que consume este aporte y cierra el balance).

Uso: python3 validate_flow.py [series.csv] [--k 90] [--warmup-year 1980]
"""

import argparse
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load_col(path, col):
    out = {}
    for r in csv.DictReader(open(path)):
        v = r.get(col, "")
        out[r["date"]] = float(v) if v not in ("", "NA") else np.nan
    return out


def linear_reservoir(x, k):
    """Mismo reservorio lineal que snowmelt-core::routing (forma exacta)."""
    s, dec = 0.0, np.exp(-1.0 / k)
    out = np.empty_like(x)
    for i, v in enumerate(x):
        if not np.isfinite(v):
            out[i] = np.nan
            continue
        s_new = s * dec + v * k * (1.0 - dec)
        out[i] = v + s - s_new
        s = s_new
    return out


def nse(sim, obs):
    return 1.0 - np.sum((sim - obs) ** 2) / np.sum((obs - np.mean(obs)) ** 2)


def metrics(sig, obs, months):
    """corr diaria, NSE de forma (escalado al volumen) y corr del ciclo anual."""
    m = np.isfinite(sig) & np.isfinite(obs)
    s, o, mo = sig[m], obs[m], months[m]
    corr = np.corrcoef(s, o)[0, 1]
    scaled = s * (o.mean() / s.mean())
    cyc_o = np.array([o[mo == i].mean() for i in range(1, 13)])
    cyc_s = np.array([scaled[mo == i].mean() for i in range(1, 13)])
    return corr, nse(scaled, o), np.corrcoef(cyc_o, cyc_s)[0, 1], cyc_o, cyc_s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("series", nargs="?", default=os.path.join(HERE, "out", "series.csv"))
    ap.add_argument("--k", type=float, default=90.0, help="recesión del reservorio [días]")
    ap.add_argument("--warmup-year", type=int, default=1980)
    args = ap.parse_args()

    melt = load_col(args.series, "melt_mm")
    runoff = load_col(args.series, "runoff_mm")
    qobs = load_col(os.path.join(HERE, "data", "qobs.csv"), "qobs_mm")

    dates = sorted(set(melt) & set(qobs))
    dates = [d for d in dates if int(d[:4]) >= args.warmup_year]
    months = np.array([int(d[5:7]) for d in dates])
    obs = np.array([qobs[d] for d in dates])
    melt_r = linear_reservoir(np.array([melt[d] for d in dates]), args.k)
    runoff_r = linear_reservoir(np.array([runoff[d] for d in dates]), args.k)

    n = np.isfinite(obs).sum()
    print(f"Río Choapa en Cuncumén (4703002) — {dates[0]}…{dates[-1]}, {n} días con qobs")
    print(f"Reservorio lineal k = {args.k:.0f} días\n")

    cm, nm, cym, cyc_o, cyc_s = metrics(melt_r, obs, months)
    cr, nr, cyr, _, _ = metrics(runoff_r, obs, months)
    print(f"{'señal':>12}  corr_diaria  NSE_forma  corr_ciclo_anual")
    print(f"{'deshielo':>12}     {cm:.3f}      {nm:+.3f}      {cym:.3f}")
    print(f"{'lluvia+desh.':>12}     {cr:.3f}      {nr:+.3f}      {cyr:.3f}")

    print("\n— Ciclo anual medio del caudal (mm/día) —")
    print("  mes  qobs  deshielo_esc")
    for i in range(12):
        bar = "█" * int(cyc_o[i] * 10)
        print(f"  {i+1:>2}  {cyc_o[i]:4.2f}  {cyc_s[i]:4.2f}  {bar}")
    print(f"  peak observado=mes {int(np.argmax(cyc_o))+1}, "
          f"deshielo=mes {int(np.argmax(cyc_s))+1}")
    print("\nEl deshielo reproduce la firma estacional del caudal; el aporte total")
    print("requiere el balance lluvia-escorrentía de rainflow (ver README).")


if __name__ == "__main__":
    main()
