#!/usr/bin/env python3
"""Adquisición de datos para validar snowmelt-rs en la cuenca alta del Maipo.

Genera en data/:
  - dem.asc                DEM Copernicus GLO-30 (Planetary Computer),
                           UTM 19S, 200 m, recortado al box.
  - forcing.csv            date,temp_c,precip_mm — temperatura diaria ERA5
                           (Open-Meteo archive, punto central del box) +
                           precipitación diaria CR2MET v2.5 (media del box).
  - modis_YYYY-MM-DD.asc   Cobertura nival binaria MOD10A1 v6.1
                           (NDSI_Snow_Cover >= 40 → 1, < 40 → 0,
                           nubes/fill → NODATA), misma grilla del DEM.

Uso: python3 fetch_data.py [dem|forcing|modis|all]
"""

import json
import os
import subprocess
import sys
import urllib.parse
import urllib.request

import numpy as np

BOX = (-70.35, -33.85, -69.95, -33.50)  # lon_min, lat_min, lon_max, lat_max
CENTER = ((BOX[0] + BOX[2]) / 2, (BOX[1] + BOX[3]) / 2)
EPSG = "EPSG:32719"
RES = "200"
START, END = "2019-04-01", "2019-12-31"
# Fechas con baja nubosidad sobre el box (sondeadas; junio 2019 quedó fuera
# por nubosidad persistente — el mejor día tuvo 57% de celdas válidas).
MODIS_DATES = [
    "2019-07-15",
    "2019-08-08",
    "2019-09-08",
    "2019-10-15",
    "2019-11-10",
]
PC_API = "https://planetarycomputer.microsoft.com/api"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
CR2MET_DIR = os.path.expanduser(
    "~/proyectos/Agentes/CR2MET_pr_v2.5_day_1960-2021_005deg/pr"
)


def http_json(url, payload=None):
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode() if payload else None,
        headers={"Content-Type": "application/json", "User-Agent": "snowmelt-rs"},
    )
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.load(r)


def sign(href):
    url = f"{PC_API}/sas/v1/sign?href=" + urllib.parse.quote(href, safe="")
    return http_json(url)["href"]


def stac_search(collection, bbox, datetime=None):
    payload = {"collections": [collection], "bbox": list(bbox), "limit": 20}
    if datetime:
        payload["datetime"] = datetime
    return http_json(f"{PC_API}/stac/v1/search", payload)["features"]


def warp_to_grid(srcs, dst_tif, resample):
    cmd = [
        "gdalwarp", "-q", "-overwrite",
        "-t_srs", EPSG, "-tr", RES, RES, "-tap",
        "-te_srs", "EPSG:4326",
        "-te", str(BOX[0]), str(BOX[1]), str(BOX[2]), str(BOX[3]),
        "-r", resample,
        *srcs, dst_tif,
    ]
    subprocess.run(cmd, check=True)


def tif_to_array(path):
    import rasterio

    with rasterio.open(path) as src:
        data = src.read(1).astype(float)
        if src.nodata is not None:
            data[data == src.nodata] = np.nan
        t = src.transform
        return data, (t.c, t.f, t.a)  # x_origin (left), y_origin (top), cellsize


def write_asc(path, data, x_left, y_top, cellsize, nodata=-9999.0):
    rows, cols = data.shape
    yll = y_top - rows * cellsize
    with open(path, "w") as f:
        f.write(f"ncols {cols}\nnrows {rows}\n")
        f.write(f"xllcorner {x_left}\nyllcorner {yll}\ncellsize {cellsize}\n")
        f.write(f"NODATA_value {nodata}\n")
        for row in data:
            vals = [(f"{v:.3f}" if np.isfinite(v) else str(nodata)) for v in row]
            f.write(" ".join(vals) + "\n")


def fetch_dem():
    print("DEM Copernicus GLO-30...")
    items = stac_search("cop-dem-glo-30", BOX)
    hrefs = ["/vsicurl/" + sign(it["assets"]["data"]["href"]) for it in items]
    print(f"  {len(hrefs)} tiles")
    tif = os.path.join(DATA, "dem.tif")
    warp_to_grid(hrefs, tif, "bilinear")
    data, (x, y, cs) = tif_to_array(tif)
    write_asc(os.path.join(DATA, "dem.asc"), data, x, y, cs)
    valid = np.isfinite(data)
    print(
        f"  dem.asc: {data.shape}, elevación {np.nanmin(data):.0f}-"
        f"{np.nanmax(data):.0f} m, media {np.nanmean(data):.0f} m"
    )
    return float(np.nanmean(data[valid]))


def fetch_forcing():
    print("Temperatura ERA5 (Open-Meteo)...")
    url = (
        "https://archive-api.open-meteo.com/v1/era5"
        f"?latitude={CENTER[1]}&longitude={CENTER[0]}"
        f"&start_date={START}&end_date={END}"
        "&daily=temperature_2m_mean&timezone=UTC"
    )
    j = http_json(url)
    z_ref = j["elevation"]
    dates = j["daily"]["time"]
    temps = j["daily"]["temperature_2m_mean"]
    print(f"  {len(dates)} días, elevación de la celda ERA5: {z_ref} m")

    print("Precipitación CR2MET (media del box)...")
    import xarray as xr

    months = sorted(
        {d[:7] for d in dates}
    )  # YYYY-MM
    precip = {}
    for ym in months:
        y, m = ym.split("-")
        path = os.path.join(CR2MET_DIR, f"CR2MET_pr_v2.5_day_{y}_{m}_005deg.nc")
        ds = xr.open_dataset(path)
        var = "pr" if "pr" in ds else list(ds.data_vars)[0]
        da = ds[var]
        lat_name = "lat" if "lat" in da.dims else "latitude"
        lon_name = "lon" if "lon" in da.dims else "longitude"
        sub = da.sel(
            {lat_name: slice(BOX[1], BOX[3]), lon_name: slice(BOX[0], BOX[2])}
        )
        if sub[lat_name].size == 0:  # latitud descendente
            sub = da.sel(
                {lat_name: slice(BOX[3], BOX[1]), lon_name: slice(BOX[0], BOX[2])}
            )
        daily = sub.mean(dim=[lat_name, lon_name]).values
        tname = "time" if "time" in da.dims else da.dims[0]
        for t, v in zip(ds[tname].values, daily):
            precip[str(t)[:10]] = float(v)
        ds.close()

    out = os.path.join(DATA, "forcing.csv")
    with open(out, "w") as f:
        f.write("date,temp_c,precip_mm\n")
        for d, t in zip(dates, temps):
            p = precip.get(d)
            if p is None or t is None:
                print(f"  WARN: sin dato para {d}, omitido")
                continue
            f.write(f"{d},{t:.2f},{max(p, 0.0):.2f}\n")
    print(f"  forcing.csv listo (z_ref sugerido: {z_ref} m)")
    return z_ref


def fetch_modis():
    print("MODIS MOD10A1 v6.1 (NDSI snow cover)...")
    for d in MODIS_DATES:
        items = stac_search(
            "modis-10A1-061", BOX, datetime=f"{d}T00:00:00Z/{d}T23:59:59Z"
        )
        if not items:
            print(f"  {d}: sin escenas, omitido")
            continue
        hrefs = [
            "/vsicurl/" + sign(it["assets"]["NDSI_Snow_Cover"]["href"])
            for it in items
        ]
        tif = os.path.join(DATA, f"modis_{d}.tif")
        warp_to_grid(hrefs, tif, "near")
        data, (x, y, cs) = tif_to_array(tif)
        # 0-100 = NDSI; >100 = nubes/fill/etc.
        snow = np.where(
            (data >= 0) & (data <= 100), (data >= 40).astype(float), np.nan
        )
        valid_frac = np.isfinite(snow).mean()
        snow_frac = np.nanmean(snow) if valid_frac > 0 else float("nan")
        write_asc(os.path.join(DATA, f"modis_{d}.asc"), snow, x, y, cs)
        print(
            f"  {d}: {len(items)} escenas, válido {valid_frac:.0%}, "
            f"nieve {snow_frac:.0%}"
        )


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("dem", "all"):
        fetch_dem()
    if what in ("forcing", "all"):
        fetch_forcing()
    if what in ("modis", "all"):
        fetch_modis()
