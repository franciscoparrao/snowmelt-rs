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

  - grids/precip_YYYY-MM-DD.asc   Precipitación diaria CR2MET v2.5
                           regrillada a la malla del DEM (mm), para el
                           forzante distribuido (--precip-grids).
  - grids/temp_YYYY-MM-DD.asc     Temperatura diaria ERA5 multi-celda
                           (Open-Meteo) con downscaling topográfico y
                           lapse rate empírico, para --temp-grids.

Uso: python3 fetch_data.py [dem|forcing|modis|grids|tempgrids|all]
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


# Gradiente orográfico fraccional de precipitación para el downscaling de
# subgrilla [1/m]. Realza CR2MET por la anomalía de elevación que su malla
# gruesa (0.05°) no resuelve. Típico andino 0.0005–0.001.
OROG_GAMMA = 0.0008


def coarse_elevation():
    """Elevación 'vista' por la malla CR2MET (0.05°): el DEM promediado a
    esa resolución y vuelto a la malla fina. La diferencia DEM − coarse es
    la anomalía de subgrilla que CR2MET no resuelve."""
    import rasterio
    from rasterio.warp import reproject, Resampling

    dem_tif = os.path.join(DATA, "dem.tif")
    coarse_tif = os.path.join(DATA, "dem_coarse005.tif")
    fine_back = os.path.join(DATA, "dem_coarse005_fine.tif")
    # DEM → 0.05° en EPSG:4326 (promedio), luego de vuelta a la malla del DEM.
    subprocess.run(
        ["gdalwarp", "-q", "-overwrite", "-t_srs", "EPSG:4326",
         "-tr", "0.05", "0.05", "-r", "average", dem_tif, coarse_tif],
        check=True,
    )
    with rasterio.open(dem_tif) as ref:
        dst = np.full((ref.height, ref.width), np.nan)
        with rasterio.open(coarse_tif) as src:
            reproject(
                source=src.read(1), destination=dst,
                src_transform=src.transform, src_crs=src.crs,
                dst_transform=ref.transform, dst_crs=ref.crs,
                resampling=Resampling.bilinear,
            )
    for f in (coarse_tif, fine_back):
        if os.path.exists(f):
            os.remove(f)
    return dst


def fetch_precip_grids():
    """Regrilla CR2MET diario a la malla del DEM con downscaling orográfico
    de subgrilla → grids/precip_DATE.asc.

    CR2MET (0.05°) se interpola bilinealmente y luego se realza por la
    anomalía de elevación de subgrilla: P = P_coarse·(1 + γ·(z − z_coarse)),
    de modo que las cumbres reciben más y los valles menos dentro de cada
    celda CR2MET (la media gruesa se preserva). Corrige el gradiente
    orográfico que la malla de 5 km no resuelve.
    """
    print("Grillas de precipitación CR2MET + downscaling orográfico → DEM...")
    import xarray as xr
    import rasterio
    from rasterio.warp import reproject, Resampling
    from rasterio.transform import from_origin

    dem_fine, _ = tif_to_array(os.path.join(DATA, "dem.tif"))
    z_coarse = coarse_elevation()
    # Factor de realce orográfico (>= 0), 1 donde no hay anomalía.
    orog = np.where(
        np.isfinite(dem_fine) & np.isfinite(z_coarse),
        np.maximum(1.0 + OROG_GAMMA * (dem_fine - z_coarse), 0.0),
        1.0,
    )

    # Geometría destino desde el DEM (.asc).
    _, (x_left, y_top, cs) = tif_to_array(os.path.join(DATA, "dem.tif"))
    with rasterio.open(os.path.join(DATA, "dem.tif")) as ref:
        dst_shape = (ref.height, ref.width)
        dst_transform = ref.transform
        dst_crs = ref.crs

    grids_dir = os.path.join(DATA, "grids")
    os.makedirs(grids_dir, exist_ok=True)

    # Fechas de la serie de forzantes.
    with open(os.path.join(DATA, "forcing.csv")) as f:
        next(f)
        dates = [line.split(",")[0] for line in f if line.strip()]
    months = sorted({d[:7] for d in dates})

    src_crs = "EPSG:4326"
    n = 0
    for ym in months:
        y, m = ym.split("-")
        path = os.path.join(CR2MET_DIR, f"CR2MET_pr_v2.5_day_{y}_{m}_005deg.nc")
        ds = xr.open_dataset(path)
        var = "pr" if "pr" in ds else list(ds.data_vars)[0]
        da = ds[var]
        lat_name = "lat" if "lat" in da.dims else "latitude"
        lon_name = "lon" if "lon" in da.dims else "longitude"
        lats = ds[lat_name].values
        lons = ds[lon_name].values
        res_lat = abs(float(lats[1] - lats[0]))
        res_lon = abs(float(lons[1] - lons[0]))
        # Transform de la grilla CR2MET (origen esquina superior izquierda).
        src_transform = from_origin(
            float(lons.min()) - res_lon / 2,
            float(lats.max()) + res_lat / 2,
            res_lon,
            res_lat,
        )
        flip = lats[0] < lats[-1]  # lat ascendente → voltear a N-arriba
        tname = "time" if "time" in da.dims else da.dims[0]
        day_strs = [str(t)[:10] for t in ds[tname].values]
        for i, d in enumerate(day_strs):
            if d not in dates:
                continue
            src = da.isel({tname: i}).values.astype("float64")
            if flip:
                src = src[::-1, :]
            dst = np.full(dst_shape, np.nan)
            reproject(
                source=src,
                destination=dst,
                src_transform=src_transform,
                src_crs=src_crs,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                resampling=Resampling.bilinear,
                src_nodata=np.nan,
                dst_nodata=np.nan,
            )
            dst = np.where(np.isfinite(dst), np.maximum(dst * orog, 0.0), np.nan)
            write_asc(
                os.path.join(grids_dir, f"precip_{d}.asc"), dst, x_left, y_top, cs
            )
            n += 1
        ds.close()
    print(f"  {n} grillas precip_*.asc en {grids_dir} (γ_orog={OROG_GAMMA})")


def fetch_temp_grids():
    """Temperatura distribuida ERA5 multi-celda → grids/temp_DATE.asc.

    Consulta una malla de puntos ERA5 (Open-Meteo), deriva un lapse rate
    EMPÍRICO diario por regresión T-vs-z de las celdas del modelo (no el
    −6.5 °C/km asumido), reduce cada punto a nivel del mar, interpola
    horizontalmente al DEM y re-extrapola con ese lapse a la elevación de
    cada celda (downscaling topográfico estándar).
    """
    print("Temperatura distribuida ERA5 multi-celda → DEM...")
    import rasterio
    from rasterio.warp import transform as warp_transform
    from scipy.interpolate import griddata

    # Malla de puntos sobre el box (6×6); Open-Meteo los snapea a celdas ERA5.
    n_side = 6
    lats = np.linspace(BOX[1], BOX[3], n_side)
    lons = np.linspace(BOX[0], BOX[2], n_side)
    grid_lats, grid_lons = np.meshgrid(lats, lons)
    lat_q = ",".join(f"{v:.4f}" for v in grid_lats.ravel())
    lon_q = ",".join(f"{v:.4f}" for v in grid_lons.ravel())
    url = (
        "https://archive-api.open-meteo.com/v1/era5"
        f"?latitude={lat_q}&longitude={lon_q}"
        f"&start_date={START}&end_date={END}"
        "&daily=temperature_2m_mean&timezone=UTC"
    )
    resp = http_json(url)
    resp = resp if isinstance(resp, list) else [resp]

    # Deduplicar celdas ERA5 (varios puntos snapean a la misma).
    cells = {}
    dates = None
    for r in resp:
        key = (round(r["latitude"], 4), round(r["longitude"], 4))
        if key in cells:
            continue
        cells[key] = (r["elevation"], np.array(r["daily"]["temperature_2m_mean"], float))
        dates = r["daily"]["time"]
    pts_lat = np.array([k[0] for k in cells])
    pts_lon = np.array([k[1] for k in cells])
    pts_z = np.array([v[0] for v in cells.values()])
    pts_t = np.array([v[1] for v in cells.values()])  # (n_cells, n_days)
    print(f"  {len(cells)} celdas ERA5 únicas, {len(dates)} días")

    # Puntos y celdas del DEM en UTM (interpolación métrica).
    px, py = warp_transform("EPSG:4326", EPSG, list(pts_lon), list(pts_lat))
    px, py = np.array(px), np.array(py)
    dem_fine, (x_left, y_top, cs) = tif_to_array(os.path.join(DATA, "dem.tif"))
    rows, cols = dem_fine.shape
    cx = x_left + (np.arange(cols) + 0.5) * cs
    cy = y_top - (np.arange(rows) + 0.5) * cs
    mesh_x, mesh_y = np.meshgrid(cx, cy)
    cell_xy = np.column_stack([mesh_x.ravel(), mesh_y.ravel()])
    src_xy = np.column_stack([px, py])

    grids_dir = os.path.join(DATA, "grids")
    os.makedirs(grids_dir, exist_ok=True)
    keep = set(read_forcing_dates())
    n = 0
    lapses = []
    for i, d in enumerate(dates):
        if d not in keep:
            continue
        t_day = pts_t[:, i]
        # Lapse empírico del día por OLS (T vs z), acotado a rango físico.
        slope = np.polyfit(pts_z, t_day, 1)[0]
        lapse = float(np.clip(slope, -0.0098, -0.004))
        lapses.append(lapse)
        t0 = t_day - lapse * pts_z  # reducción a nivel del mar
        lin = griddata(src_xy, t0, cell_xy, method="linear")
        near = griddata(src_xy, t0, cell_xy, method="nearest")
        t0_grid = np.where(np.isfinite(lin), lin, near).reshape(rows, cols)
        t_grid = np.where(np.isfinite(dem_fine), t0_grid + lapse * dem_fine, np.nan)
        write_asc(os.path.join(grids_dir, f"temp_{d}.asc"), t_grid, x_left, y_top, cs)
        n += 1
    print(f"  {n} grillas temp_*.asc (lapse empírico medio "
          f"{np.mean(lapses) * 1000:.2f} °C/km, rango {min(lapses) * 1000:.1f}"
          f"…{max(lapses) * 1000:.1f})")


def read_forcing_dates():
    with open(os.path.join(DATA, "forcing.csv")) as f:
        next(f)
        return [line.split(",")[0] for line in f if line.strip()]


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("dem", "all"):
        fetch_dem()
    if what in ("forcing", "all"):
        fetch_forcing()
    if what in ("modis", "all"):
        fetch_modis()
    if what in ("grids", "all"):
        fetch_precip_grids()
    if what in ("tempgrids", "all"):
        fetch_temp_grids()
