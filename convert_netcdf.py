#!/usr/bin/env python3
"""
Convert translated GRIB1 files (grb2d / grbsanl) to NetCDF.

Reads GRIB1 records using eccodes, builds xarray datasets with proper
variable names, coordinates, and dimensions, then writes compressed NetCDF4.

Usage:
    python convert_netcdf.py grb2d 2026 1              # Full month
    python convert_netcdf.py grbsanl 2026 1             # Full month
    python convert_netcdf.py grb2d 2026 1 --day 1 --hour 0  # Single timestep
    python convert_netcdf.py grb2d 2026 1 -o output.nc  # Custom output path
"""

import argparse
import os
import struct
import sys
import logging
from datetime import datetime
from pathlib import Path

import numpy as np

try:
    import eccodes
    import xarray as xr
except ImportError:
    sys.exit("Required: pip install eccodes xarray netcdf4")

from config import (
    GRB2D_VARS, GRBSANL_SIGMA_VARS, GRBSANL_DERIVED_VARS, GRBSANL_SURFACE_VARS,
    SIGMA_LEVELS, SIGMA_KPDS7, GAUSS_LATS_94, NLON, NLAT, DLON,
    WORK_DIR, sixhourly_timestamps,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)


# ── GRIB1 parameter ID → NetCDF variable name ───────────────────────────
# Key: (kpds5, kpds6, kpds7) → (nc_name, long_name, units)

GRB2D_NAMES = {}
for var in GRB2D_VARS:
    key = (var["kpds5"], var["kpds6"], var["kpds7"])
    name = var["core_name"]
    # Disambiguate by level type
    if var["kpds6"] == 8:
        nc_name = f"{name}_toa"
    elif var["kpds6"] == 105 and var["kpds7"] == 2:
        nc_name = f"{name}_2m"
    elif var["kpds6"] == 105 and var["kpds7"] == 10:
        nc_name = f"{name}_10m"
    elif var["kpds6"] == 112:
        top = var["kpds7"] >> 8
        bot = var["kpds7"] & 0xFF
        nc_name = f"{name}_{top}_{bot}cm"
    elif var["kpds6"] == 111:
        nc_name = f"{name}_{var['kpds7']}cm"
    elif var["kpds6"] in (200,):
        nc_name = f"{name}_col"
    elif var["kpds6"] in (214, 224, 234):
        layer = {214: "low", 224: "mid", 234: "high"}[var["kpds6"]]
        nc_name = f"{name}_{layer}"
    else:
        nc_name = name
    # Handle TSOIL→TMP_soil, SSRUN→RUNOF, SBSNO→SNOEV, CPOFP→SRWEQ
    renames = {"TSOIL": "TMP", "SSRUN": "RUNOF", "SBSNO": "SNOEV", "CPOFP": "SRWEQ"}
    if name in renames:
        nc_name = nc_name.replace(name, renames[name])
    GRB2D_NAMES[key] = nc_name.lower()

# Special case: TMAX/TMIN use same kpds6/kpds7 as TMP_2m but different kpds5
GRB2D_NAMES[(15, 105, 2)] = "tmax_2m"
GRB2D_NAMES[(16, 105, 2)] = "tmin_2m"

GRBSANL_NAMES = {}
# Sigma-level vars
for varname, kpds5 in GRBSANL_SIGMA_VARS:
    GRBSANL_NAMES[kpds5] = varname.lower()
for varname, kpds5 in GRBSANL_DERIVED_VARS.items():
    GRBSANL_NAMES[kpds5] = varname.lower()
# Surface vars
for varname, (kpds5, kpds6, kpds7) in GRBSANL_SURFACE_VARS.items():
    GRBSANL_NAMES[(kpds5, kpds6)] = varname.lower()


# ── Coordinate arrays ────────────────────────────────────────────────────

LATS = np.array(GAUSS_LATS_94)
LONS = np.arange(0, 360, DLON)


# ── GRIB1 readers ────────────────────────────────────────────────────────

def read_grb2d_file(filepath):
    """Read a grb2d GRIB1 file, return dict of {nc_name: 2D array}."""
    fields = {}
    with open(filepath, "rb") as f:
        while True:
            msgid = eccodes.codes_grib_new_from_file(f)
            if msgid is None:
                break
            kpds5 = eccodes.codes_get_long(msgid, "indicatorOfParameter")
            kpds6 = eccodes.codes_get_long(msgid, "indicatorOfTypeOfLevel")
            kpds7 = eccodes.codes_get_long(msgid, "level")
            if kpds6 in (112, 116):
                top = eccodes.codes_get_long(msgid, "topLevel")
                bot = eccodes.codes_get_long(msgid, "bottomLevel")
                kpds7 = top * 256 + bot

            key = (kpds5, kpds6, kpds7)
            nc_name = GRB2D_NAMES.get(key, f"var{kpds5}_{kpds6}_{kpds7}")
            values = eccodes.codes_get_values(msgid)
            eccodes.codes_release(msgid)

            fields[nc_name] = values.reshape(NLAT, NLON)
    return fields


def read_grbsanl_file(filepath):
    """
    Read a grbsanl prefixed GRIB1 file.
    Returns dict of sigma-level fields {name: (nsigma, nlat, nlon)}
    and surface fields {name: (nlat, nlon)}.
    """
    sigma_fields = {}  # name → list of (sigma_idx, 2D array)
    surface_fields = {}

    with open(filepath, "rb") as f:
        data = f.read()

    pos = 0
    while pos < len(data):
        if pos + 4 > len(data):
            break
        reclen = struct.unpack(">I", data[pos:pos + 4])[0]
        pos += 4
        if pos + reclen > len(data):
            break
        msg_data = data[pos:pos + reclen]
        pos += reclen
        # Skip Fortran trailer if present
        if pos + 4 <= len(data):
            trailer = struct.unpack(">I", data[pos:pos + 4])[0]
            if trailer == reclen:
                pos += 4

        if msg_data[:4] != b"GRIB":
            continue

        try:
            msgid = eccodes.codes_new_from_message(msg_data)
        except Exception:
            continue

        kpds5 = eccodes.codes_get_long(msgid, "indicatorOfParameter")
        kpds6 = eccodes.codes_get_long(msgid, "indicatorOfTypeOfLevel")
        kpds7 = eccodes.codes_get_long(msgid, "level")
        values = eccodes.codes_get_values(msgid).reshape(NLAT, NLON)
        eccodes.codes_release(msgid)

        if kpds6 == 107:  # sigma level
            nc_name = GRBSANL_NAMES.get(kpds5, f"var{kpds5}")
            if nc_name not in sigma_fields:
                sigma_fields[nc_name] = {}
            try:
                sigma_idx = SIGMA_KPDS7.index(kpds7)
            except ValueError:
                sigma_idx = kpds7
            sigma_fields[nc_name][sigma_idx] = values
        else:  # surface
            nc_name = GRBSANL_NAMES.get((kpds5, kpds6), f"var{kpds5}_sfc")
            surface_fields[nc_name] = values

    # Assemble sigma fields into 3D arrays
    nsigma = len(SIGMA_LEVELS)
    sigma_3d = {}
    for name, level_dict in sigma_fields.items():
        arr = np.full((nsigma, NLAT, NLON), np.nan)
        for idx, vals in level_dict.items():
            if isinstance(idx, int) and idx < nsigma:
                arr[idx] = vals
        sigma_3d[name] = arr

    return sigma_3d, surface_fields


# ── Dataset builders ─────────────────────────────────────────────────────

def build_grb2d_dataset(year, month, timestamps=None, indir=None):
    """Build xarray Dataset for grb2d from translated files."""
    if indir is None:
        indir = WORK_DIR / "grb2d" / str(year) / f"{month:02d}"
    else:
        indir = Path(indir)

    if timestamps is None:
        timestamps = sixhourly_timestamps(year, month)

    times = []
    all_fields = []

    for y, m, d, h in timestamps:
        yy = y % 100
        fname = f"grb2d{yy:02d}{m:02d}{d:02d}{h:02d}"
        fpath = indir / fname
        if not fpath.exists():
            log.warning(f"Missing: {fpath}")
            continue

        fields = read_grb2d_file(str(fpath))
        all_fields.append(fields)
        times.append(datetime(y, m, d, h))

    if not all_fields:
        log.error("No files found")
        return None

    # Collect all variable names across timesteps
    all_varnames = set()
    for fields in all_fields:
        all_varnames.update(fields.keys())

    # Build data arrays
    time_coord = np.array(times, dtype="datetime64[ns]")
    data_vars = {}
    for varname in sorted(all_varnames):
        arr = np.full((len(times), NLAT, NLON), np.nan, dtype=np.float32)
        for t, fields in enumerate(all_fields):
            if varname in fields:
                arr[t] = fields[varname]
        data_vars[varname] = (["time", "latitude", "longitude"], arr)

    ds = xr.Dataset(
        data_vars,
        coords={
            "time": time_coord,
            "latitude": LATS,
            "longitude": LONS,
        },
    )
    ds.attrs["title"] = f"CDAS grb2d translated from CORe, {year}-{month:02d}"
    ds.attrs["source"] = "CORe reanalysis via CDAS-CORe translation pipeline"
    ds.attrs["Conventions"] = "CF-1.8"

    return ds


def build_grbsanl_dataset(year, month, timestamps=None, indir=None):
    """Build xarray Dataset for grbsanl from translated files."""
    if indir is None:
        indir = WORK_DIR / "grbsanl" / str(year) / f"{month:02d}"
    else:
        indir = Path(indir)

    if timestamps is None:
        timestamps = sixhourly_timestamps(year, month)

    times = []
    all_sigma = []
    all_surface = []

    for y, m, d, h in timestamps:
        yy = y % 100
        fname = f"grbsanl{yy:02d}{m:02d}{d:02d}{h:02d}"
        fpath = indir / fname
        if not fpath.exists():
            log.warning(f"Missing: {fpath}")
            continue

        sigma_3d, surface = read_grbsanl_file(str(fpath))
        all_sigma.append(sigma_3d)
        all_surface.append(surface)
        times.append(datetime(y, m, d, h))

    if not all_sigma:
        log.error("No files found")
        return None

    time_coord = np.array(times, dtype="datetime64[ns]")
    sigma_coord = np.array(SIGMA_LEVELS, dtype=np.float64)
    nt = len(times)
    nsigma = len(SIGMA_LEVELS)

    data_vars = {}

    # Sigma-level variables: (time, sigma, lat, lon)
    sigma_varnames = set()
    for fields in all_sigma:
        sigma_varnames.update(fields.keys())
    for varname in sorted(sigma_varnames):
        arr = np.full((nt, nsigma, NLAT, NLON), np.nan, dtype=np.float32)
        for t, fields in enumerate(all_sigma):
            if varname in fields:
                arr[t] = fields[varname]
        data_vars[varname] = (["time", "sigma", "latitude", "longitude"], arr)

    # Surface variables: (time, lat, lon)
    sfc_varnames = set()
    for fields in all_surface:
        sfc_varnames.update(fields.keys())
    for varname in sorted(sfc_varnames):
        arr = np.full((nt, NLAT, NLON), np.nan, dtype=np.float32)
        for t, fields in enumerate(all_surface):
            if varname in fields:
                arr[t] = fields[varname]
        data_vars[varname] = (["time", "latitude", "longitude"], arr)

    ds = xr.Dataset(
        data_vars,
        coords={
            "time": time_coord,
            "sigma": sigma_coord,
            "latitude": LATS,
            "longitude": LONS,
        },
    )
    ds.sigma.attrs["long_name"] = "sigma level"
    ds.sigma.attrs["positive"] = "down"
    ds.sigma.attrs["units"] = "1"
    ds.attrs["title"] = f"CDAS grbsanl translated from CORe, {year}-{month:02d}"
    ds.attrs["source"] = "CORe reanalysis via CDAS-CORe translation pipeline"
    ds.attrs["Conventions"] = "CF-1.8"

    return ds


# ── Write NetCDF ─────────────────────────────────────────────────────────

ENCODING_DEFAULTS = dict(
    zlib=True,
    complevel=4,
    dtype="float32",
)


def write_netcdf(ds, outpath):
    """Write dataset to compressed NetCDF4."""
    encoding = {}
    for var in ds.data_vars:
        encoding[var] = ENCODING_DEFAULTS.copy()

    ds.to_netcdf(outpath, encoding=encoding, format="NETCDF4")
    size_mb = os.path.getsize(outpath) / 1e6
    log.info(f"Wrote {outpath} ({size_mb:.1f} MB)")


# ── CLI ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Convert translated GRIB1 to NetCDF"
    )
    parser.add_argument("filetype", choices=["grb2d", "grbsanl"])
    parser.add_argument("year", type=int)
    parser.add_argument("month", type=int)
    parser.add_argument("--day", type=int, help="Single day")
    parser.add_argument("--hour", type=int, help="Single hour (requires --day)")
    parser.add_argument("-o", "--output", help="Output NetCDF path")
    parser.add_argument("--indir", help="Input directory override")
    args = parser.parse_args()

    year, month = args.year, args.month

    if args.day is not None and args.hour is not None:
        timestamps = [(year, month, args.day, args.hour)]
        default_out = f"{args.filetype}_{year}{month:02d}{args.day:02d}{args.hour:02d}.nc"
    else:
        timestamps = None
        default_out = f"{args.filetype}_{year}{month:02d}.nc"

    outpath = args.output or str(WORK_DIR / default_out)
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)

    if args.filetype == "grb2d":
        ds = build_grb2d_dataset(year, month, timestamps, args.indir)
    else:
        ds = build_grbsanl_dataset(year, month, timestamps, args.indir)

    if ds is None:
        sys.exit(1)

    log.info(f"Dataset: {dict(ds.dims)}, {len(ds.data_vars)} variables")
    write_netcdf(ds, outpath)


if __name__ == "__main__":
    main()
