#!/usr/bin/env python3
"""
Translate CORe flx (GRIB2, 512x256, 3-hourly) → GDEX grb2d (GRIB1, 192x94, 6-hourly).

For each 6-hour window:
  - Instantaneous fields: use the synoptic-time file (00/06/12/18Z)
  - Time-averaged fields: average the two 3-hr CORe files spanning the window
  - TMAX/TMIN: take max/min across the two 3-hr periods

Pipeline per variable: wgrib2 extract → cdo remapbil → grib_set metadata fix
"""

import os
import sys
import shutil
import tempfile
import struct
import logging
from pathlib import Path
from datetime import datetime, timedelta

import numpy as np

try:
    import eccodes
except ImportError:
    sys.exit("eccodes Python package required: pip install eccodes")

from config import (
    GRB2D_VARS, WORK_DIR, TARGET_GRID,
    CENTER, SUBCENTER, PROCESS, TABLE_VERSION,
    core_flx_path, sixhourly_timestamps,
    run_cmd, wgrib2_extract, cdo_remap, grib_set,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)


def build_wgrib2_pattern(var):
    """Build a wgrib2 -match regex for a variable definition."""
    name = var["core_name"]
    level = var["core_level"]
    ctype = var["core_type"]
    if ctype == "anl":
        return f":{name}:{level}:anl:"
    elif ctype == "ave":
        return f":{name}:{level}:.*ave fcst:"
    elif ctype == "max":
        return f":{name}:{level}:.*max fcst:"
    elif ctype == "min":
        return f":{name}:{level}:.*min fcst:"
    elif ctype == "acc":
        return f":{name}:{level}:.*acc fcst:"
    else:
        raise ValueError(f"Unknown core_type: {ctype}")


def read_grib_data(filepath):
    """Read all data values from a single-record GRIB file. Returns numpy array."""
    with open(filepath, "rb") as f:
        msgid = eccodes.codes_grib_new_from_file(f)
        if msgid is None:
            raise RuntimeError(f"No GRIB message in {filepath}")
        values = eccodes.codes_get_values(msgid)
        eccodes.codes_release(msgid)
    return values


def write_grib_data(template_path, values, output_path):
    """Write data values to a GRIB file using a template for metadata."""
    with open(template_path, "rb") as f:
        msgid = eccodes.codes_grib_new_from_file(f)
    eccodes.codes_set_values(msgid, values)
    with open(output_path, "wb") as f:
        eccodes.codes_write(msgid, f)
    eccodes.codes_release(msgid)


def process_variable(var, flx_t0, flx_t1, tmpdir, year, month, day, hour):
    """
    Process a single variable for one 6-hour window.

    Args:
        var: variable definition dict from GRB2D_VARS
        flx_t0: path to CORe flx file at window start (e.g., 00Z)
        flx_t1: path to CORe flx file at window start + 3h (e.g., 03Z)
        tmpdir: temporary directory for intermediate files
        year, month, day, hour: synoptic time for this window

    Returns:
        Path to final GRIB1 file with correct metadata, or None on failure.
    """
    pattern = build_wgrib2_pattern(var)
    ctype = var["core_type"]
    label = f"{var['core_name']}_{var['core_level']}_{ctype}".replace(" ", "_")

    # Step 1: Extract from CORe GRIB2
    if ctype == "anl":
        # Instantaneous: use only the synoptic-time file
        extract_path = os.path.join(tmpdir, f"{label}_extract.grb2")
        try:
            wgrib2_extract(flx_t0, pattern, extract_path)
        except RuntimeError as e:
            log.warning(f"Extract failed for {label}: {e}")
            return None
        if not os.path.exists(extract_path) or os.path.getsize(extract_path) == 0:
            log.warning(f"No data extracted for {label}")
            return None
        remap_input = extract_path

    elif ctype in ("ave", "acc"):
        # Time-averaged/accumulated: average the two 3-hr files
        ext0 = os.path.join(tmpdir, f"{label}_ext0.grb2")
        ext1 = os.path.join(tmpdir, f"{label}_ext1.grb2")
        try:
            wgrib2_extract(flx_t0, pattern, ext0)
            wgrib2_extract(flx_t1, pattern, ext1)
        except RuntimeError as e:
            log.warning(f"Extract failed for {label}: {e}")
            return None
        if not all(os.path.exists(p) and os.path.getsize(p) > 0 for p in [ext0, ext1]):
            log.warning(f"Missing data for averaging {label}")
            return None
        # Average via cdo
        avg_path = os.path.join(tmpdir, f"{label}_avg.grb2")
        try:
            run_cmd(f"cdo ensmean {ext0} {ext1} {avg_path}")
        except RuntimeError as e:
            log.warning(f"cdo ensmean failed for {label}: {e}")
            return None
        remap_input = avg_path

    elif ctype == "max":
        # Take max across two 3-hr periods
        ext0 = os.path.join(tmpdir, f"{label}_ext0.grb2")
        ext1 = os.path.join(tmpdir, f"{label}_ext1.grb2")
        try:
            wgrib2_extract(flx_t0, pattern, ext0)
            wgrib2_extract(flx_t1, pattern, ext1)
        except RuntimeError as e:
            log.warning(f"Extract failed for {label}: {e}")
            return None
        if not all(os.path.exists(p) and os.path.getsize(p) > 0 for p in [ext0, ext1]):
            log.warning(f"Missing data for max {label}")
            return None
        max_path = os.path.join(tmpdir, f"{label}_max.grb2")
        try:
            run_cmd(f"cdo ensmax {ext0} {ext1} {max_path}")
        except RuntimeError as e:
            log.warning(f"cdo ensmax failed for {label}: {e}")
            return None
        remap_input = max_path

    elif ctype == "min":
        # Take min across two 3-hr periods
        ext0 = os.path.join(tmpdir, f"{label}_ext0.grb2")
        ext1 = os.path.join(tmpdir, f"{label}_ext1.grb2")
        try:
            wgrib2_extract(flx_t0, pattern, ext0)
            wgrib2_extract(flx_t1, pattern, ext1)
        except RuntimeError as e:
            log.warning(f"Extract failed for {label}: {e}")
            return None
        if not all(os.path.exists(p) and os.path.getsize(p) > 0 for p in [ext0, ext1]):
            log.warning(f"Missing data for min {label}")
            return None
        min_path = os.path.join(tmpdir, f"{label}_min.grb2")
        try:
            run_cmd(f"cdo ensmin {ext0} {ext1} {min_path}")
        except RuntimeError as e:
            log.warning(f"cdo ensmin failed for {label}: {e}")
            return None
        remap_input = min_path
    else:
        log.error(f"Unknown core_type: {ctype}")
        return None

    # Step 2: Remap to 192x94 Gaussian + convert to GRIB1
    remap_path = os.path.join(tmpdir, f"{label}_remap.grb1")
    try:
        cdo_remap(remap_input, remap_path)
    except RuntimeError as e:
        log.warning(f"Remap failed for {label}: {e}")
        return None

    # Step 3: Fix GRIB1 metadata
    tr = 10 if var["time_type"] == "inst" else 3
    p1 = 0
    p2 = 6
    final_path = os.path.join(tmpdir, f"{label}_final.grb1")

    # For layer-type levels (kpds6=112), kpds7 encodes top*256+bottom
    # which exceeds 255 and must be set via topLevel/bottomLevel separately
    kpds6 = var["kpds6"]
    kpds7 = var["kpds7"]
    level_keys = {}
    if kpds6 in (112, 116):
        # Layer types encode topLevel*256+bottomLevel in kpds7
        level_keys["topLevel"] = kpds7 >> 8
        level_keys["bottomLevel"] = kpds7 & 0xFF
    else:
        level_keys["level"] = kpds7

    try:
        grib_set(
            remap_path, final_path,
            indicatorOfParameter=var["kpds5"],
            indicatorOfTypeOfLevel=kpds6,
            **level_keys,
            table2Version=TABLE_VERSION,
            centre=CENTER,
            subCentre=SUBCENTER,
            generatingProcessIdentifier=PROCESS,
            timeRangeIndicator=tr,
            P1=p1,
            P2=p2,
            unitOfTimeRange=1,
            numberIncludedInAverage=0,
            yearOfCentury=year % 100,
            month=month,
            day=day,
            hour=hour,
            minute=0,
        )
    except RuntimeError as e:
        log.warning(f"grib_set failed for {label}: {e}")
        return None

    return final_path


def translate_timestep(year, month, day, hour, outdir, tmpdir=None):
    """
    Translate one 6-hourly timestep of CORe flx to GDEX grb2d format.

    Args:
        year, month, day, hour: the 6-hourly synoptic time
        outdir: directory to write output file
        tmpdir: temp directory (created if None)

    Returns:
        Path to output grb2d file, or None on failure.
    """
    # The two CORe 3-hr files spanning this 6-hr window
    t0 = datetime(year, month, day, hour)
    t1 = t0 + timedelta(hours=3)

    flx_t0 = core_flx_path(t0.year, t0.month, t0.day, t0.hour)
    flx_t1 = core_flx_path(t1.year, t1.month, t1.day, t1.hour)

    if not flx_t0.exists():
        log.error(f"Missing CORe flx file: {flx_t0}")
        return None
    if not flx_t1.exists():
        log.error(f"Missing CORe flx file: {flx_t1}")
        return None

    yy = year % 100
    outname = f"grb2d{yy:02d}{month:02d}{day:02d}{hour:02d}"
    outpath = os.path.join(outdir, outname)

    cleanup = tmpdir is None
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix="grb2d_")

    try:
        record_files = []
        for var in GRB2D_VARS:
            result = process_variable(var, flx_t0, flx_t1, tmpdir, year, month, day, hour)
            if result is not None:
                record_files.append(result)

        if not record_files:
            log.error(f"No records produced for {outname}")
            return None

        # Concatenate all GRIB1 records into one file
        with open(outpath, "wb") as outf:
            for rf in record_files:
                with open(rf, "rb") as inf:
                    outf.write(inf.read())

        log.info(f"Wrote {outname} with {len(record_files)} records")
        return outpath

    finally:
        if cleanup and os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)


def translate_month(year, month, outdir=None):
    """Translate all 6-hourly timesteps for a month."""
    if outdir is None:
        outdir = str(WORK_DIR / "grb2d" / f"{year}" / f"{month:02d}")
    os.makedirs(outdir, exist_ok=True)

    timestamps = sixhourly_timestamps(year, month)
    results = []
    for y, m, d, h in timestamps:
        path = translate_timestep(y, m, d, h, outdir)
        if path:
            results.append(path)
        else:
            log.warning(f"Failed timestep: {y}-{m:02d}-{d:02d} {h:02d}Z")

    log.info(f"Completed {len(results)}/{len(timestamps)} timesteps for {year}-{month:02d}")
    return results


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} YYYY MM [outdir]")
        sys.exit(1)
    year = int(sys.argv[1])
    month = int(sys.argv[2])
    outdir = sys.argv[3] if len(sys.argv) > 3 else None
    translate_month(year, month, outdir)
