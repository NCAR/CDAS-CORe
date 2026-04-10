#!/usr/bin/env python3
"""
Translate CORe pgb (GRIB2, 512x256, 3-hourly pressure levels) →
GDEX grbsanl (GRIB1, 192x94, 6-hourly sigma levels).

Pipeline per 6-hourly timestep:
1. Extract UGRD, VGRD, TMP, SPFH, HGT, PRES from CORe pgb on pressure levels
2. Remap spatially from 512x256 → 192x94 Gaussian (bilinear via cdo)
3. Vertical interpolation: pressure levels → 28 CDAS sigma levels
   using log-pressure linear interpolation (p_sigma = sigma * Ps)
4. Compute derived fields: RELV, RELD (vorticity, divergence),
   HGTX, HGTY (terrain gradients), LPSX, LPSY (log-pressure gradients)
5. Write GRIB1 with correct metadata
6. Prepend 4-byte big-endian length prefix (Fortran sequential format)
"""

import os
import sys
import shutil
import struct
import tempfile
import logging
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

try:
    import eccodes
except ImportError:
    sys.exit("eccodes Python package required: pip install eccodes")

from config import (
    SIGMA_LEVELS, SIGMA_KPDS7, PRESSURE_LEVELS,
    GRBSANL_SIGMA_VARS, GRBSANL_DERIVED_VARS, GRBSANL_SURFACE_VARS,
    GAUSS_LATS_94, NLON, NLAT, DLON, EARTH_RADIUS,
    CENTER, SUBCENTER, PROCESS, TABLE_VERSION,
    WORK_DIR, TARGET_GRID,
    core_pgb_path, sixhourly_timestamps,
    run_cmd, wgrib2_extract, cdo_remap, grib_set,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)


def extract_and_remap_pressure_field(pgb_file, varname, level_mb, tmpdir):
    """
    Extract a single variable at a single pressure level from CORe pgb,
    remap to 192x94 Gaussian, return as numpy array (NLAT*NLON,).
    """
    label = f"{varname}_{level_mb}mb"
    pattern = f":{varname}:{level_mb} mb:anl:"
    ext_path = os.path.join(tmpdir, f"{label}_ext.grb2")
    remap_path = os.path.join(tmpdir, f"{label}_remap.grb1")

    wgrib2_extract(pgb_file, pattern, ext_path)
    if not os.path.exists(ext_path) or os.path.getsize(ext_path) == 0:
        raise RuntimeError(f"No data extracted for {label}")

    cdo_remap(ext_path, remap_path)

    # Read remapped data
    with open(remap_path, "rb") as f:
        msgid = eccodes.codes_grib_new_from_file(f)
        values = eccodes.codes_get_values(msgid)
        eccodes.codes_release(msgid)
    return values


def extract_and_remap_surface(pgb_file, varname, level_str, tmpdir):
    """Extract and remap a surface-level field."""
    label = f"{varname}_sfc"
    pattern = f":{varname}:{level_str}:anl:"
    ext_path = os.path.join(tmpdir, f"{label}_ext.grb2")
    remap_path = os.path.join(tmpdir, f"{label}_remap.grb1")

    wgrib2_extract(pgb_file, pattern, ext_path)
    if not os.path.exists(ext_path) or os.path.getsize(ext_path) == 0:
        raise RuntimeError(f"No data extracted for {label}")

    cdo_remap(ext_path, remap_path)

    with open(remap_path, "rb") as f:
        msgid = eccodes.codes_grib_new_from_file(f)
        values = eccodes.codes_get_values(msgid)
        eccodes.codes_release(msgid)
    return values


def get_grib1_template(tmpdir):
    """
    Get a GRIB1 template file path from a remapped field.
    Returns the path to a single-record GRIB1 file suitable as template.
    """
    # We'll create this during processing and cache it
    return os.path.join(tmpdir, "template.grb1")


def write_grib1_record(values, template_path, output_path,
                       kpds5, kpds6, kpds7, tr=10, p1=0, p2=0,
                       year=2026, month=1, day=1, hour=0):
    """
    Write a single GRIB1 record using eccodes, based on a template.
    """
    with open(template_path, "rb") as f:
        msgid = eccodes.codes_grib_new_from_file(f)

    eccodes.codes_set(msgid, "indicatorOfParameter", kpds5)
    eccodes.codes_set(msgid, "indicatorOfTypeOfLevel", kpds6)
    eccodes.codes_set(msgid, "level", kpds7)
    eccodes.codes_set(msgid, "table2Version", TABLE_VERSION)
    eccodes.codes_set(msgid, "centre", CENTER)
    eccodes.codes_set(msgid, "subCentre", SUBCENTER)
    eccodes.codes_set(msgid, "generatingProcessIdentifier", PROCESS)
    eccodes.codes_set(msgid, "timeRangeIndicator", tr)
    eccodes.codes_set(msgid, "P1", p1)
    eccodes.codes_set(msgid, "P2", p2)
    eccodes.codes_set(msgid, "unitOfTimeRange", 1)
    eccodes.codes_set(msgid, "yearOfCentury", year % 100)
    eccodes.codes_set(msgid, "month", month)
    eccodes.codes_set(msgid, "day", day)
    eccodes.codes_set(msgid, "hour", hour)
    eccodes.codes_set(msgid, "minute", 0)
    eccodes.codes_set(msgid, "numberIncludedInAverage", 0)

    eccodes.codes_set_values(msgid, values.astype(np.float64))

    with open(output_path, "wb") as f:
        eccodes.codes_write(msgid, f)
    eccodes.codes_release(msgid)


def vertical_interp_log_pressure(var_on_plevs, ps, pressure_levels, sigma_levels):
    """
    Interpolate a field from pressure levels to sigma levels using
    log-pressure linear interpolation.

    Args:
        var_on_plevs: (nlevels, ngridpoints) array of variable on pressure levels
        ps: (ngridpoints,) surface pressure in Pa
        pressure_levels: list of pressure levels in mb, descending
        sigma_levels: list of target sigma values

    Returns:
        (nsigma, ngridpoints) array on sigma levels
    """
    npts = ps.shape[0]
    nsigma = len(sigma_levels)
    plevs_pa = np.array(pressure_levels) * 100.0  # mb → Pa
    log_plevs = np.log(plevs_pa)

    result = np.zeros((nsigma, npts), dtype=np.float64)

    for i in range(npts):
        # Target pressures at this gridpoint
        target_p = np.array(sigma_levels) * ps[i]
        log_target = np.log(np.maximum(target_p, 1.0))  # avoid log(0)

        # Interpolate in log-pressure space
        # pressure_levels are sorted descending (1000→1 mb), so log_plevs is descending
        # interp1d needs monotonically increasing x
        f = interp1d(
            log_plevs[::-1], var_on_plevs[::-1, i],
            kind="linear", bounds_error=False, fill_value="extrapolate"
        )
        result[:, i] = f(log_target)

    return result


def compute_vorticity_divergence(u_sigma, v_sigma, lats, nlon):
    """
    Compute relative vorticity and divergence on a Gaussian grid.

    RELV = (1/(a cos φ)) ∂v/∂λ - (1/a) ∂(u cos φ)/∂φ
    RELD = (1/(a cos φ)) ∂u/∂λ + (1/(a cos φ)) ∂(v cos φ)/∂φ

    Args:
        u_sigma: (nlat, nlon) u-wind on a single sigma level
        v_sigma: (nlat, nlon) v-wind
        lats: latitude array (degrees, N→S)
        nlon: number of longitude points

    Returns:
        (relv, reld) each (nlat, nlon)
    """
    nlat = len(lats)
    a = EARTH_RADIUS
    dlon_rad = np.radians(DLON)
    lat_rad = np.radians(lats)
    cos_lat = np.cos(lat_rad)

    # Avoid division by zero at poles
    cos_lat = np.where(np.abs(cos_lat) < 1e-10, 1e-10, cos_lat)

    # ∂v/∂λ and ∂u/∂λ using centered differences (periodic in longitude)
    dvdlon = np.zeros_like(v_sigma)
    dudlon = np.zeros_like(u_sigma)
    dvdlon[:, 1:-1] = (v_sigma[:, 2:] - v_sigma[:, :-2]) / (2 * dlon_rad)
    dvdlon[:, 0] = (v_sigma[:, 1] - v_sigma[:, -1]) / (2 * dlon_rad)
    dvdlon[:, -1] = (v_sigma[:, 0] - v_sigma[:, -2]) / (2 * dlon_rad)
    dudlon[:, 1:-1] = (u_sigma[:, 2:] - u_sigma[:, :-2]) / (2 * dlon_rad)
    dudlon[:, 0] = (u_sigma[:, 1] - u_sigma[:, -1]) / (2 * dlon_rad)
    dudlon[:, -1] = (u_sigma[:, 0] - u_sigma[:, -2]) / (2 * dlon_rad)

    # ∂(u cos φ)/∂φ and ∂(v cos φ)/∂φ using centered differences
    u_cosph = u_sigma * cos_lat[:, np.newaxis]
    v_cosph = v_sigma * cos_lat[:, np.newaxis]

    # dlat varies for Gaussian grid - compute from actual latitudes
    ducosph_dphi = np.zeros_like(u_sigma)
    dvcosph_dphi = np.zeros_like(v_sigma)
    for j in range(nlat):
        if j == 0:
            dphi = np.radians(lats[0] - lats[1])
            ducosph_dphi[j, :] = (u_cosph[0, :] - u_cosph[1, :]) / dphi
            dvcosph_dphi[j, :] = (v_cosph[0, :] - v_cosph[1, :]) / dphi
        elif j == nlat - 1:
            dphi = np.radians(lats[-2] - lats[-1])
            ducosph_dphi[j, :] = (u_cosph[-2, :] - u_cosph[-1, :]) / dphi
            dvcosph_dphi[j, :] = (v_cosph[-2, :] - v_cosph[-1, :]) / dphi
        else:
            dphi = np.radians(lats[j - 1] - lats[j + 1])
            ducosph_dphi[j, :] = (u_cosph[j - 1, :] - u_cosph[j + 1, :]) / dphi
            dvcosph_dphi[j, :] = (v_cosph[j - 1, :] - v_cosph[j + 1, :]) / dphi

    relv = (dvdlon / (a * cos_lat[:, np.newaxis])
            - ducosph_dphi / a)
    reld = (dudlon / (a * cos_lat[:, np.newaxis])
            + dvcosph_dphi / (a * cos_lat[:, np.newaxis]))

    return relv, reld


def compute_gradients(field_2d, lats, nlon):
    """
    Compute x and y gradients of a 2D field on a Gaussian grid.

    ∂f/∂x = (1/(a cos φ)) ∂f/∂λ
    ∂f/∂y = (1/a) ∂f/∂φ

    Returns (dfdx, dfdy) each (nlat, nlon).
    """
    nlat = len(lats)
    a = EARTH_RADIUS
    dlon_rad = np.radians(DLON)
    lat_rad = np.radians(lats)
    cos_lat = np.cos(lat_rad)
    cos_lat = np.where(np.abs(cos_lat) < 1e-10, 1e-10, cos_lat)

    # ∂f/∂λ (periodic)
    dfdlon = np.zeros_like(field_2d)
    dfdlon[:, 1:-1] = (field_2d[:, 2:] - field_2d[:, :-2]) / (2 * dlon_rad)
    dfdlon[:, 0] = (field_2d[:, 1] - field_2d[:, -1]) / (2 * dlon_rad)
    dfdlon[:, -1] = (field_2d[:, 0] - field_2d[:, -2]) / (2 * dlon_rad)

    dfdx = dfdlon / (a * cos_lat[:, np.newaxis])

    # ∂f/∂φ
    dfdphi = np.zeros_like(field_2d)
    for j in range(nlat):
        if j == 0:
            dphi = np.radians(lats[0] - lats[1])
            dfdphi[j, :] = (field_2d[0, :] - field_2d[1, :]) / dphi
        elif j == nlat - 1:
            dphi = np.radians(lats[-2] - lats[-1])
            dfdphi[j, :] = (field_2d[-2, :] - field_2d[-1, :]) / dphi
        else:
            dphi = np.radians(lats[j - 1] - lats[j + 1])
            dfdphi[j, :] = (field_2d[j - 1, :] - field_2d[j + 1, :]) / dphi

    dfdy = dfdphi / a

    return dfdx, dfdy


def translate_timestep(year, month, day, hour, outdir, tmpdir=None):
    """
    Translate one 6-hourly timestep of CORe pgb to GDEX grbsanl format.

    Returns path to output file, or None on failure.
    """
    pgb_file = core_pgb_path(year, month, day, hour)
    if not pgb_file.exists():
        log.error(f"Missing CORe pgb file: {pgb_file}")
        return None

    yy = year % 100
    outname = f"grbsanl{yy:02d}{month:02d}{day:02d}{hour:02d}"
    outpath = os.path.join(outdir, outname)

    if os.path.exists(outpath) and os.path.getsize(outpath) > 0:
        log.info(f"Skipping {outname} (already exists)")
        return outpath

    cleanup = tmpdir is None
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix="grbsanl_")

    try:
        log.info(f"Processing {outname}...")

        # ── Step 1: Extract and remap all needed fields ──────────────
        plevs = PRESSURE_LEVELS
        nlevs = len(plevs)
        ngrid = NLAT * NLON

        # Surface pressure (Pa)
        ps = extract_and_remap_surface(pgb_file, "PRES", "surface", tmpdir)

        # Surface geopotential height
        hgt_sfc = extract_and_remap_surface(pgb_file, "HGT", "surface", tmpdir)

        # Create GRIB1 template from one of the remapped files
        template_path = get_grib1_template(tmpdir)
        # Use the PRES surface remap as template
        pres_remap = os.path.join(tmpdir, "PRES_sfc_remap.grb1")
        if os.path.exists(pres_remap):
            shutil.copy2(pres_remap, template_path)
        else:
            # Fallback: create template from any remapped file
            for f in Path(tmpdir).glob("*_remap.grb1"):
                shutil.copy2(f, template_path)
                break

        # Pressure-level fields: shape (nlevs, ngrid)
        fields = {}
        for varname, _kpds5 in GRBSANL_SIGMA_VARS:
            data = np.zeros((nlevs, ngrid))
            for k, plev in enumerate(plevs):
                try:
                    data[k, :] = extract_and_remap_pressure_field(
                        pgb_file, varname, plev, tmpdir
                    )
                except RuntimeError as e:
                    log.warning(f"Missing {varname} at {plev} mb: {e}")
                    data[k, :] = np.nan
            fields[varname] = data

        # ── Step 2: Vertical interpolation to sigma levels ───────────
        sigma_fields = {}
        for varname, _kpds5 in GRBSANL_SIGMA_VARS:
            sigma_fields[varname] = vertical_interp_log_pressure(
                fields[varname], ps, plevs, SIGMA_LEVELS
            )

        # ── Step 3: Compute derived fields ───────────────────────────
        lats = np.array(GAUSS_LATS_94)
        nsigma = len(SIGMA_LEVELS)

        relv_all = np.zeros((nsigma, ngrid))
        reld_all = np.zeros((nsigma, ngrid))
        for k in range(nsigma):
            u_2d = sigma_fields["UGRD"][k, :].reshape(NLAT, NLON)
            v_2d = sigma_fields["VGRD"][k, :].reshape(NLAT, NLON)
            rv, rd = compute_vorticity_divergence(u_2d, v_2d, lats, NLON)
            relv_all[k, :] = rv.ravel()
            reld_all[k, :] = rd.ravel()

        # Terrain height gradients
        hgt_2d = hgt_sfc.reshape(NLAT, NLON)
        hgtx, hgty = compute_gradients(hgt_2d, lats, NLON)

        # Log-pressure gradients
        lnps_2d = np.log(np.maximum(ps.reshape(NLAT, NLON), 1.0))
        lpsx, lpsy = compute_gradients(lnps_2d, lats, NLON)

        # ── Step 4: Write all GRIB1 records ──────────────────────────
        # GDEX record order: RELV, RELD, TMP, SPFH, LPSX, LPSY,
        #                     UGRD, VGRD, PRES, HGT, HGTX, HGTY
        record_files = []

        # Derived fields on sigma levels: RELV, RELD
        for varname, kpds5 in GRBSANL_DERIVED_VARS.items():
            data = relv_all if varname == "RELV" else reld_all
            for k, sigma_k7 in enumerate(SIGMA_KPDS7):
                recpath = os.path.join(tmpdir, f"{varname}_{sigma_k7}.grb1")
                write_grib1_record(
                    data[k, :], template_path, recpath,
                    kpds5=kpds5, kpds6=107, kpds7=sigma_k7,
                    year=year, month=month, day=day, hour=hour,
                )
                record_files.append(recpath)

        # Interpolated fields on sigma levels (ordered: TMP, SPFH, UGRD, VGRD)
        # LPSX/LPSY inserted after SPFH to match GDEX
        for varname, kpds5 in GRBSANL_SIGMA_VARS:
            for k, sigma_k7 in enumerate(SIGMA_KPDS7):
                recpath = os.path.join(tmpdir, f"{varname}_{sigma_k7}.grb1")
                write_grib1_record(
                    sigma_fields[varname][k, :], template_path, recpath,
                    kpds5=kpds5, kpds6=107, kpds7=sigma_k7,
                    year=year, month=month, day=day, hour=hour,
                )
                record_files.append(recpath)

            # Insert LPSX, LPSY after SPFH
            if varname == "SPFH":
                for sfc_name, (sk5, sk6, sk7) in [
                    ("LPSX", GRBSANL_SURFACE_VARS["LPSX"]),
                    ("LPSY", GRBSANL_SURFACE_VARS["LPSY"]),
                ]:
                    sfc_data = lpsx.ravel() if sfc_name == "LPSX" else lpsy.ravel()
                    recpath = os.path.join(tmpdir, f"{sfc_name}.grb1")
                    write_grib1_record(
                        sfc_data, template_path, recpath,
                        kpds5=sk5, kpds6=sk6, kpds7=sk7,
                        year=year, month=month, day=day, hour=hour,
                    )
                    record_files.append(recpath)

        # Surface fields: PRES, HGT, HGTX, HGTY
        surface_data = {
            "PRES": ps,
            "HGT": hgt_sfc,
            "HGTX": hgtx.ravel(),
            "HGTY": hgty.ravel(),
        }
        for sfc_name in ["PRES", "HGT", "HGTX", "HGTY"]:
            kpds5, kpds6, kpds7 = GRBSANL_SURFACE_VARS[sfc_name]
            recpath = os.path.join(tmpdir, f"{sfc_name}_sfc.grb1")
            write_grib1_record(
                surface_data[sfc_name], template_path, recpath,
                kpds5=kpds5, kpds6=kpds6, kpds7=kpds7,
                year=year, month=month, day=day, hour=hour,
            )
            record_files.append(recpath)

        if not record_files:
            log.error(f"No records produced for {outname}")
            return None

        # ── Step 5: Concatenate with Fortran unformatted I/O framing ─
        # Each record: 4-byte BE length header + GRIB data + 4-byte BE length trailer
        with open(outpath, "wb") as outf:
            for rf in record_files:
                with open(rf, "rb") as inf:
                    data = inf.read()
                outf.write(struct.pack(">I", len(data)))
                outf.write(data)
                outf.write(struct.pack(">I", len(data)))

        log.info(f"Wrote {outname} with {len(record_files)} records")
        return outpath

    finally:
        if cleanup and os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)


def translate_month(year, month, outdir=None):
    """Translate all 6-hourly timesteps for a month."""
    if outdir is None:
        outdir = str(WORK_DIR / "grbsanl" / f"{year}" / f"{month:02d}")
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
