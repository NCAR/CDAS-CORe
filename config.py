"""
Configuration for CORe → GDEX translation pipeline.

Paths, constants, variable mappings, and helper functions for
translating CPC CORe reanalysis (GRIB2, 512x256 Gaussian, 3-hourly)
to GDEX/NCAR archive format (GRIB1, 192x94 Gaussian, 6-hourly).
"""

import os
import subprocess
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────

CORE_DIR = Path(os.path.expanduser("~/Data/CORe"))
GDEX_DIR = Path(os.path.expanduser("~/Data/GDEX"))
WORK_DIR = Path(os.path.expanduser("~/Data/work"))
REPO_DIR = Path(__file__).parent
TARGET_GRID = REPO_DIR / "target_grid_192x94.txt"

CONDA_ENV = "sarb"

# ── GRIB1 metadata constants ──────────────────────────────────────────────

CENTER = 7        # NCEP
SUBCENTER = 1     # NCEP sub
PROCESS = 180     # Generating process ID
TABLE_VERSION = 2 # GRIB1 parameter table version

# ── 28 CDAS sigma levels ─────────────────────────────────────────────────

SIGMA_LEVELS = [
    0.9950, 0.9821, 0.9644, 0.9425, 0.9159, 0.8838, 0.8458, 0.8014,
    0.7508, 0.6943, 0.6329, 0.5681, 0.5017, 0.4357, 0.3720, 0.3125,
    0.2582, 0.2101, 0.1682, 0.1326, 0.1028, 0.0782, 0.0580, 0.0418,
    0.0288, 0.0183, 0.0101, 0.0027,
]

# GRIB1 kpds7 encoding for sigma: sigma * 10000, as integer (round to match GDEX)
SIGMA_KPDS7 = [round(s * 10000) for s in SIGMA_LEVELS]

# CORe pressure levels in pgb files (mb), sorted descending for interpolation
PRESSURE_LEVELS = [
    1000, 925, 850, 800, 750, 700, 600, 500, 400, 300,
    250, 200, 150, 100, 70, 50, 30, 20, 10, 5, 2, 1,
]

# ── grb2d variable mappings (CORe flx → GDEX grb2d) ─────────────────────
#
# Each entry is a dict with:
#   core_name: wgrib2 variable name in CORe flx
#   core_level: wgrib2 level string to match
#   core_type: "anl" (instantaneous) or "ave" (0-3 hour average) or
#              "max" / "min" / "acc" (3-hour extremes/accumulations)
#   kpds5: GRIB1 parameter number
#   kpds6: GRIB1 level type
#   kpds7: GRIB1 level value
#   time_type: "inst" (TR=10) or "avg" (TR=3) in output
#
# Omitted from GDEX: CDCON, CFNLF, CFNSF, SGLYR, RHCLD, SRWEQ,
# cloud top/bot PRES (kpds6=212-233), cloud top TMP (kpds6=213/223/233)
# CPOFP substitutes for SRWEQ (kpds5=64)

GRB2D_VARS = [
    # ── Radiation (time-averaged) ──
    dict(core_name="DLWRF", core_level="surface", core_type="ave",
         kpds5=205, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="ULWRF", core_level="surface", core_type="ave",
         kpds5=212, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="ULWRF", core_level="top of atmosphere", core_type="ave",
         kpds5=212, kpds6=8, kpds7=0, time_type="avg"),
    dict(core_name="DSWRF", core_level="surface", core_type="ave",
         kpds5=204, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="DSWRF", core_level="top of atmosphere", core_type="ave",
         kpds5=204, kpds6=8, kpds7=0, time_type="avg"),
    dict(core_name="USWRF", core_level="surface", core_type="ave",
         kpds5=211, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="USWRF", core_level="top of atmosphere", core_type="ave",
         kpds5=211, kpds6=8, kpds7=0, time_type="avg"),
    # ── Clear-sky radiation (time-averaged) ──
    dict(core_name="CSDLF", core_level="surface", core_type="ave",
         kpds5=163, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="CSDSF", core_level="surface", core_type="ave",
         kpds5=161, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="CSULF", core_level="top of atmosphere", core_type="ave",
         kpds5=162, kpds6=8, kpds7=0, time_type="avg"),
    dict(core_name="CSUSF", core_level="surface", core_type="ave",
         kpds5=160, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="CSUSF", core_level="top of atmosphere", core_type="ave",
         kpds5=160, kpds6=8, kpds7=0, time_type="avg"),
    # ── Shortwave components (time-averaged) ──
    dict(core_name="NBDSF", core_level="surface", core_type="ave",
         kpds5=168, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="NDDSF", core_level="surface", core_type="ave",
         kpds5=169, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="VBDSF", core_level="surface", core_type="ave",
         kpds5=166, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="VDDSF", core_level="surface", core_type="ave",
         kpds5=167, kpds6=1, kpds7=0, time_type="avg"),
    # ── Surface fluxes (time-averaged) ──
    dict(core_name="SHTFL", core_level="surface", core_type="ave",
         kpds5=122, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="LHTFL", core_level="surface", core_type="ave",
         kpds5=121, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="UFLX", core_level="surface", core_type="ave",
         kpds5=124, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="VFLX", core_level="surface", core_type="ave",
         kpds5=125, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="GFLUX", core_level="surface", core_type="ave",
         kpds5=155, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="SNOHF", core_level="surface", core_type="ave",
         kpds5=229, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="SBSNO", core_level="surface", core_type="ave",
         kpds5=230, kpds6=1, kpds7=0, time_type="avg"),  # SBSNO→SNOEV
    # ── Precipitation (time-averaged) ──
    dict(core_name="PRATE", core_level="surface", core_type="ave",
         kpds5=59, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="CPRAT", core_level="surface", core_type="ave",
         kpds5=214, kpds6=1, kpds7=0, time_type="avg"),
    # ── Gravity wave drag (time-averaged) ──
    dict(core_name="U-GWD", core_level="surface", core_type="ave",
         kpds5=147, kpds6=1, kpds7=0, time_type="avg"),
    dict(core_name="V-GWD", core_level="surface", core_type="ave",
         kpds5=148, kpds6=1, kpds7=0, time_type="avg"),
    # ── Cloud cover (time-averaged) ──
    dict(core_name="TCDC", core_level="high cloud layer", core_type="ave",
         kpds5=71, kpds6=234, kpds7=0, time_type="avg"),
    dict(core_name="TCDC", core_level="middle cloud layer", core_type="ave",
         kpds5=71, kpds6=224, kpds7=0, time_type="avg"),
    dict(core_name="TCDC", core_level="low cloud layer", core_type="ave",
         kpds5=71, kpds6=214, kpds7=0, time_type="avg"),
    # ── Cloud water (time-averaged) ──
    dict(core_name="CWORK", core_level="atmos col", core_type="ave",
         kpds5=146, kpds6=200, kpds7=0, time_type="inst"),
    # ── Temperature extremes ──
    dict(core_name="TMP", core_level="2 m above ground", core_type="max",
         kpds5=15, kpds6=105, kpds7=2, time_type="inst"),
    dict(core_name="TMP", core_level="2 m above ground", core_type="min",
         kpds5=16, kpds6=105, kpds7=2, time_type="inst"),
    # ── Instantaneous surface/near-surface fields ──
    dict(core_name="PRES", core_level="surface", core_type="anl",
         kpds5=1, kpds6=1, kpds7=0, time_type="inst"),
    dict(core_name="PWAT", core_level="atmos col", core_type="anl",
         kpds5=54, kpds6=200, kpds7=0, time_type="inst"),
    dict(core_name="ICEC", core_level="surface", core_type="anl",
         kpds5=91, kpds6=1, kpds7=0, time_type="inst"),
    dict(core_name="LAND", core_level="surface", core_type="anl",
         kpds5=81, kpds6=1, kpds7=0, time_type="inst"),
    dict(core_name="SFCR", core_level="surface", core_type="anl",
         kpds5=83, kpds6=1, kpds7=0, time_type="inst"),
    dict(core_name="WEASD", core_level="surface", core_type="anl",
         kpds5=65, kpds6=1, kpds7=0, time_type="inst"),
    dict(core_name="PEVPR", core_level="surface", core_type="anl",
         kpds5=145, kpds6=1, kpds7=0, time_type="inst"),
    dict(core_name="SSRUN", core_level="surface", core_type="acc",
         kpds5=90, kpds6=1, kpds7=0, time_type="inst"),  # SSRUN→RUNOF
    dict(core_name="TMP", core_level="surface", core_type="anl",
         kpds5=11, kpds6=1, kpds7=0, time_type="inst"),
    dict(core_name="TMP", core_level="2 m above ground", core_type="anl",
         kpds5=11, kpds6=105, kpds7=2, time_type="inst"),
    dict(core_name="SPFH", core_level="2 m above ground", core_type="anl",
         kpds5=51, kpds6=105, kpds7=2, time_type="inst"),
    dict(core_name="UGRD", core_level="10 m above ground", core_type="anl",
         kpds5=33, kpds6=105, kpds7=10, time_type="inst"),
    dict(core_name="VGRD", core_level="10 m above ground", core_type="anl",
         kpds5=34, kpds6=105, kpds7=10, time_type="inst"),
    # ── Soil (instantaneous) ──
    dict(core_name="SOILW", core_level="0-0.1 m below ground", core_type="anl",
         kpds5=144, kpds6=112, kpds7=10, time_type="inst"),
    dict(core_name="SOILW", core_level="0.1-0.4 m below ground", core_type="anl",
         kpds5=144, kpds6=112, kpds7=2760, time_type="inst"),
    dict(core_name="TSOIL", core_level="0-0.1 m below ground", core_type="anl",
         kpds5=11, kpds6=112, kpds7=10, time_type="inst"),
    dict(core_name="TSOIL", core_level="0.1-0.4 m below ground", core_type="anl",
         kpds5=11, kpds6=112, kpds7=2760, time_type="inst"),
    dict(core_name="TSOIL", core_level="1-2 m below ground", core_type="anl",
         kpds5=11, kpds6=111, kpds7=300, time_type="inst"),
    # ── CPOFP substituting for SRWEQ ──
    dict(core_name="CPOFP", core_level="surface", core_type="anl",
         kpds5=64, kpds6=1, kpds7=0, time_type="avg"),
]

# ── grbsanl variable mappings ────────────────────────────────────────────

# Variables interpolated from pressure levels to sigma levels
# (ordered to match GDEX record layout):
GRBSANL_SIGMA_VARS = [
    ("TMP", 11),
    ("SPFH", 51),
    ("UGRD", 33),
    ("VGRD", 34),
]

# Derived fields computed on sigma levels:
GRBSANL_DERIVED_VARS = {
    "RELV": 43,  # relative vorticity
    "RELD": 44,  # divergence
}

# Surface fields in grbsanl:
GRBSANL_SURFACE_VARS = {
    "PRES": (1, 1, 0),     # surface pressure
    "HGT":  (7, 1, 0),     # surface geopotential height
    "LPSX": (181, 1, 0),   # ln(Ps) x-gradient
    "LPSY": (182, 1, 0),   # ln(Ps) y-gradient
    "HGTX": (183, 1, 0),   # terrain height x-gradient
    "HGTY": (184, 1, 0),   # terrain height y-gradient
}

# ── 192x94 Gaussian grid parameters ─────────────────────────────────────

GAUSS_LATS_94 = [
    88.5420, 86.6532, 84.7532, 82.8508, 80.9474, 79.0435, 77.1393,
    75.2351, 73.3307, 71.4262, 69.5217, 67.6171, 65.7125, 63.8079,
    61.9033, 59.9986, 58.0940, 56.1893, 54.2846, 52.3799, 50.4752,
    48.5705, 46.6658, 44.7611, 42.8564, 40.9517, 39.0470, 37.1422,
    35.2375, 33.3328, 31.4281, 29.5234, 27.6186, 25.7139, 23.8092,
    21.9044, 19.9997, 18.0950, 16.1902, 14.2855, 12.3808, 10.4760,
    8.5713, 6.6666, 4.7618, 2.8571, 0.9524,
    -0.9524, -2.8571, -4.7618, -6.6666, -8.5713, -10.4760,
    -12.3808, -14.2855, -16.1902, -18.0950, -19.9997, -21.9044,
    -23.8092, -25.7139, -27.6186, -29.5234, -31.4281, -33.3328,
    -35.2375, -37.1422, -39.0470, -40.9517, -42.8564, -44.7611,
    -46.6658, -48.5705, -50.4752, -52.3799, -54.2846, -56.1893,
    -58.0940, -59.9986, -61.9033, -63.8079, -65.7125, -67.6171,
    -69.5217, -71.4262, -73.3307, -75.2351, -77.1393, -79.0435,
    -80.9474, -82.8508, -84.7532, -86.6532, -88.5420,
]

NLON = 192
NLAT = 94
DLON = 1.875

# Earth radius in meters
EARTH_RADIUS = 6.371e6


# ── Helper functions ─────────────────────────────────────────────────────

def run_cmd(cmd, check=True):
    """Run a shell command, returning CompletedProcess."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {cmd}\nstderr: {result.stderr}\nstdout: {result.stdout}"
        )
    return result


def wgrib2_extract(infile, match_pattern, outfile):
    """Extract GRIB2 records matching pattern via wgrib2."""
    cmd = f"wgrib2 {infile} -match '{match_pattern}' -grib {outfile}"
    return run_cmd(cmd)


def cdo_remap(infile, outfile, grid_file=None):
    """Bilinear remap to target grid and convert to GRIB1."""
    if grid_file is None:
        grid_file = TARGET_GRID
    cmd = f"cdo -f grb1 remapbil,{grid_file} {infile} {outfile}"
    return run_cmd(cmd)


def grib_set(infile, outfile, **keys):
    """Set GRIB1 keys via grib_set."""
    pairs = ",".join(f"{k}={v}" for k, v in keys.items())
    cmd = f"grib_set -s {pairs} {infile} {outfile}"
    return run_cmd(cmd)


def core_flx_path(year, month, day, hour):
    """Path to a CORe flx file."""
    return (CORE_DIR / "flx" / f"{year}" / f"{month:02d}" /
            f"flx.{year}{month:02d}{day:02d}{hour:02d}.grb")


def core_pgb_path(year, month, day, hour):
    """Path to a CORe pgb file."""
    return (CORE_DIR / "pgb" / f"{year}" / f"{month:02d}" /
            f"pgb.{year}{month:02d}{day:02d}{hour:02d}.grb")


def sixhourly_timestamps(year, month):
    """Generate all 6-hourly (year, month, day, hour) tuples for a month."""
    import calendar
    ndays = calendar.monthrange(year, month)[1]
    timestamps = []
    for day in range(1, ndays + 1):
        for hour in (0, 6, 12, 18):
            timestamps.append((year, month, day, hour))
    return timestamps
