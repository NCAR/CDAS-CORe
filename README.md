# CDAS-CORe

Translation pipeline for converting CPC CORe reanalysis data (GRIB2, 512x256 Gaussian, 3-hourly) to CDAS/GDEX archive format (GRIB1, 192x94 Gaussian, 6-hourly) for use as MATCH transport model inputs.

## Overview

The pipeline produces two output products per timestep:

- **grb2d** — 52-record surface/flux file translated from CORe `flx` files. Includes radiation, surface fluxes, precipitation, cloud cover, soil fields, and near-surface meteorology.
- **grbsanl** — 174-record sigma-level file translated from CORe `pgb` files. Includes T, Q, U, V interpolated from 22 pressure levels to 28 CDAS sigma levels, plus derived vorticity, divergence, and surface gradient fields.

## Setup

### Requirements

- Python >= 3.12
- [CDO](https://code.mpimet.mpg.de/projects/cdo) (Climate Data Operators)
- [wgrib2](https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/)
- [eccodes](https://confluence.ecmwf.int/display/ECC) Python bindings

### Install

```bash
conda env create -f environment.yml
conda activate sarb
```

### Data directories

Edit `config.py` to set paths for your system. Defaults:

| Variable   | Default              | Description                         |
|------------|----------------------|-------------------------------------|
| `CORE_DIR` | `~/Data/CORe`       | CORe input files (`flx/` and `pgb/` subdirs) |
| `GDEX_DIR` | `~/Data/GDEX`       | GDEX reference archives (for validation) |
| `WORK_DIR` | `~/Data/work`       | Working/output directory            |

CORe input files are expected at:
```
$CORE_DIR/flx/{YYYY}/{MM}/flx.{YYYYMMDD}{HH}.grb
$CORE_DIR/pgb/{YYYY}/{MM}/pgb.{YYYYMMDD}{HH}.grb
```

## Usage

### Translation

Run the full pipeline (translate + package) for a given month:

```bash
python run_translation.py 2026 1 --all
```

Or run individual stages:

```bash
python run_translation.py 2026 1 --grb2d      # Translate flx → grb2d only
python run_translation.py 2026 1 --grbsanl     # Translate pgb → grbsanl only
python run_translation.py 2026 1 --package      # Package into monthly tar archives
```

Translated files are written to `$WORK_DIR/{grb2d,grbsanl}/{YYYY}/{MM}/`.

### Validation

Compare translated output against GDEX reference data:

```bash
# Single timestep
python validate.py grb2d 2026 1 1 0
python validate.py grbsanl 2026 1 1 0 -v       # Verbose per-record stats

# Full month
python validate.py grb2d 2026 1 --month
python validate.py grbsanl 2026 1 --month
```

Reports matched/unmatched record counts, RMSD, NRMSD, correlation, and bias for each variable.

### NetCDF conversion

Convert translated GRIB1 files to CF-compliant NetCDF4:

```bash
python convert_netcdf.py grb2d 2026 1              # Full month
python convert_netcdf.py grbsanl 2026 1             # Full month
python convert_netcdf.py grb2d 2026 1 --day 1 --hour 0   # Single timestep
python convert_netcdf.py grb2d 2026 1 -o output.nc  # Custom output path
```

### Diagnostic plots

Generate zonal cross-sections, sigma-level maps, and surface pressure maps from grbsanl NetCDF output:

```bash
python plot_diagnostics.py 2026 1                    # All plots
python plot_diagnostics.py 2026 1 -o ~/plots         # Custom output dir
python plot_diagnostics.py 2026 1 --input file.nc    # Custom input NetCDF
python plot_diagnostics.py 2026 1 --ref ref.nc       # Include difference plots vs reference
```

Outputs PNG (300 DPI) and PDF for each plot.

## Pipeline details

### grb2d translation (`translate_grb2d.py`)

For each 6-hour window:
- **Instantaneous fields**: extracted from the synoptic-time CORe file
- **Time-averaged fields**: averaged across two 3-hr CORe files spanning the window
- **TMAX/TMIN**: max/min across the two 3-hr periods

Per-variable pipeline: wgrib2 extract → CDO bilinear remap (512x256 → 192x94) → grib_set metadata correction.

### grbsanl translation (`translate_grbsanl.py`)

Per 6-hourly timestep:
1. Extract U, V, T, Q, HGT, PRES from CORe `pgb` on pressure levels
2. Bilinear spatial remap (512x256 → 192x94 Gaussian via CDO)
3. Log-pressure vertical interpolation from 22 pressure levels to 28 sigma levels
4. Compute derived fields (vorticity, divergence, terrain and ln(Ps) gradients)
5. Write GRIB1 with Fortran sequential (4-byte length-prefixed) record format

## Project structure

```
config.py               Configuration, paths, variable mappings, grid constants
run_translation.py      CLI entry point for the translation pipeline
translate_grb2d.py      CORe flx → GDEX grb2d translation
translate_grbsanl.py    CORe pgb → GDEX grbsanl translation
package_gdex.py         Package translated files into monthly tar archives
validate.py             Validate translated output against GDEX reference data
convert_netcdf.py       Convert translated GRIB1 to CF-compliant NetCDF4
plot_diagnostics.py     Diagnostic plots (zonal cross-sections, maps)
target_grid_192x94.txt  CDO grid description for the 192x94 Gaussian target
environment.yml         Conda environment specification
```
