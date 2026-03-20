#!/usr/bin/env python3
"""
Diagnostic plots for CORe → GDEX translated data (MATCH transport model inputs).

Produces zonal mean cross-sections, sigma-level maps, and surface pressure maps
from the grbsanl monthly NetCDF. Uses NCAR brand styling from DAVINCI-MONET.

Usage:
    python plot_diagnostics.py 2026 1                    # All plots
    python plot_diagnostics.py 2026 1 -o ~/plots         # Custom output dir
    python plot_diagnostics.py 2026 1 --input file.nc    # Custom input
"""

import argparse
import os
import sys
import logging
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import xarray as xr

from config import WORK_DIR, SIGMA_LEVELS

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# ── NCAR style (adapted from DAVINCI-MONET) ─────────────────────────────

NCAR_COLORS = {
    "space": "#011837",
    "dark_blue": "#00357A",
    "ncar_blue": "#0A5DDA",
    "aqua": "#00A2B4",
    "light_blue": "#CEDFF8",
    "light_gray": "#F1F0EE",
    "orange": "#FF8C00",
    "yellow": "#FFDD31",
    "gray": "#58595B",
    "red": "#D62839",
    "green": "#2E8B57",
    "purple": "#7B68EE",
}

NCAR_PALETTE = [
    NCAR_COLORS["ncar_blue"], NCAR_COLORS["aqua"], NCAR_COLORS["orange"],
    NCAR_COLORS["purple"], NCAR_COLORS["green"], NCAR_COLORS["red"],
    NCAR_COLORS["yellow"], NCAR_COLORS["dark_blue"],
]


def apply_ncar_style():
    """Apply NCAR brand styling globally."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Poppins", "Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.titlesize": 16,
        "axes.prop_cycle": plt.cycler(color=NCAR_PALETTE),
        "lines.linewidth": 1.5,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "-",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


def save_figure(fig, name, output_dir):
    """Save figure as PNG (300 DPI) and PDF (rasterized base layer)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    png_path = output_dir / f"{name}.png"
    pdf_path = output_dir / f"{name}.pdf"

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    log.info(f"Saved: {png_path}")
    log.info(f"Saved: {pdf_path}")


# ── Plot functions ───────────────────────────────────────────────────────

def plot_zonal_cross_section(ds, varname, title, units, cmap, vmin, vmax,
                             output_dir, filename, nlevels=20):
    """
    Zonal mean cross-section: latitude (x) vs sigma (y, inverted).
    Monthly mean, averaged over longitude and time.
    """
    data = ds[varname].mean(dim=["time", "longitude"])
    lats = ds.latitude.values
    sigma = ds.sigma.values

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)

    levels = np.linspace(vmin, vmax, nlevels + 1)
    cf = ax.contourf(lats, sigma, data.values, levels=levels, cmap=cmap, extend="both")
    cf.set_rasterized(True)
    cs = ax.contour(lats, sigma, data.values, levels=levels[::2],
                    colors="k", linewidths=0.4, alpha=0.5)
    ax.clabel(cs, cs.levels[::2], fontsize=7, fmt="%.4g")

    ax.set_ylim(1.0, 0.0)
    ax.set_xlim(-90, 90)
    ax.set_xlabel("Latitude")
    ax.set_ylabel("Sigma")
    ax.set_title(title)

    cbar = fig.colorbar(cf, ax=ax, pad=0.02, aspect=30)
    cbar.set_label(units)

    save_figure(fig, filename, output_dir)
    plt.close(fig)


def plot_sigma_level_map(ds, varname, sigma_val, title, units, cmap, vmin, vmax,
                         output_dir, filename, nlevels=20,
                         overlay_winds=False, ds_u=None, ds_v=None):
    """
    Map at a specific sigma level. Monthly mean over time.
    Optionally overlay wind vectors.
    """
    # Find nearest sigma level
    sigma_idx = int(np.argmin(np.abs(np.array(SIGMA_LEVELS) - sigma_val)))
    actual_sigma = SIGMA_LEVELS[sigma_idx]

    data = ds[varname].mean(dim="time").isel(sigma=sigma_idx)
    lats = ds.latitude.values
    lons = ds.longitude.values

    # Add cyclic point to close the longitude seam
    data_cyc, lons_cyc = add_cyclic_point(data.values, coord=lons)

    fig = plt.figure(figsize=(12, 5), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

    levels = np.linspace(vmin, vmax, nlevels + 1)
    cf = ax.contourf(lons_cyc, lats, data_cyc, levels=levels, cmap=cmap,
                     extend="both", transform=ccrs.PlateCarree())
    cf.set_rasterized(True)

    if overlay_winds and ds_u is not None and ds_v is not None:
        u = ds_u["ugrd"].mean(dim="time").isel(sigma=sigma_idx).values
        v = ds_v["vgrd"].mean(dim="time").isel(sigma=sigma_idx).values
        # Subsample for legibility
        skip_lat = max(1, len(lats) // 20)
        skip_lon = max(1, len(lons) // 30)
        ax.quiver(lons[::skip_lon], lats[::skip_lat],
                  u[::skip_lat, ::skip_lon], v[::skip_lat, ::skip_lon],
                  transform=ccrs.PlateCarree(), scale=400, width=0.002,
                  color=NCAR_COLORS["gray"], alpha=0.7, zorder=5)

    ax.coastlines(linewidth=0.5, color=NCAR_COLORS["gray"])
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor=NCAR_COLORS["gray"])
    ax.set_global()
    ax.set_title(f"{title} (σ={actual_sigma:.4f})")

    cbar = fig.colorbar(cf, ax=ax, orientation="horizontal", pad=0.05,
                        shrink=0.7, aspect=35)
    cbar.set_label(units)

    save_figure(fig, filename, output_dir)
    plt.close(fig)


def plot_surface_map(ds, varname, title, units, cmap, vmin, vmax,
                     output_dir, filename, nlevels=20, scale_factor=1.0):
    """Map of a surface field. Monthly mean over time."""
    data = ds[varname].mean(dim="time") * scale_factor
    lats = ds.latitude.values
    lons = ds.longitude.values

    # Add cyclic point to close the longitude seam
    data_cyc, lons_cyc = add_cyclic_point(data.values, coord=lons)

    fig = plt.figure(figsize=(12, 5), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

    levels = np.linspace(vmin, vmax, nlevels + 1)
    cf = ax.contourf(lons_cyc, lats, data_cyc, levels=levels, cmap=cmap,
                     extend="both", transform=ccrs.PlateCarree())
    cf.set_rasterized(True)

    ax.coastlines(linewidth=0.5, color=NCAR_COLORS["gray"])
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor=NCAR_COLORS["gray"])
    ax.set_global()
    ax.set_title(title)

    cbar = fig.colorbar(cf, ax=ax, orientation="horizontal", pad=0.05,
                        shrink=0.7, aspect=35)
    cbar.set_label(units)

    save_figure(fig, filename, output_dir)
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────

def plot_diff_zonal_cross_section(ds_test, ds_ref, varname, title, units, vmax_diff,
                                  output_dir, filename, nlevels=20):
    """
    Zonal mean cross-section of (test - reference) difference.
    """
    diff = (ds_test[varname] - ds_ref[varname]).mean(dim=["time", "longitude"])
    lats = ds_test.latitude.values
    sigma = ds_test.sigma.values

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)

    levels = np.linspace(-vmax_diff, vmax_diff, nlevels + 1)
    cf = ax.contourf(lats, sigma, diff.values, levels=levels, cmap="RdBu_r", extend="both")
    cf.set_rasterized(True)
    cs = ax.contour(lats, sigma, diff.values, levels=levels[::2],
                    colors="k", linewidths=0.4, alpha=0.5)

    ax.set_ylim(1.0, 0.0)
    ax.set_xlim(-90, 90)
    ax.set_xlabel("Latitude")
    ax.set_ylabel("Sigma")
    ax.set_title(title)

    cbar = fig.colorbar(cf, ax=ax, pad=0.02, aspect=30)
    cbar.set_label(units)

    save_figure(fig, filename, output_dir)
    plt.close(fig)


def plot_diff_sigma_level_map(ds_test, ds_ref, varname, sigma_val, title, units, vmax_diff,
                              output_dir, filename, nlevels=20):
    """
    Map of (test - reference) difference at a specific sigma level.
    """
    sigma_idx = int(np.argmin(np.abs(np.array(SIGMA_LEVELS) - sigma_val)))
    actual_sigma = SIGMA_LEVELS[sigma_idx]

    diff = (ds_test[varname] - ds_ref[varname]).mean(dim="time").isel(sigma=sigma_idx)
    lats = ds_test.latitude.values
    lons = ds_test.longitude.values
    diff_cyc, lons_cyc = add_cyclic_point(diff.values, coord=lons)

    fig = plt.figure(figsize=(12, 5), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

    levels = np.linspace(-vmax_diff, vmax_diff, nlevels + 1)
    cf = ax.contourf(lons_cyc, lats, diff_cyc, levels=levels, cmap="RdBu_r",
                     extend="both", transform=ccrs.PlateCarree())
    cf.set_rasterized(True)

    ax.coastlines(linewidth=0.5, color=NCAR_COLORS["gray"])
    ax.set_global()
    ax.set_title(f"{title} (σ={actual_sigma:.4f})")

    cbar = fig.colorbar(cf, ax=ax, orientation="horizontal", pad=0.05,
                        shrink=0.7, aspect=35)
    cbar.set_label(units)

    save_figure(fig, filename, output_dir)
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────

def generate_all_plots(year, month, input_path=None, ref_path=None, output_dir=None):
    """Generate all diagnostic plots for a month of grbsanl data."""
    if input_path is None:
        input_path = WORK_DIR / f"grbsanl_{year}{month:02d}.nc"
    if ref_path is None:
        ref_path = WORK_DIR / f"grbsanl_ref_{year}{month:02d}.nc"
    if output_dir is None:
        output_dir = WORK_DIR / "plots" / f"{year}{month:02d}"

    log.info(f"Reading {input_path}")
    ds = xr.open_dataset(input_path)

    has_ref = Path(ref_path).exists()
    if has_ref:
        log.info(f"Reading reference {ref_path}")
        ds_ref = xr.open_dataset(ref_path)
    else:
        log.warning(f"No reference file at {ref_path}, skipping diff plots")

    label = f"CORe→GDEX Jan {year}"

    # ── Zonal mean cross-sections ────────────────────────────────────
    log.info("Plotting zonal mean cross-sections...")

    plot_zonal_cross_section(
        ds, "tmp",
        title=f"Zonal Mean Temperature — {label}",
        units="K", cmap="RdYlBu_r", vmin=190, vmax=300,
        output_dir=output_dir, filename="zonal_T",
    )

    plot_zonal_cross_section(
        ds, "ugrd",
        title=f"Zonal Mean Zonal Wind (U) — {label}",
        units="m/s", cmap="RdBu_r", vmin=-40, vmax=40,
        output_dir=output_dir, filename="zonal_U",
    )

    plot_zonal_cross_section(
        ds, "spfh",
        title=f"Zonal Mean Specific Humidity — {label}",
        units="kg/kg", cmap="YlGnBu", vmin=0, vmax=0.016,
        output_dir=output_dir, filename="zonal_Q",
    )

    plot_zonal_cross_section(
        ds, "vgrd",
        title=f"Zonal Mean Meridional Wind (V) — {label}",
        units="m/s", cmap="RdBu_r", vmin=-4, vmax=4,
        output_dir=output_dir, filename="zonal_V",
    )

    plot_zonal_cross_section(
        ds, "relv",
        title=f"Zonal Mean Relative Vorticity — {label}",
        units="s⁻¹", cmap="RdBu_r", vmin=-3e-5, vmax=3e-5,
        output_dir=output_dir, filename="zonal_RELV",
    )

    plot_zonal_cross_section(
        ds, "reld",
        title=f"Zonal Mean Divergence — {label}",
        units="s⁻¹", cmap="RdBu_r", vmin=-1e-6, vmax=1e-6,
        output_dir=output_dir, filename="zonal_RELD",
    )

    # ── Sigma-level maps: boundary layer (σ≈0.85) ───────────────────
    log.info("Plotting sigma-level maps (boundary layer)...")

    plot_sigma_level_map(
        ds, "tmp", sigma_val=0.8458,
        title=f"Temperature — {label}",
        units="K", cmap="RdYlBu_r", vmin=230, vmax=310,
        output_dir=output_dir, filename="map_T_sigma085",
        overlay_winds=True, ds_u=ds, ds_v=ds,
    )

    plot_sigma_level_map(
        ds, "spfh", sigma_val=0.8458,
        title=f"Specific Humidity — {label}",
        units="kg/kg", cmap="YlGnBu", vmin=0, vmax=0.016,
        output_dir=output_dir, filename="map_Q_sigma085",
    )

    # ── Sigma-level maps: jet level (σ≈0.25) ────────────────────────
    log.info("Plotting sigma-level maps (jet level)...")

    plot_sigma_level_map(
        ds, "tmp", sigma_val=0.2582,
        title=f"Temperature — {label}",
        units="K", cmap="RdYlBu_r", vmin=210, vmax=240,
        output_dir=output_dir, filename="map_T_sigma025",
        overlay_winds=True, ds_u=ds, ds_v=ds,
    )

    plot_sigma_level_map(
        ds, "ugrd", sigma_val=0.2582,
        title=f"Zonal Wind (U) — {label}",
        units="m/s", cmap="RdBu_r", vmin=-60, vmax=60,
        output_dir=output_dir, filename="map_U_sigma025",
    )

    # ── Surface pressure ─────────────────────────────────────────────
    log.info("Plotting surface pressure...")

    plot_surface_map(
        ds, "pres",
        title=f"Surface Pressure — {label}",
        units="hPa", cmap="viridis", vmin=500, vmax=1050,
        output_dir=output_dir, filename="map_Ps",
        scale_factor=0.01,  # Pa → hPa
    )

    plot_surface_map(
        ds, "hgt",
        title=f"Surface Geopotential Height — {label}",
        units="m", cmap="terrain", vmin=-200, vmax=5000,
        output_dir=output_dir, filename="map_HGT",
    )

    # ── Difference plots vs CDAS reference ───────────────────────────
    if has_ref:
        log.info("Plotting differences vs CDAS reference...")

        plot_diff_zonal_cross_section(
            ds, ds_ref, "tmp",
            title=f"Zonal Mean ΔT (CORe − CDAS) — Jan {year}",
            units="K", vmax_diff=5,
            output_dir=output_dir, filename="diff_zonal_T",
        )

        plot_diff_zonal_cross_section(
            ds, ds_ref, "ugrd",
            title=f"Zonal Mean ΔU (CORe − CDAS) — Jan {year}",
            units="m/s", vmax_diff=8,
            output_dir=output_dir, filename="diff_zonal_U",
        )

        plot_diff_zonal_cross_section(
            ds, ds_ref, "spfh",
            title=f"Zonal Mean ΔQ (CORe − CDAS) — Jan {year}",
            units="kg/kg", vmax_diff=0.002,
            output_dir=output_dir, filename="diff_zonal_Q",
        )

        plot_diff_zonal_cross_section(
            ds, ds_ref, "relv",
            title=f"Zonal Mean ΔRELV (CORe − CDAS) — Jan {year}",
            units="s⁻¹", vmax_diff=2e-5,
            output_dir=output_dir, filename="diff_zonal_RELV",
        )

        plot_diff_zonal_cross_section(
            ds, ds_ref, "reld",
            title=f"Zonal Mean ΔRELD (CORe − CDAS) — Jan {year}",
            units="s⁻¹", vmax_diff=5e-7,
            output_dir=output_dir, filename="diff_zonal_RELD",
        )

        plot_diff_sigma_level_map(
            ds, ds_ref, "tmp", sigma_val=0.8458,
            title=f"ΔT (CORe − CDAS) — Jan {year}",
            units="K", vmax_diff=8,
            output_dir=output_dir, filename="diff_map_T_sigma085",
        )

        plot_diff_sigma_level_map(
            ds, ds_ref, "ugrd", sigma_val=0.2582,
            title=f"ΔU (CORe − CDAS) — Jan {year}",
            units="m/s", vmax_diff=15,
            output_dir=output_dir, filename="diff_map_U_sigma025",
        )

        plot_diff_sigma_level_map(
            ds, ds_ref, "spfh", sigma_val=0.8458,
            title=f"ΔQ (CORe − CDAS) — Jan {year}",
            units="kg/kg", vmax_diff=0.003,
            output_dir=output_dir, filename="diff_map_Q_sigma085",
        )

        ds_ref.close()

    ds.close()
    log.info(f"All plots saved to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="Diagnostic plots for MATCH inputs")
    parser.add_argument("year", type=int)
    parser.add_argument("month", type=int)
    parser.add_argument("--input", help="Input NetCDF path")
    parser.add_argument("--ref", help="Reference NetCDF path")
    parser.add_argument("-o", "--output", help="Output directory")
    args = parser.parse_args()

    apply_ncar_style()
    generate_all_plots(args.year, args.month, args.input, args.ref, args.output)


if __name__ == "__main__":
    main()
