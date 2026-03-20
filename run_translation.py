#!/usr/bin/env python3
"""
CLI entry point for CORe → GDEX translation pipeline.

Usage:
    python run_translation.py YYYY MM [--grb2d] [--grbsanl] [--package] [--all]

Examples:
    python run_translation.py 2026 1 --all          # Full pipeline for Jan 2026
    python run_translation.py 2026 1 --grb2d         # Only translate flx → grb2d
    python run_translation.py 2026 1 --grbsanl       # Only translate pgb → grbsanl
    python run_translation.py 2026 1 --package        # Only package (assumes translated files exist)
"""

import argparse
import logging
import sys

from translate_grb2d import translate_month as translate_grb2d_month
from translate_grbsanl import translate_month as translate_grbsanl_month
from package_gdex import package_month

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="CORe → GDEX translation pipeline"
    )
    parser.add_argument("year", type=int, help="4-digit year")
    parser.add_argument("month", type=int, help="Month (1-12)")
    parser.add_argument("--grb2d", action="store_true",
                        help="Translate CORe flx → GDEX grb2d")
    parser.add_argument("--grbsanl", action="store_true",
                        help="Translate CORe pgb → GDEX grbsanl")
    parser.add_argument("--package", action="store_true",
                        help="Package translated files into monthly tar archives")
    parser.add_argument("--all", action="store_true",
                        help="Run full pipeline (translate + package)")
    args = parser.parse_args()

    if not any([args.grb2d, args.grbsanl, args.package, args.all]):
        parser.print_help()
        sys.exit(1)

    do_grb2d = args.grb2d or args.all
    do_grbsanl = args.grbsanl or args.all
    do_package = args.package or args.all

    year, month = args.year, args.month

    if do_grb2d:
        log.info(f"=== Translating CORe flx → grb2d for {year}-{month:02d} ===")
        grb2d_files = translate_grb2d_month(year, month)
        log.info(f"grb2d: {len(grb2d_files)} timesteps translated")

    if do_grbsanl:
        log.info(f"=== Translating CORe pgb → grbsanl for {year}-{month:02d} ===")
        grbsanl_files = translate_grbsanl_month(year, month)
        log.info(f"grbsanl: {len(grbsanl_files)} timesteps translated")

    if do_package:
        log.info(f"=== Packaging monthly archives for {year}-{month:02d} ===")
        if do_grb2d or args.package:
            result = package_month(year, month, "grb2d")
            if result:
                log.info(f"grb2d archive: {result}")
        if do_grbsanl or args.package:
            result = package_month(year, month, "grbsanl")
            if result:
                log.info(f"grbsanl archive: {result}")

    log.info("Done.")


if __name__ == "__main__":
    main()
