#!/usr/bin/env python3
"""
Package translated grb2d/grbsanl files into monthly tar archives
matching GDEX naming conventions.

Archive naming: A#####-YYYYMM.grb2d / A#####-YYYYMM.grbsanl
Inner file naming: grb2dYYMMDDHH / grbsanlYYMMDDHH
A##### numbering continues from the last existing archive in GDEX_DIR.
"""

import os
import sys
import glob
import tarfile
import logging
import re
from pathlib import Path

from config import GDEX_DIR, WORK_DIR, sixhourly_timestamps

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)


def find_last_archive_number(gdex_subdir):
    """
    Find the highest A##### number in existing GDEX archives.
    Returns the integer portion (e.g., 27337 from A27337-202602.grb2d).
    """
    max_num = 0
    pattern = re.compile(r"A(\d+)-")
    if gdex_subdir.is_dir():
        for entry in os.listdir(gdex_subdir):
            m = pattern.match(entry)
            if m:
                max_num = max(max_num, int(m.group(1)))
    # Also search subdirectories (year folders)
    for yeardir in gdex_subdir.glob("*"):
        if yeardir.is_dir():
            for entry in os.listdir(yeardir):
                m = pattern.match(entry)
                if m:
                    max_num = max(max_num, int(m.group(1)))
    return max_num


def package_month(year, month, filetype, indir=None, outdir=None):
    """
    Package all 6-hourly files for a month into a tar archive.

    Args:
        year: 4-digit year
        month: month (1-12)
        filetype: "grb2d" or "grbsanl"
        indir: directory containing translated files (default: WORK_DIR/filetype/year/month)
        outdir: directory for output tar (default: GDEX_DIR/filetype/year/)

    Returns:
        Path to created tar archive, or None on failure.
    """
    if filetype not in ("grb2d", "grbsanl"):
        raise ValueError(f"filetype must be 'grb2d' or 'grbsanl', got {filetype}")

    if indir is None:
        indir = WORK_DIR / filetype / f"{year}" / f"{month:02d}"
    else:
        indir = Path(indir)

    if outdir is None:
        outdir = GDEX_DIR / filetype / f"{year}"
    else:
        outdir = Path(outdir)
    os.makedirs(outdir, exist_ok=True)

    # Find expected files
    yy = year % 100
    timestamps = sixhourly_timestamps(year, month)
    expected_files = []
    for y, m, d, h in timestamps:
        fname = f"{filetype}{yy:02d}{m:02d}{d:02d}{h:02d}"
        fpath = indir / fname
        if fpath.exists():
            expected_files.append((fname, fpath))
        else:
            log.warning(f"Missing file: {fpath}")

    if not expected_files:
        log.error(f"No files found in {indir}")
        return None

    # Determine A##### number
    gdex_subdir = GDEX_DIR / filetype
    last_num = find_last_archive_number(gdex_subdir)
    new_num = last_num + 1

    archive_name = f"A{new_num:05d}-{year}{month:02d}.{filetype}"
    archive_path = outdir / archive_name

    # Create tar archive
    with tarfile.open(archive_path, "w") as tar:
        for fname, fpath in sorted(expected_files):
            tar.add(str(fpath), arcname=fname)

    log.info(f"Created {archive_path} with {len(expected_files)} files")
    return str(archive_path)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} YYYY MM grb2d|grbsanl [indir] [outdir]")
        sys.exit(1)
    year = int(sys.argv[1])
    month = int(sys.argv[2])
    filetype = sys.argv[3]
    indir = sys.argv[4] if len(sys.argv) > 4 else None
    outdir = sys.argv[5] if len(sys.argv) > 5 else None
    package_month(year, month, filetype, indir, outdir)
