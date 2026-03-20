#!/usr/bin/env python3
"""
Validate translated CORe output against GDEX reference data.

Reads GRIB1 records from both files, matches by variable identifier
(kpds5, kpds6, kpds7), and computes comparison statistics.

Usage:
    # Compare single timestep (grb2d)
    python validate.py grb2d 2026 1 1 0

    # Compare single timestep (grbsanl)
    python validate.py grbsanl 2026 1 1 0

    # Compare full month
    python validate.py grb2d 2026 1 --month
    python validate.py grbsanl 2026 1 --month

    # Verbose: print per-record stats
    python validate.py grb2d 2026 1 1 0 -v
"""

import argparse
import io
import os
import struct
import sys
import tarfile
import logging
from collections import OrderedDict
from pathlib import Path

import numpy as np

try:
    import eccodes
except ImportError:
    sys.exit("eccodes Python package required: pip install eccodes")

from config import (
    GDEX_DIR, WORK_DIR, sixhourly_timestamps,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)


# ── GRIB1 reading ────────────────────────────────────────────────────────

def read_grib1_records(filepath):
    """
    Read all GRIB1 records from a file.
    Returns list of dicts with keys: kpds5, kpds6, kpds7, tr, p1, p2, label, values
    """
    records = []
    with open(filepath, "rb") as f:
        while True:
            msgid = eccodes.codes_grib_new_from_file(f)
            if msgid is None:
                break
            rec = _extract_record(msgid)
            records.append(rec)
            eccodes.codes_release(msgid)
    return records


def read_grib1_records_prefixed(filepath):
    """
    Read GRIB1 records from a file with 4-byte BE length prefix per record
    (Fortran sequential format, used by grbsanl).
    """
    records = []
    with open(filepath, "rb") as f:
        data = f.read()

    pos = 0
    rec_num = 0
    while pos < len(data):
        if pos + 4 > len(data):
            break
        reclen = struct.unpack(">I", data[pos:pos + 4])[0]
        pos += 4
        if pos + reclen > len(data):
            break
        msg_data = data[pos:pos + reclen]
        pos += reclen
        rec_num += 1

        # Skip trailing length marker if present (Fortran unformatted I/O)
        if pos + 4 <= len(data):
            trailer = struct.unpack(">I", data[pos:pos + 4])[0]
            if trailer == reclen:
                pos += 4

        # Verify GRIB magic number
        if msg_data[:4] != b"GRIB":
            log.warning(f"Record {rec_num}: missing GRIB magic, skipping")
            continue

        try:
            msgid = eccodes.codes_new_from_message(msg_data)
            rec = _extract_record(msgid)
            records.append(rec)
            eccodes.codes_release(msgid)
        except Exception as e:
            log.warning(f"Record {rec_num}: failed to decode ({e}), skipping")

    return records


def _extract_record(msgid):
    """Extract metadata and values from an eccodes message handle."""
    kpds5 = eccodes.codes_get_long(msgid, "indicatorOfParameter")
    kpds6 = eccodes.codes_get_long(msgid, "indicatorOfTypeOfLevel")
    kpds7 = eccodes.codes_get_long(msgid, "level")
    tr = eccodes.codes_get_long(msgid, "timeRangeIndicator")
    p1 = eccodes.codes_get_long(msgid, "P1")
    p2 = eccodes.codes_get_long(msgid, "P2")
    values = eccodes.codes_get_values(msgid)

    # For layer types, reconstruct kpds7 from top/bottom
    if kpds6 in (112, 116):
        top = eccodes.codes_get_long(msgid, "topLevel")
        bot = eccodes.codes_get_long(msgid, "bottomLevel")
        kpds7 = top * 256 + bot

    # Build human-readable label
    try:
        name = eccodes.codes_get(msgid, "shortName")
    except Exception:
        name = f"p{kpds5}"
    label = f"{name}:k5={kpds5}:k6={kpds6}:k7={kpds7}:TR={tr}"

    return dict(
        kpds5=kpds5, kpds6=kpds6, kpds7=kpds7,
        tr=tr, p1=p1, p2=p2,
        label=label, values=values,
    )


def record_key(rec):
    """Unique key for matching records between files."""
    return (rec["kpds5"], rec["kpds6"], rec["kpds7"], rec["tr"])


# ── Reference data extraction ────────────────────────────────────────────

def find_gdex_archive(filetype, year, month):
    """Find the GDEX reference tar archive for a given month."""
    gdex_subdir = GDEX_DIR / filetype / str(year)
    if not gdex_subdir.is_dir():
        return None
    import re
    pattern = re.compile(rf"A\d+-{year}{month:02d}\.{filetype}$")
    for entry in sorted(os.listdir(gdex_subdir)):
        if pattern.match(entry):
            return gdex_subdir / entry
    return None


def extract_from_tar(archive_path, member_name, tmpdir="/tmp"):
    """Extract a single member from a tar archive, return path."""
    outpath = os.path.join(tmpdir, member_name)
    with tarfile.open(archive_path, "r") as tar:
        tar.extract(member_name, path=tmpdir)
    return outpath


# ── Comparison ────────────────────────────────────────────────────────────

def compare_records(test_records, ref_records, verbose=False):
    """
    Compare matched records between test and reference.
    Returns (stats_list, summary_dict).
    """
    # Index reference by key
    ref_by_key = OrderedDict()
    for rec in ref_records:
        key = record_key(rec)
        ref_by_key[key] = rec

    # Index test by key
    test_by_key = OrderedDict()
    for rec in test_records:
        key = record_key(rec)
        test_by_key[key] = rec

    matched = set(test_by_key.keys()) & set(ref_by_key.keys())
    test_only = set(test_by_key.keys()) - set(ref_by_key.keys())
    ref_only = set(ref_by_key.keys()) - set(test_by_key.keys())

    stats = []
    for key in sorted(matched):
        test_vals = test_by_key[key]["values"]
        ref_vals = ref_by_key[key]["values"]
        label = test_by_key[key]["label"]

        if len(test_vals) != len(ref_vals):
            stats.append(dict(
                label=label, key=key, error="size mismatch",
                test_size=len(test_vals), ref_size=len(ref_vals),
            ))
            continue

        diff = test_vals - ref_vals
        abs_diff = np.abs(diff)
        ref_range = np.ptp(ref_vals)  # range of reference values

        rmsd = np.sqrt(np.mean(diff ** 2))
        mae = np.mean(abs_diff)
        max_abs = np.max(abs_diff)
        mean_bias = np.mean(diff)

        # Correlation (handle constant fields)
        if np.std(test_vals) > 0 and np.std(ref_vals) > 0:
            corr = np.corrcoef(test_vals, ref_vals)[0, 1]
        else:
            corr = 1.0 if np.allclose(test_vals, ref_vals) else 0.0

        # Normalized RMSD (relative to reference range, skip if range ~0)
        nrmsd = rmsd / ref_range if ref_range > 1e-15 else 0.0

        s = dict(
            label=label, key=key,
            rmsd=rmsd, mae=mae, max_abs=max_abs, mean_bias=mean_bias,
            corr=corr, nrmsd=nrmsd, ref_range=ref_range,
            ref_mean=np.mean(ref_vals), test_mean=np.mean(test_vals),
        )
        stats.append(s)

        if verbose:
            print(f"  {label:50s}  RMSD={rmsd:12.5g}  MaxAbs={max_abs:12.5g}  "
                  f"NRMSD={nrmsd:.4f}  Corr={corr:.6f}  Bias={mean_bias:+.5g}")

    summary = dict(
        n_matched=len(matched),
        n_test_only=len(test_only),
        n_ref_only=len(ref_only),
        test_only_keys=sorted(test_only),
        ref_only_keys=sorted(ref_only),
    )

    # Aggregate stats across matched fields
    if stats and all("rmsd" in s for s in stats):
        nrmsds = [s["nrmsd"] for s in stats if "nrmsd" in s]
        corrs = [s["corr"] for s in stats if "corr" in s]
        summary["median_nrmsd"] = np.median(nrmsds) if nrmsds else None
        summary["min_corr"] = np.min(corrs) if corrs else None
        summary["max_nrmsd"] = np.max(nrmsds) if nrmsds else None
        summary["mean_corr"] = np.mean(corrs) if corrs else None

    return stats, summary


# ── Timestep comparison ──────────────────────────────────────────────────

def validate_timestep(filetype, year, month, day, hour, verbose=False):
    """
    Compare one translated timestep against GDEX reference.
    Returns (stats, summary) or None if files not found.
    """
    yy = year % 100
    fname = f"{filetype}{yy:02d}{month:02d}{day:02d}{hour:02d}"

    # Translated output
    test_path = WORK_DIR / filetype / str(year) / f"{month:02d}" / fname
    if not test_path.exists():
        log.error(f"Translated file not found: {test_path}")
        return None

    # Reference from GDEX tar
    archive = find_gdex_archive(filetype, year, month)
    if archive is None:
        log.error(f"No GDEX reference archive for {filetype} {year}-{month:02d}")
        return None

    ref_path = extract_from_tar(str(archive), fname)

    # Read records
    prefixed = (filetype == "grbsanl")
    reader = read_grib1_records_prefixed if prefixed else read_grib1_records

    test_records = reader(str(test_path))
    ref_records = reader(ref_path)

    # Clean up extracted reference
    os.unlink(ref_path)

    log.info(f"{fname}: {len(test_records)} test records, {len(ref_records)} ref records")

    stats, summary = compare_records(test_records, ref_records, verbose=verbose)
    return stats, summary


def validate_month(filetype, year, month, verbose=False):
    """Compare all translated timesteps for a month against GDEX reference."""
    timestamps = sixhourly_timestamps(year, month)
    all_stats = []
    summaries = []
    n_ok = 0
    n_fail = 0

    for y, m, d, h in timestamps:
        result = validate_timestep(filetype, y, m, d, h, verbose=verbose)
        if result is None:
            n_fail += 1
            continue
        stats, summary = result
        all_stats.extend(stats)
        summaries.append((f"{y}{m:02d}{d:02d}{h:02d}", summary))
        n_ok += 1

    log.info(f"Month {year}-{month:02d} {filetype}: {n_ok} timesteps compared, {n_fail} missing")

    # Aggregate monthly stats
    if all_stats:
        nrmsds = [s["nrmsd"] for s in all_stats if "nrmsd" in s]
        corrs = [s["corr"] for s in all_stats if "corr" in s]
        max_abs_all = [s["max_abs"] for s in all_stats if "max_abs" in s]

        print(f"\n{'='*70}")
        print(f"MONTHLY SUMMARY: {filetype} {year}-{month:02d}")
        print(f"{'='*70}")
        print(f"  Timesteps compared:   {n_ok}")
        print(f"  Total record pairs:   {len(all_stats)}")
        if nrmsds:
            print(f"  Median NRMSD:         {np.median(nrmsds):.6f}")
            print(f"  95th pctile NRMSD:    {np.percentile(nrmsds, 95):.6f}")
            print(f"  Max NRMSD:            {np.max(nrmsds):.6f}")
        if corrs:
            print(f"  Mean correlation:     {np.mean(corrs):.6f}")
            print(f"  Min correlation:      {np.min(corrs):.6f}")
        if max_abs_all:
            print(f"  Max abs difference:   {np.max(max_abs_all):.6g}")

        # Flag worst fields
        worst = sorted(
            [s for s in all_stats if "nrmsd" in s],
            key=lambda s: s["nrmsd"], reverse=True
        )[:10]
        if worst:
            print(f"\n  Worst 10 fields by NRMSD:")
            for s in worst:
                print(f"    {s['label']:50s}  NRMSD={s['nrmsd']:.4f}  "
                      f"Corr={s['corr']:.4f}  MaxAbs={s['max_abs']:.5g}")

        # Report unmatched records (from first timestep as representative)
        if summaries:
            _, first_summary = summaries[0]
            if first_summary["ref_only_keys"]:
                print(f"\n  Reference-only records (omitted vars, per timestep): "
                      f"{first_summary['n_ref_only']}")
            if first_summary["test_only_keys"]:
                print(f"  Test-only records (extra vars, per timestep): "
                      f"{first_summary['n_test_only']}")
                for key in first_summary["test_only_keys"][:5]:
                    print(f"    {key}")

        print(f"{'='*70}\n")

    return all_stats, summaries


# ── CLI ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Validate translated CORe output against GDEX reference"
    )
    parser.add_argument("filetype", choices=["grb2d", "grbsanl"],
                        help="File type to validate")
    parser.add_argument("year", type=int, help="4-digit year")
    parser.add_argument("month", type=int, help="Month (1-12)")
    parser.add_argument("day", type=int, nargs="?", help="Day (for single timestep)")
    parser.add_argument("hour", type=int, nargs="?", help="Hour (for single timestep)")
    parser.add_argument("--month", dest="full_month", action="store_true",
                        help="Compare all timesteps in the month")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print per-record statistics")
    args = parser.parse_args()

    if args.full_month:
        validate_month(args.filetype, args.year, args.month, verbose=args.verbose)
    elif args.day is not None and args.hour is not None:
        result = validate_timestep(
            args.filetype, args.year, args.month, args.day, args.hour,
            verbose=args.verbose,
        )
        if result is None:
            sys.exit(1)
        stats, summary = result
        print(f"\nMatched: {summary['n_matched']}  "
              f"Test-only: {summary['n_test_only']}  "
              f"Ref-only: {summary['n_ref_only']}")
        if summary.get("median_nrmsd") is not None:
            print(f"Median NRMSD: {summary['median_nrmsd']:.6f}  "
                  f"Min corr: {summary['min_corr']:.6f}")
        if summary["ref_only_keys"]:
            print(f"Reference-only (omitted vars): {summary['n_ref_only']} records")
        if summary["test_only_keys"]:
            print(f"Test-only (extra vars):")
            for key in summary["test_only_keys"]:
                print(f"  {key}")
    else:
        parser.error("Provide day and hour for single timestep, or use --month")


if __name__ == "__main__":
    main()
