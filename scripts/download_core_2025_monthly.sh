#!/usr/bin/env bash
set -euo pipefail

BASE_URL="${BASE_URL:-https://storage.googleapis.com/noaa-nws-ncep-core/grib/month}"
OUT_ROOT="${OUT_ROOT:-$HOME/Data/CORe}"
YEAR=2025

usage() {
  cat <<'EOF'
Usage: download_core_2025_monthly.sh [flx] [pgb]

Download all monthly 2025 CORe archives into:
  $OUT_ROOT/month/{flx,pgb}/2025/{flx,pgb}.2025MM.grb

Environment:
  OUT_ROOT   Destination root. Default: $HOME/Data/CORe
  BASE_URL   Source root. Default: NOAA NODD CORe monthly archive

Reruns skip archives already verified as complete. Partial files are resumed.
EOF
}

if (($#)) && [[ "${1:-}" =~ ^(-h|--help)$ ]]; then
  usage
  exit 0
fi

if (($#)); then
  KINDS=("$@")
else
  KINDS=(flx pgb)
fi

for kind in "${KINDS[@]}"; do
  case "$kind" in
    flx|pgb) ;;
    *)
      printf 'Unsupported kind: %s\n' "$kind" >&2
      usage >&2
      exit 2
      ;;
  esac
done

remote_size() {
  curl -fsSI --retry 5 --retry-delay 5 --retry-connrefused "$1" |
    awk 'BEGIN{IGNORECASE=1} /^content-length:/ {gsub("\r", "", $2); print $2; exit}'
}

archive_complete() {
  local outfile=$1
  local url=$2
  local marker=$3
  local local_size expected_size

  if [[ -f "$marker" && -s "$outfile" ]]; then
    return 0
  fi

  if [[ ! -f "$outfile" ]]; then
    return 1
  fi

  local_size=$(wc -c < "$outfile")
  expected_size=$(remote_size "$url") || return 1

  if [[ "$local_size" -ne "$expected_size" ]]; then
    return 1
  fi

  touch "$marker"
  return 0
}

for month in {01..12}; do
  for kind in "${KINDS[@]}"; do
    outdir="$OUT_ROOT/month/$kind/$YEAR"
    outfile="$outdir/$kind.$YEAR$month.grb"
    url="$BASE_URL/$kind/$YEAR/$kind.$YEAR$month.grb"
    marker="$outdir/.$kind.$YEAR$month.complete"

    mkdir -p "$outdir"

    if archive_complete "$outfile" "$url" "$marker"; then
      printf 'Skipping complete %s archive %s-%s\n' "$kind" "$YEAR" "$month"
      continue
    fi

    printf 'Downloading %s\n' "$url"
    curl -fL -C - --retry 5 --retry-delay 5 --retry-connrefused \
      -o "$outfile" "$url"
    touch "$marker"
  done
done
