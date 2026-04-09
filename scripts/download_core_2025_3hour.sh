#!/usr/bin/env bash
set -euo pipefail

BASE_URL="${BASE_URL:-https://storage.googleapis.com/noaa-nws-ncep-core/grib/3hour}"
OUT_ROOT="${OUT_ROOT:-$HOME/Data/CORe}"
START_TS=2025010100
END_TS=2026010100

usage() {
  cat <<'EOF'
Usage: download_core_2025_3hour.sh [flx] [pgb]

Download all 3-hourly 2025 CORe files into:
  $OUT_ROOT/{flx,pgb}/YYYY/MM/{flx,pgb}.YYYYMMDDHH.grb

Environment:
  OUT_ROOT   Destination root. Default: $HOME/Data/CORe
  BASE_URL   Source root. Default: NOAA NODD CORe 3-hour archive

Reruns skip months already verified as complete. Incomplete months are resumed.
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

next_ts() {
  date -u -d "${1:0:4}-${1:4:2}-${1:6:2} ${1:8:2}:00 UTC +3 hours" +%Y%m%d%H
}

next_month_ts() {
  date -u -d "${1}-${2}-01 +1 month" +%Y%m0100
}

expected_files_in_month() {
  local ndays
  ndays=$(date -u -d "${1}-${2}-01 +1 month -1 day" +%d)
  printf '%d\n' "$((10#$ndays * 8))"
}

last_ts_in_month() {
  date -u -d "${1}-${2}-01 +1 month -3 hours" +%Y%m%d%H
}

remote_size() {
  curl -fsSI --retry 5 --retry-delay 5 --retry-connrefused "$1" |
    awk 'BEGIN{IGNORECASE=1} /^content-length:/ {gsub("\r", "", $2); print $2; exit}'
}

month_complete() {
  local kind=$1
  local year=$2
  local month=$3
  local outdir="$OUT_ROOT/$kind/$year/$month"
  local marker="$outdir/.complete"
  local actual expected last_ts last_file last_url local_size expected_size

  if [[ -f "$marker" ]]; then
    return 0
  fi

  if [[ ! -d "$outdir" ]]; then
    return 1
  fi

  expected=$(expected_files_in_month "$year" "$month")
  actual=$(find "$outdir" -maxdepth 1 -type f -name "$kind.$year$month*.grb" | wc -l)
  if [[ "$actual" -ne "$expected" ]]; then
    return 1
  fi

  last_ts=$(last_ts_in_month "$year" "$month")
  last_file="$outdir/$kind.$last_ts.grb"
  last_url="$BASE_URL/$kind/$year/$month/$kind.$last_ts.grb"

  if [[ ! -f "$last_file" ]]; then
    return 1
  fi

  local_size=$(wc -c < "$last_file")
  expected_size=$(remote_size "$last_url") || return 1

  if [[ "$local_size" -ne "$expected_size" ]]; then
    return 1
  fi

  touch "$marker"
  return 0
}

month_ts=$START_TS
while [[ "$month_ts" < "$END_TS" ]]; do
  year=${month_ts:0:4}
  month=${month_ts:4:2}
  next_month=$(next_month_ts "$year" "$month")
  pending=()

  for kind in "${KINDS[@]}"; do
    if month_complete "$kind" "$year" "$month"; then
      printf 'Skipping complete %s month %s-%s\n' "$kind" "$year" "$month"
    else
      pending+=("$kind")
    fi
  done

  if ((${#pending[@]} == 0)); then
    month_ts=$next_month
    continue
  fi

  ts=$month_ts
  while [[ "$ts" < "$next_month" ]]; do
    for kind in "${pending[@]}"; do
      outdir="$OUT_ROOT/$kind/$year/$month"
      outfile="$outdir/$kind.$ts.grb"
      url="$BASE_URL/$kind/$year/$month/$kind.$ts.grb"

      mkdir -p "$outdir"
      printf 'Downloading %s\n' "$url"
      curl -fL -C - --retry 5 --retry-delay 5 --retry-connrefused \
        -o "$outfile" "$url"
    done

    ts=$(next_ts "$ts")
  done

  for kind in "${pending[@]}"; do
    touch "$OUT_ROOT/$kind/$year/$month/.complete"
  done

  month_ts=$next_month
done
