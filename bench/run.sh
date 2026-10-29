#!/usr/bin/env bash
# Run one (configuration, dataset) pair through run_signals.R, store the output
# under a key that identifies the exact code, and score it.
#
#   bench/run.sh <config> <dataset> [--image IMG] [--script PATH] [--dry-run]
#
#   <config>   a name from bench/configs.tsv
#   <dataset>  a name from bench/datasets.tsv
#   --image    container image (default signals:dev, else signals:latest)
#   --script   run_signals.R to use (default signals_isabl/run_signals.R)
#   --dry-run  print what would run, and the provenance, without running
#
# Output: bench/results/<dataset>/<config>/<signals-sha>_<script-sha>/
#   hscn.csv.gz, qc.csv, confidence.csv, heatmap.png, ...   the run
#   meta.tsv      flags, image, both shas, wall time, exit code
#   metrics.csv   bench/R/score.R on this run
#   run.log       stdout/stderr
#
# The key includes BOTH shas because the same flags have meant different things
# at different script versions - see the r14 note in history.tsv.
set -euo pipefail

REPO="${SIGNALS_REPO:-$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)}"
B="$REPO/bench"

[ $# -ge 2 ] || { sed -n '2,20p' "$0"; exit 2; }
CONFIG="$1"; DATASET="$2"; shift 2
IMAGE=""; SCRIPT="$REPO/signals_isabl/run_signals.R"; DRY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --image)   IMAGE="$2"; shift 2 ;;
    --script)  SCRIPT="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

# --- look up the configuration and dataset ----------------------------------
tsv_row() { awk -F'\t' -v k="$2" '!/^#/ && $1==k' "$1"; }
CROW="$(tsv_row "$B/configs.tsv" "$CONFIG")"
[ -n "$CROW" ] || { echo "no configuration '$CONFIG' in configs.tsv" >&2; exit 2; }
FLAGS="$(printf '%s' "$CROW" | cut -f2)"
[ -f "$B/datasets.tsv" ] || { echo "no bench/datasets.tsv - copy datasets.example.tsv and fill it in" >&2; exit 2; }
DROW="$(tsv_row "$B/datasets.tsv" "$DATASET")"
[ -n "$DROW" ] || { echo "no dataset '$DATASET' in datasets.tsv" >&2; exit 2; }
IFS=$'\t' read -r _ DDIR QC READS COUNTS SEX VERDICT <<< "$DROW"
DDIR="$REPO/$DDIR"
for f in "$QC" "$READS" "$COUNTS"; do
  [ -f "$DDIR/$f" ] || { echo "missing input: $DDIR/$f" >&2; exit 1; }
done
[ -f "$SCRIPT" ] || { echo "missing script: $SCRIPT" >&2; exit 1; }

if [ -z "$IMAGE" ]; then
  if docker image inspect signals:dev >/dev/null 2>&1; then IMAGE=signals:dev
  else IMAGE=signals:latest; fi
fi

# --- provenance ---------------------------------------------------------------
SIG_SHA="$(git -C "$REPO" rev-parse --short HEAD)"
git -C "$REPO" diff --quiet HEAD -- R src DESCRIPTION NAMESPACE 2>/dev/null || SIG_SHA="${SIG_SHA}-dirty"
SCRIPT_SHA="$(shasum -a 256 "$SCRIPT" | cut -c1-8)"
IMAGE_ID="$(docker image inspect --format '{{.Id}}' "$IMAGE" 2>/dev/null || true)"
IMAGE_ID="${IMAGE_ID:7:12}"; IMAGE_ID="${IMAGE_ID:-missing}"
OUT="$B/results/$DATASET/$CONFIG/${SIG_SHA}_${SCRIPT_SHA}"

SHARED="--seed 42 --viterbiver cpp --mincells 2 --ncores 5 --sex $SEX"
CMD=(docker run --rm --platform linux/x86_64 -w "$OUT" -v "$REPO:$REPO" "$IMAGE"
     Rscript "$SCRIPT" $SHARED $FLAGS
     --hmmcopyqc "$DDIR/$QC" --hmmcopyreads "$DDIR/$READS" --allelecounts "$DDIR/$COUNTS"
     --rephase_report "$OUT/confidence.csv"
     --qcplot "$OUT/qcplot.png" --heatmap "$OUT/heatmap.png" --heatmapraw "$OUT/heatmapraw.png"
     --csvfile "$OUT/hscn.csv.gz" --qccsvfile "$OUT/qc.csv" --Rdatafile "$OUT/signals.Rdata")

echo "config    $CONFIG   [$FLAGS]"
echo "dataset   $DATASET  ($DDIR)"
echo "image     $IMAGE  ($IMAGE_ID)"
echo "signals   $SIG_SHA    script $SCRIPT_SHA  ($SCRIPT)"
echo "output    $OUT"
if [ "$DRY" = 1 ]; then printf '\n'; printf '%q ' "${CMD[@]}"; printf '\n'; exit 0; fi

[ "$IMAGE_ID" != "missing" ] || { echo "image $IMAGE not found - build it first (dev/build-dev-image.sh)" >&2; exit 1; }
if [ -f "$OUT/hscn.csv.gz" ]; then
  echo "already run for this code version - delete $OUT to rerun" >&2; exit 0
fi
mkdir -p "$OUT"

# --- run ------------------------------------------------------------------------
T0=$(date +%s); set +e
"${CMD[@]}" > "$OUT/run.log" 2>&1; RC=$?
set -e; T1=$(date +%s)

{
  printf 'key\tvalue\n'
  printf 'config\t%s\nflags\t%s\nshared_flags\t%s\n' "$CONFIG" "$FLAGS" "$SHARED"
  printf 'dataset\t%s\nimage\t%s\nimage_id\t%s\n' "$DATASET" "$IMAGE" "$IMAGE_ID"
  printf 'signals_sha\t%s\nscript\t%s\nscript_sha256_8\t%s\n' "$SIG_SHA" "$SCRIPT" "$SCRIPT_SHA"
  printf 'started\t%s\nwall_seconds\t%s\nexit_code\t%s\n' "$(date -r "$T0" '+%F %T')" "$((T1 - T0))" "$RC"
} > "$OUT/meta.tsv"

echo "exit $RC after $((T1 - T0))s"
[ "$RC" = 0 ] || { echo "run failed - see $OUT/run.log" >&2; exit "$RC"; }

if [ -n "$VERDICT" ] && [ -f "$REPO/$VERDICT" ]; then
  Rscript "$B/R/score.R" --verdict "$REPO/$VERDICT" --input "$DDIR/$READS" \
    --out "$OUT/metrics.csv" "$CONFIG=$OUT/hscn.csv.gz"
fi
