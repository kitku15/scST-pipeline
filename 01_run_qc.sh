#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=04:00:00
#PBS -N run_qc
#PBS -J 1-4
#PBS -j oe

set -euo pipefail
cd "$PBS_O_WORKDIR"

CONFIG_NAME="configs/config_pCosMx.toml"
SIF_IMAGE="$(readlink -f kitku.sif)"

[[ -f "$SIF_IMAGE" ]] || { echo "Missing SIF image: $SIF_IMAGE" >&2; exit 1; }
[[ -f "$CONFIG_NAME" ]] || { echo "Missing config: $CONFIG_NAME" >&2; exit 1; }

# get base directory from the TOML file
BASE_DIR=$(grep -oP 'base_raw_dir\s*=\s*"\K[^"]+' "$CONFIG_NAME")

# Count how many dataset directories exist inside the base directory
NUM_SLIDES=$(find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)

# If this job's array index is greater than the slides we actually have, exit immediately
if [ "$PBS_ARRAY_INDEX" -gt "$NUM_SLIDES" ]; then
    echo "Skipping... Array index $PBS_ARRAY_INDEX is greater than the $NUM_SLIDES slide(s) provided."
    exit 0
fi

#########################################
## Output some useful job information. ##
#########################################

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: job identifier is $PBS_JOBID
echo PBS: array index is $PBS_ARRAY_INDEX
echo ------------------------------------------------------

APPTAINER_TMPDIR_DEFAULT="$PBS_O_WORKDIR/tmp/apptainer-tmp-${PBS_ARRAY_INDEX}"
APPTAINER_WORKDIR_DEFAULT="$PBS_O_WORKDIR/tmp/apptainer-work-${PBS_ARRAY_INDEX}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$APPTAINER_TMPDIR_DEFAULT}"
export APPTAINER_WORKDIR="${APPTAINER_WORKDIR:-$APPTAINER_WORKDIR_DEFAULT}"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_WORKDIR"

KEEP_APPTAINER_TMP="${KEEP_APPTAINER_TMP:-0}"
cleanup_apptainer_tmp() {
  [[ "$KEEP_APPTAINER_TMP" == "1" ]] && return 0
  local base="$(readlink -f "$PBS_O_WORKDIR/tmp" 2>/dev/null || true)"
  [[ -n "$base" && -d "$base" ]] || return 0
  local tmpdir="$(readlink -f "$APPTAINER_TMPDIR" 2>/dev/null || true)"
  local workdir="$(readlink -f "$APPTAINER_WORKDIR" 2>/dev/null || true)"
  case "$tmpdir" in "$base"/*) rm -rf --one-file-system "$tmpdir" || true ;; esac
  case "$workdir" in "$base"/*) rm -rf --one-file-system "$workdir" || true ;; esac
}
trap cleanup_apptainer_tmp EXIT
echo "Running Spatial Pipeline: Format & QC for Slide ${PBS_ARRAY_INDEX}" >&2

apptainer run --writable-tmpfs -W "$APPTAINER_WORKDIR" \
  --bind "$PBS_O_WORKDIR:/app" \
  "$SIF_IMAGE" \
  "$CONFIG_NAME" --modules 0 1 --sample_index "$PBS_ARRAY_INDEX"

echo "QC finished for Slide ${PBS_ARRAY_INDEX}!" >&2