#!/bin/bash
#PBS -l select=1:ncpus=32:mem=512gb
#PBS -l walltime=48:00:00
#PBS -N run_ds
#PBS -j oe

set -euo pipefail
cd "$PBS_O_WORKDIR"

#########################################
## Output some useful job information. ##
#########################################

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: job identifier is $PBS_JOBID
echo ------------------------------------------------------

APPTAINER_TMPDIR_DEFAULT="$PBS_O_WORKDIR/tmp/apptainer-tmp"
APPTAINER_WORKDIR_DEFAULT="$PBS_O_WORKDIR/tmp/apptainer-work"
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

# USE THE MASTER DOWNSTREAM CONFIG
CONFIG_NAME="config_tyler_downstream_manual.toml"
SIF_IMAGE="$(readlink -f kitku.sif)"

[[ -f "$SIF_IMAGE" ]] || { echo "Missing SIF image: $SIF_IMAGE" >&2; exit 1; }
[[ -f "$CONFIG_NAME" ]] || { echo "Missing config: $CONFIG_NAME" >&2; exit 1; }

echo "Running Downstream Spatial Pipeline: Merge -> Analysis" >&2

apptainer run --writable-tmpfs -W "$APPTAINER_WORKDIR" \
  --bind "$PBS_O_WORKDIR:/app" \
  "$SIF_IMAGE" \
  "$CONFIG_NAME" --modules 1b 2 3 4 5 6 7 8 8b 8c 9 10
echo "Pipeline finished!" >&2 