#!/bin/bash

# Get the parent directory of this script
SCRIPT_PATH=$(readlink -f "$0")
PARENT_DIR=$(dirname "$SCRIPT_PATH")

# Use parent directory to mount docker volume
HOST_PATH=$PARENT_DIR
VOLUME_PATH="/out"

for TRIAL_NUM in {0..24}; do
  echo ""
  echo "working on ${TRIAL_NUM}"

  TRIAL_DIR="$VOLUME_PATH"/reference/samples/trial_"$TRIAL_NUM"
  OUTPUT_DIR="$VOLUME_PATH"/optitype_output_paired/optitype_results_run_"$TRIAL_NUM"

  docker run --rm \
    -v "$HOST_PATH":"$VOLUME_PATH" \
    -t \
    fred2/optitype \
    -i "$TRIAL_DIR"/out.Unmapped.out.mate1.fastq \
       "$TRIAL_DIR"/out.Unmapped.out.mate2.fastq \
    --verbose \
    -r \
    -o "$OUTPUT_DIR"/optitype_results_run_"$TRIAL_NUM"

done
