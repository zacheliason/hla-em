#!/bin/bash

for TRIAL_NUM in {0..24}; do
  echo ""
  echo "working on ${TRIAL_NUM}"
  VOLUME_PATH="/out"
  HOST_PATH=$(pwd)

  TRIAL_DIR = "$VOLUME_PATH"/reference/samples/trial_"$TRIAL_NUM"
  OUTPUT_DIR = "$VOLUME_PATH"/optitype_output_paired/optitype_results_run_"$TRIAL_NUM"

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
