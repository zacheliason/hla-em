#!/bin/bash

##TRIAL_NUM=${1}
#if [ -z "$1" ]; then
#    echo "Error: TRIAL_NUM not provided. Please provide TRIAL_NUM as a command-line argument."
#    exit 1
#fi
#
##TRIAL_NUM=$1
#echo $TRIAL_NUM

for TRIAL_NUM in {0..31}; do
  echo "working on ${TRIAL_NUM}"
  VOLUME_PATH="/out"
  HOST_PATH=$(pwd)
#  echo $HOST_PATH
#  echo $VOLUME_PATH
#
#  echo "mounting $HOST_PATH":"$VOLUME_PATH"

  # Run the Docker command with the specified volume path
#  docker run --rm -v "$HOST_PATH":"$VOLUME_PATH" -t fred2/optitype -i "$VOLUME_PATH"/reference_paired/samples/trial_"$TRIAL_NUM"/out.Unmapped.out.mate1.fastq --verbose -r -o "$VOLUME_PATH"/optitype_output_single/optitype_results_run_"$TRIAL_NUM"
  docker run --rm -v "$HOST_PATH":"$VOLUME_PATH" -t fred2/optitype -i "$VOLUME_PATH"/reference_paired/samples/trial_"$TRIAL_NUM"/out.Unmapped.out.mate1.fastq "$VOLUME_PATH"/reference_paired/samples/trial_"$TRIAL_NUM"/out.Unmapped.out.mate2.fastq --verbose -r -o "$VOLUME_PATH"/optitype_output_paired/optitype_results_run_"$TRIAL_NUM"
done