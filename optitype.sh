#!/bin/bash

##TRIAL_NUM=${1}
#if [ -z "$1" ]; then
#    echo "Error: TRIAL_NUM not provided. Please provide TRIAL_NUM as a command-line argument."
#    exit 1
#fi
#
##TRIAL_NUM=$1
#echo $TRIAL_NUM

for TRIAL_NUM in {1..1}; do
  echo "working on ${TRIAL_NUM}"
  DEFAULT_VOLUME_PATH="/out"
  VOLUME_PATH=${2:-$DEFAULT_VOLUME_PATH}
  HOST_PATH=$(pwd)
  echo $HOST_PATH
  echo $VOLUME_PATH

  echo "mounting $HOST_PATH":"$VOLUME_PATH"

  # Run the Docker command with the specified volume path
  #docker run --rm -v "$HOST_PATH":"$VOLUME_PATH" -t fred2/optitype -i "$VOLUME_PATH"/reference/samples/trial_"$TRIAL_NUM"/unmapped_reads1.fastq "$VOLUME_PATH"/reference/samples/trial_"$TRIAL_NUM"/unmapped_reads2.fastq  --verbose -r -o "$VOLUME_PATH"/optitype_results_run_"$TRIAL_NUM"
  docker run --rm -v "$HOST_PATH":"$VOLUME_PATH" -t fred2/optitype -i "$VOLUME_PATH"/reference/samples/trial_"$TRIAL_NUM"/out.Unmapped.out.mate1.fastq  --verbose -r -o "$VOLUME_PATH"/output/optitype_results_run_"$TRIAL_NUM"
done