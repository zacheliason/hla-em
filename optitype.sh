#!/bin/bash

TRIAL_NUM=${1}
if [ -z "$1" ]; then
    echo "Error: TRIAL_NUM not provided. Please provide TRIAL_NUM as a command-line argument."
    exit 1
fi

TRIAL_NUM=$1
DEFAULT_VOLUME_PATH="/src"
VOLUME_PATH=${2:-$DEFAULT_VOLUME_PATH}
HOST_PATH=$(pwd)/src

echo "mounting $HOST_PATH":"$VOLUME_PATH"

# Run the Docker command with the specified volume path
docker run --rm -v "$HOST_PATH":"$VOLUME_PATH" -t fred2/optitype -i "$VOLUME_PATH"/reference/samples/trial_"$TRIAL_NUM"/unmapped_reads1.fastq "$VOLUME_PATH"/reference/samples/trial_"$TRIAL_NUM"/unmapped_reads2.fastq  --verbose -r -o "$VOLUME_PATH"/optitype_results_run_"$TRIAL_NUM"

