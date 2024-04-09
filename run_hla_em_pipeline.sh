#!/bin/bash

# Get the parent directory of this script
SCRIPT_PATH=$(readlink -f "$0")
PARENT_DIR=$(dirname "$SCRIPT_PATH")

# Use parent directory to mount docker volume
HOST_PATH=$PARENT_DIR
VOLUME_PATH=/data

# Change number range as needed for trials in reference directory
for trial in {0..24}; do
    TRIAL_DIR="${VOLUME_PATH}/reference/samples/trial_${trial}"

    # Change output directory as needed
    OUTPUT_DIR="${VOLUME_PATH}/output/trial_${trial}"

    echo ""
    echo "now working on ${TRIAL_DIR}"

    # Disable dust for now, I haven't returned dust filtering to the pipeline yet
    docker run -it --rm -v "$HOST_PATH":"$VOLUME_PATH" zeliason/hla-em:latest \
        -t 4 \
        -d \
        --suppress_figs \
        -o "${OUTPUT_DIR}" \
        -s "${VOLUME_PATH}/EnsembleGenome_STAR_without_scaffolds" \
        -r "${VOLUME_PATH}/hla_gen.fasta" \
        "${TRIAL_DIR}/sim.HLA.reads_01.fq" \
        "${TRIAL_DIR}/sim.HLA.reads_02.fq"

done
