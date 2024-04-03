#!/bin/bash

HOST_PATH=$(pwd)
VOLUME_PATH=/data

echo "$HOST_PATH":"$VOLUME_PATH"

docker run -it --rm -v "$HOST_PATH":"$VOLUME_PATH" zeliason/hla-em:latest -t 4 -d -o "${VOLUME_PATH}/output" -s "${VOLUME_PATH}/hla_gen.fasta_STAR" -r "${VOLUME_PATH}/hla_gen.fasta" "${VOLUME_PATH}/reference/samples/trial_0/sim.HLA.reads_01.fq" "${VOLUME_PATH}/reference/samples/trial_0/sim.HLA.reads_02.fq"
