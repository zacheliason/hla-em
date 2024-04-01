#VOLUME_PATH="/data"
#VOLUME_PATH="/Users/zeliason/Downloads/hla-em"

HOST_PATH=$(pwd)
VOLUME_PATH=/data

# Install required Python packages
#pip3 install matplotlib
#pip3 install whichcraft
#pip3 install matplotlib==3.4.3

for trial in {0..47}; do
    TRIAL_DIR="${VOLUME_PATH}/reference/samples/trial_${trial}"

    echo ""
    echo "now working on ${TRIAL_DIR}"
    OUTPUT_DIR="${VOLUME_PATH}/output_training/trial_${trial}"

    echo "$HOST_PATH":"$VOLUME_PATH"

#    docker run -it --rm -v "$HOST_PATH":"$VOLUME_PATH" zeliason/hla-em:latest -t 4 -o "${VOLUME_PATH}/output" -s "${VOLUME_PATH}/hla_gen.fasta_STAR" -r "${VOLUME_PATH}/hla_gen.fasta" "${VOLUME_PATH}/reference/samples/trial_0/sim.HLA.reads_01.fq" "${VOLUME_PATH}/reference/samples/trial_0/sim.HLA.reads_02.fq"
    docker run -it --rm -v "$HOST_PATH":"$VOLUME_PATH" zeliason/hla-em:latest
#    docker run -it --rm -v "$HOST_PATH":"$VOLUME_PATH" zeliason/hla-em:latest -t 4 -o "${VOLUME_PATH}/output" -s "${VOLUME_PATH}/EnsembleGenome_STAR_without_scaffolds" "${VOLUME_PATH}/reference/samples/trial_0/sim.HLA.reads_01.fq" "${VOLUME_PATH}/reference/samples/trial_0/sim.HLA.reads_02.fq"

done
