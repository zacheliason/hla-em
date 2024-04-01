VOLUME_PATH="/data"
#VOLUME_PATH="/Users/zeliason/Desktop/hla-em"

# Install required Python packages
#pip3 install matplotlib
#pip3 install whichcraft
#pip3 install matplotlib==3.4.3

for trial in {0..47}; do
    TRIAL_DIR="${VOLUME_PATH}/reference/samples/trial_${trial}"

    echo ""
    echo "now working on ${TRIAL_DIR}"
    OUTPUT_DIR="${VOLUME_PATH}/output_training/trial_${trial}"

    python3 "${VOLUME_PATH}/HLA_EM.py" \
        -t 4 \
        --shortcut \
        --suppress_figs \
        --training ${VOLUME_PATH}/output_training/training.tsv \
        -o "${OUTPUT_DIR}" \
        -s "${VOLUME_PATH}/EnsembleGenome_STAR_without_scaffolds" \
        -r "${VOLUME_PATH}/hla_gen.fasta" \
        "${TRIAL_DIR}/sim.HLA.reads_01.fq" \
        "${TRIAL_DIR}/sim.HLA.reads_02.fq"
#            "${TRIAL_DIR}/sim.HLA.reads_02.fq"
#            --starHLA "${VOLUME_PATH}/hla_gen.fasta_STAR" \
#            "${TRIAL_DIR}/sim.HLA.reads_02.fq"

##    OUTPUT_DIR="${VOLUME_PATH}/output/trial_${trial}_dd"
#    OUTPUT_DIR_OLD="${VOLUME_PATH}/output/trial_${trial}"
#
#    mkdir $OUTPUT_DIR
#    cp "$OUTPUT_DIR_OLD/trial_${trial}.Log.final.out" "$OUTPUT_DIR/trial_${trial}_dd.Log.final.out"
#
#    python3 "${VOLUME_PATH}/HLA_EM.py" \
#        -t 4 \
#        -d \
#        --shortcut \
#        -o "${OUTPUT_DIR}" \
#        -s "${VOLUME_PATH}/EnsembleGenome_STAR_without_scaffolds" \
#        -r "${VOLUME_PATH}/hla_gen.fasta" \
#        "${TRIAL_DIR}/sim.HLA.reads_01.fq"
done
