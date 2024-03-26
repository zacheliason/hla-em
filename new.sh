VOLUME_PATH="/data"

# Install required Python packages
pip3 install matplotlib
pip3 install whichcraft

for trial in {10..19}; do
    TRIAL_DIR="${VOLUME_PATH}/reference/samples/trial_${trial}"

    echo ""
    echo "now working on ${TRIAL_DIR}"
    OUTPUT_DIR="${VOLUME_PATH}/output/trial_${trial}"

    python3 "${VOLUME_PATH}/HLA_EM.py" \
        -t 4 \
        --shortcut \
        -o "${OUTPUT_DIR}" \
        -s "${VOLUME_PATH}/EnsembleGenome_STAR_without_scaffolds" \
        -r "${VOLUME_PATH}/hla_gen.fasta" \
        "${TRIAL_DIR}/sim.HLA.reads_01.fq"
#            --shortcut \
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
