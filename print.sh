VOLUME_PATH="/Users/zeliason/Desktop/hla-em"

for trial in {0..19}; do
    TRIAL_DIR="${VOLUME_PATH}/reference/samples/trial_${trial}"

    echo ""
    OUTPUT_DIR="${VOLUME_PATH}/output_paired/trial_${trial}"
    echo "-t 4 --shortcut -o ${OUTPUT_DIR} -s ${VOLUME_PATH}/EnsembleGenome_STAR_without_scaffolds -r ${VOLUME_PATH}/hla_gen.fasta ${TRIAL_DIR}/sim.HLA.reads_01.fq ${TRIAL_DIR}/sim.HLA.reads_02.fq"
#"${TRIAL_DIR}/sim.HLA.reads_02.fq"
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
