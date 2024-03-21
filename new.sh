VOLUME_PATH="/data"

# Install required Python packages
pip3 install matplotlib
pip3 install whichcraft

outFilterScoreMinOverLreads_values=(0.66 0.7 0.75)
outFilterMatchNminOverLreads_values=(0.66 0.7 0.75)
outFilterMultimapNmaxs_values=(50 60 70)
winAnchorMultimapNmaxs_values=(50 60 70)

array_length=${#outFilterScoreMinOverLreads_values[@]}

# Define the trials to iterate over
for trial in {0..4}; do
    TRIAL_DIR="${VOLUME_PATH}/src/reference/samples/trial_${trial}"

    echo ""
    echo "now working on ${TRIAL_DIR}"

    # Loop over the indices of one of the arrays
#    for ((i = 0; i < array_length; i++)); do
#        outFilterScoreMinOverLreads=${outFilterScoreMinOverLreads_values[i]}
#        outFilterMatchNminOverLreads=${outFilterMatchNminOverLreads_values[i]}
#        outFilterMultimapNmaxs=${outFilterMultimapNmaxs_values[i]}
#        winAnchorMultimapNmaxs=${winAnchorMultimapNmaxs_values[i]}
#
#        echo "Parameter set $((i + 1))"
#        echo "outFilterScoreMinOverLreads: $outFilterScoreMinOverLreads"
#        echo "outFilterMatchNminOverLreads: $outFilterMatchNminOverLreads"
#        echo "outFilterMultimapNmaxs: $outFilterMultimapNmaxs"
#        echo "winAnchorMultimapNmaxs: $winAnchorMultimapNmaxs"

#        OUTPUT_DIR="${VOLUME_PATH}/src/output/trial_${trial}_$((i+1))"
        OUTPUT_DIR="${VOLUME_PATH}/src/output/trial_${trial}"

        python3 "${VOLUME_PATH}/HLA_EM.py" \
            -t 4 \
            -o "${OUTPUT_DIR}" \
            -s "${VOLUME_PATH}/src/EnsembleGenome_STAR_without_scaffolds" \
            --starHLA "${VOLUME_PATH}/src/hla_gen.fasta_STAR" \
            --shortcut \
            -r "${VOLUME_PATH}/src/hla_gen.fasta" \
            "${TRIAL_DIR}/sim.HLA.reads_01.fq"
#            --outFilterScoreMinOverLreads "$outFilterScoreMinOverLreads" \
#            --outFilterMatchNminOverLreads "$outFilterMatchNminOverLreads" \
#            --outFilterMultimapNmaxs "$outFilterMultimapNmaxs" \
#            --winAnchorMultimapNmaxs "$winAnchorMultimapNmaxs" \
#            "${TRIAL_DIR}/sim.HLA.reads_01.fq" #\
#            "${TRIAL_DIR}/sim.HLA.reads_02.fq"
#    done
done
