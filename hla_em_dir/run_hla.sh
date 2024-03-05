
echo "here"
DEFAULT_VOLUME_PATH="/data"
VOLUME_PATH=${1:-$DEFAULT_VOLUME_PATH}

echo "running ${VOLUME_PATH}/HLA_EM.py"
echo ""

pip3 install matplotlib
pip3 install whichcraft

python3 "${VOLUME_PATH}/HLA_EM.py" -t 4 -o "${VOLUME_PATH}/output" -s "${VOLUME_PATH}/EnsembleGenome_STAR_without_scaffolds" -r "${VOLUME_PATH}/hla_gen.fasta" "${VOLUME_PATH}/reference/samples/trial_0/sim.HLA.reads_01.fq" "${VOLUME_PATH}/reference/samples/trial_0/sim.HLA.reads_02.fq"
