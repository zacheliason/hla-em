#!/bin/bash

ALLELE_RECORD_PATH="$(pwd)/reference/allele_record.csv"

OUTPUT_PATH="$(pwd)/output"
NAME="HLA-EM_paired"

CSV_PATH="$(pwd)"/"$NAME"_scores.csv

python3 src/evaluate_trials.py -n "$NAME" -o "$OUTPUT_PATH" -r "$ALLELE_RECORD_PATH" --outname "$CSV_PATH"

OUTPUT_PATH="$(pwd)/optitype_output_paired"
NAME="optitype_single"

CSV_PATH="$(pwd)"/"$NAME"_scores.csv

python3 src/evaluate_trials.py -n "$NAME" -o "$OUTPUT_PATH" -r "$ALLELE_RECORD_PATH" --outname "$CSV_PATH"
