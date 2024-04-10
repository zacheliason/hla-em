#!/bin/bash

ALLELE_RECORD_PATH="$(pwd)/reference/allele_record.csv"

# Evaluate HLA-EM with paired reads
OUTPUT_PATH="$(pwd)/output"
NAME="HLA-EM_paired"

CSV_PATH="$(pwd)"/"$NAME"_scores.csv

python3 src/evaluate_trials.py -n "$NAME" -o "$OUTPUT_PATH" -r "$ALLELE_RECORD_PATH" --outname "$CSV_PATH"

# Evaluate Optitype with paired reads
OUTPUT_PATH="$(pwd)/optitype_output_paired"
NAME="optitype_paired"

CSV_PATH="$(pwd)"/"$NAME"_scores.csv

python3 src/evaluate_trials.py -n "$NAME" -o "$OUTPUT_PATH" -r "$ALLELE_RECORD_PATH" --outname "$CSV_PATH"

# Evaluate Optitype with single reads
OUTPUT_PATH="$(pwd)/optitype_output_single"
NAME="optitype_single"

CSV_PATH="$(pwd)"/"$NAME"_scores.csv

python3 src/evaluate_trials.py -n "$NAME" -o "$OUTPUT_PATH" -r "$ALLELE_RECORD_PATH" --outname "$CSV_PATH"
