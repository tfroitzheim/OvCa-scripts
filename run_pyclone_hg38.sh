#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="/wrk/natgen_reviews/sarcoma_vcfs/workspace_hg38/pyclone_input"
OUTPUT_DIR="/wrk/natgen_reviews/sarcoma_vcfs/workspace_hg38/pyclone_output"
mkdir -p "$OUTPUT_DIR"

for src in src16 src20 src21 src23 src26 src27 src32 src35 src39; do
    infile="$INPUT_DIR/${src}.tsv"
    [[ -f "$infile" ]] || { echo "!! missing $infile, skipping $src" >&2; continue; }

    pyclone-vi fit \
        -i "$infile" \
        -o "$OUTPUT_DIR/${src}.h5" \
        -d beta-binomial -r 10 -t 4 --seed 42

    pyclone-vi write-results-file \
        -i "$OUTPUT_DIR/${src}.h5" \
        -o "$OUTPUT_DIR/${src}_results.tsv"

    echo "$src: done"
done
