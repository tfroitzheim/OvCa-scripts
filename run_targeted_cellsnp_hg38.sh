#!/usr/bin/env bash
set -euo pipefail
### this is run on cancergenome, not meg 
DNA_ALIGNED="/wrk/data/human_genomic/aligned"
RNA_ALIGNED="/wrk/data/human/aligned"
CANDIDATE_VCF_DIR="/wrk/data/tim/natgen_revisions_cg/cellsnp_sarcomas_bulkWGS/candidate_snps_hg38"
OUTPUT_BASE="/wrk/data/tim/natgen_revisions_cg/cellsnp_sarcomas_bulkWGS/output_hg38"
THREADS_PER_CELL=2
PARALLEL_CELLS=6
MIN_COUNT=1
MIN_MAF=0

declare -A SRC_PLATES=(
  [src16]="LBA00702 LBA00703"
  [src20]="LBA01102"
  [src21]="LBA01202"
  [src23]="LBA01402 LBA01403"
  [src26]="LBA01701 LBA01702"
  [src27]="LBA01601"
  [src32]="LBA02102"
  [src35]="LBA02501 LBA02502"
  [src39]="LBA02602 LBA02603"
)

find_dna_cells() {
  local plate="$1"
  find "$DNA_ALIGNED" -maxdepth 1 -type d -iregex ".*/${plate}[Dd][-_.]?.*" 2>/dev/null | sort
}

find_rna_cells() {
  local plate="$1"
  find "$RNA_ALIGNED" -maxdepth 1 -type d -iregex ".*/${plate}[Rr][-_.]?.*" 2>/dev/null | sort
}

run_one_cell() {
  local cell_dir="$1" candidate_vcf="$2" cohort_out="$3" modality="$4"
  local cell=$(basename "$cell_dir")
  local bam

  if [[ "$modality" == "dna" ]]; then
    bam="$cell_dir/${cell}.dedup.bam"
  else
    bam="$cell_dir/${cell}.sorted.bam"
  fi

  if [[ ! -f "$bam" || ! -f "${bam}.bai" ]]; then
    echo "[MISSING] $cell ($modality): no bam/.bai at $bam" >&2
    return 0
  fi

  local cell_out="$cohort_out/$cell"
  if [[ -s "$cell_out/cellSNP.base.vcf.gz" ]]; then
    echo "[SKIP] $cell already done"
    return 0
  fi
  mkdir -p "$cell_out"

  cellsnp-lite \
    -s "$bam" \
    -R "$candidate_vcf" \
    -O "$cell_out" \
    -p "$THREADS_PER_CELL" \
    --minCOUNT "$MIN_COUNT" \
    --minMAF "$MIN_MAF" \
    --UMItag None \
    --cellTAG None \
    --gzip \
    > "$cell_out/cellsnp.log" 2>&1
  echo "[DONE] $cell ($modality)"
}
export -f run_one_cell
export THREADS_PER_CELL MIN_COUNT MIN_MAF

for src_id in "${!SRC_PLATES[@]}"; do
  candidate_vcf="$CANDIDATE_VCF_DIR/${src_id}_hg38.vcf.gz"
  if [[ ! -f "$candidate_vcf" ]]; then
    echo "!! No hg38 candidate VCF for $src_id — skipping" >&2
    continue
  fi

  for modality in dna rna; do
    cohort_out="$OUTPUT_BASE/${modality}_data/$src_id"
    mkdir -p "$cohort_out"
    echo "=== $src_id [$modality] (candidates: $candidate_vcf) ==="

    for plate in ${SRC_PLATES[$src_id]}; do
      echo "  -- plate $plate --"
      if [[ "$modality" == "dna" ]]; then
        find_dna_cells "$plate"
      else
        find_rna_cells "$plate"
      fi | xargs -P "$PARALLEL_CELLS" -I{} bash -c \
        'run_one_cell "$@"' _ {} "$candidate_vcf" "$cohort_out" "$modality"
    done
  done
done


