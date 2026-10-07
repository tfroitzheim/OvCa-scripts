#!/usr/bin/env Rscript
# DPClust subclonal reconstruction, hg38 pipeline.
# Reuses workspace_hg38/pyclone_input/<src>.tsv (same filtered mutation set as PyClone-VI,
# hg38 coordinates) so DPClust and PyClone-VI results are directly comparable.

library(DPClust)

PYCLONE_INPUT_DIR <- "/wrk/natgen_reviews/sarcoma_vcfs/workspace_hg38/pyclone_input"
OUT_BASE <- "/wrk/natgen_reviews/sarcoma_vcfs/workspace_hg38/dpclust_output"
CLINICAL_DIR <- "/wrk/natgen_reviews/sarcoma_vcfs/data/clinical_bulk_data"

SRC_TO_K <- c(src16="K19085-21", src20="K1190-22", src21="K1554-22", src23="K1862-22",
              src26="K2208-22", src27="K2765-22", src32="K11788-22", src35="K742-23",
              src39="K6634-23")

# --- purity lookup, pulled directly from each patient's purple.purity.tsv ---
get_purity <- function(k_id) {
  path <- file.path(CLINICAL_DIR, k_id, "purple", paste0("DNA-T-", k_id, ".purple.purity.tsv"))
  if (!file.exists(path)) return(NA)
  row <- read.delim(path, sep = "\t", nrows = 1)
  row$purity
}
purity_lookup <- setNames(sapply(SRC_TO_K, get_purity), names(SRC_TO_K))

# --- build DPClust's required input columns from our filtered PyClone-VI input tsv ---
# CRITICAL FIX: no.chrs.bearing.mut (multiplicity) must be estimated per-mutation,

build_dpclust_input <- function(pyclone_input_tsv) {
  df <- read.delim(pyclone_input_tsv, sep = "\t", stringsAsFactors = FALSE)

  parts <- strsplit(df$mutation_id, "_")
  df$chr <- sapply(parts, `[`, 1)
  df$end <- as.integer(sapply(parts, `[`, 2))

  df$WT.count  <- df$ref_counts
  df$mut.count <- df$alt_counts
  df$subclonal.CN <- df$major_cn + df$minor_cn

  vaf <- df$mut.count / (df$mut.count + df$WT.count)
  purity <- df$tumour_content

  raw_mcn <- (vaf / purity) * df$subclonal.CN

  # multiplicity: round to nearest integer, bounded [1, major_cn] --> not sure if this is absolutely the best way to go about this. From what I saw both methods need integers to run, but maybe you have better ideas how to approach this??
  
  df$no.chrs.bearing.mut <- pmin(pmax(round(raw_mcn), 1), df$major_cn)

  df$mutation.copy.number <- raw_mcn
  df$subclonal.fraction <- pmin(df$mutation.copy.number / df$no.chrs.bearing.mut, 1)

  df[, c("chr", "end", "WT.count", "mut.count", "subclonal.CN",
         "mutation.copy.number", "subclonal.fraction", "no.chrs.bearing.mut")]
}

# run DPClust for one patient 
run_dpclust_fixed <- function(src_id, purity, no.iters = 1250, no.iters.burn.in = 250) {
  input_tsv <- file.path(PYCLONE_INPUT_DIR, paste0(src_id, ".tsv"))
  if (!file.exists(input_tsv)) {
    cat(sprintf("%s: no filtered input found, skipping\n", src_id))
    return(NULL)
  }

  dp_input <- build_dpclust_input(input_tsv)
  out_dir <- file.path(OUT_BASE, src_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  as_col <- function(x) matrix(x, ncol = 1)  # single-sample: DPClust needs 1-column matrices

  result <- tryCatch({
    DPClust:::DirichletProcessClustering(
      mutCount = as_col(dp_input$mut.count),
      WTCount = as_col(dp_input$WT.count),
      totalCopyNumber = as_col(dp_input$subclonal.CN),
      copyNumberAdjustment = as_col(dp_input$no.chrs.bearing.mut),
      mutation.copy.number = as_col(dp_input$mutation.copy.number),
      cellularity = purity,
      output_folder = out_dir,
      no.iters = no.iters,
      no.iters.burn.in = no.iters.burn.in,
      subsamplesrun = src_id,
      samplename = src_id,
      conc_param = 1,
      cluster_conc = 5,
      mut.assignment.type = 1,
      most.similar.mut = NA,
      mutationTypes = "SNV",
      max.considered.clusters = 10
    )
  }, error = function(e) {
    cat(sprintf("%s FAILED: %s\n", src_id, e$message))
    NULL
  })

  if (!is.null(result)) {
    mutation_ids <- read.delim(input_tsv, sep = "\t")$mutation_id
    enriched <- data.frame(
      mutation_id = mutation_ids,
      dpclust_cluster = result$best.node.assignments,
      dpclust_confidence = result$best.assignment.likelihoods,
      subclonal_fraction = dp_input$subclonal.fraction
    )
    write.csv(enriched, file.path(out_dir, paste0(src_id, "_mutation_cluster_map.csv")), row.names = FALSE)
    saveRDS(result, file.path(out_dir, paste0(src_id, "_dpclust_result.Rds")))
    cat(sprintf("%s: done, %d clusters used\n", src_id, length(unique(result$best.node.assignments))))
  }
  result
}

# --- run across all 9 patients ---
srcs <- names(SRC_TO_K)
dp_results <- list()
for (s in srcs) {
  dp_results[[s]] <- run_dpclust_fixed(s, purity_lookup[[s]])
}

cat("\n=== Cluster locations per patient ===\n")
for (s in srcs) {
  if (is.null(dp_results[[s]])) next
  cat(sprintf("\n%s:\n", s))
  print(dp_results[[s]]$cluster.locations)
}
