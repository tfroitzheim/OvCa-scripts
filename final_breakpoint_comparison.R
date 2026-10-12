#!/usr/bin/env Rscript
# Canonical bulk-vs-scWGS breakpoint comparison, consolidating every filter

#   - both flanking segments >= min_bulk_segment_size
#   - real total-CN change >= min_cn_diff
#   - not centromere/telomere-supported
#   - near-zero-evidence micro-segments are MERGED into their larger
#     neighbor before breakpoints are defined (rather than just excluding
#     the boundary), recovering real matches that were previously dropped
#     purely because of an SV-marker artifact segment with 0 evidence
# Falls back to the best available ASCENT slot per patient and flags any
# patient not on the standard `masked` slot (currently: src35).
# still somewhat in progress so if u have better ideas let me know! 

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(purrr); library(tibble)
  library(GenomicRanges); library(rtracklayer); library(optparse)
})

sample_key <- tribble(
  ~src,    ~knum,
  "src16", "K19085-21", "src20", "K1190-22", "src21", "K1554-22",
  "src23", "K1862-22",  "src24", "K1883-22", "src26", "K2208-22",
  "src27", "K2765-22",  "src29", "K3784-22", "src32", "K11788-22",
  "src35", "K742-23",   "src39", "K6634-23"
)
bf_only <- c("src16", "src32") #src16 and 32 only had the "bf" (before filtering) RDS files available, so I used these

option_list <- list(
  make_option("--purple_root", type="character"),
  make_option("--ascent_root", type="character"),
  make_option("--chain", type="character"),
  make_option("--src", type="character"),
  make_option("--min_bulk_segment_size", type="integer", default=100000),
  make_option("--min_cn_diff", type="double", default=0.3),
  make_option("--min_evidence_count", type="integer", default=50),
  make_option("--tolerance", type="integer", default=100000),
  make_option("--n_perm", type="integer", default=1000),
  make_option("--outdir", type="character", default="/wrk/natgen_reviews/sarcoma_vcfs/bp_concordance")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
chain <- import.chain(opt$chain)

liftover_points <- function(chr, pos, chain) {
  if (length(pos) == 0) return(tibble(orig_idx = integer(0), chr = character(0), pos = numeric(0)))
  gr <- GRanges(chr, IRanges(pos, pos))
  lifted <- liftOver(gr, chain)
  n_per <- lengths(lifted)
  keep_idx <- which(n_per >= 1)
  if (length(keep_idx) == 0) return(tibble(orig_idx = integer(0), chr = character(0), pos = numeric(0)))
  first_each <- do.call(c, lapply(lifted[keep_idx], function(x) x[1]))
  tibble(orig_idx = keep_idx, chr = as.character(seqnames(first_each)), pos = start(first_each))
}

# Merge low-evidence micro-segments into their larger neighbor, per chromosome,
# before? computing any breakpoints. Repeats until no weak segments remain. Maybe re-evaluate?
merge_weak_segments <- function(purple, min_evidence_count, max_iter = 5000) {
  purple <- purple %>% arrange(chr, start) %>% as.data.frame()
  iter <- 0
  repeat {
    iter <- iter + 1
    if (iter > max_iter) { cat("  [merge] hit max_iter, stopping early\n"); break }
    weak <- which(purple$bafCount < min_evidence_count | purple$depthWindowCount < min_evidence_count)
    if (length(weak) == 0) break
    i <- weak[1]
    chr_i <- purple$chr[i]
    prev_i <- if (i > 1 && purple$chr[i-1] == chr_i) i - 1 else NA
    next_i <- if (i < nrow(purple) && purple$chr[i+1] == chr_i) i + 1 else NA
    if (is.na(prev_i) && is.na(next_i)) { purple <- purple[-i, ]; next }
    prev_size <- if (!is.na(prev_i)) purple$end[prev_i] - purple$start[prev_i] + 1 else -1
    next_size <- if (!is.na(next_i)) purple$end[next_i] - purple$start[next_i] + 1 else -1
    if (next_size >= prev_size) {
      purple$start[next_i] <- purple$start[i]
      purple$bafCount[next_i] <- purple$bafCount[next_i] + purple$bafCount[i]
      purple$depthWindowCount[next_i] <- purple$depthWindowCount[next_i] + purple$depthWindowCount[i]
      purple$segmentStartSupport[next_i] <- purple$segmentStartSupport[i]
    } else {
      purple$end[prev_i] <- purple$end[i]
      purple$bafCount[prev_i] <- purple$bafCount[prev_i] + purple$bafCount[i]
      purple$depthWindowCount[prev_i] <- purple$depthWindowCount[prev_i] + purple$depthWindowCount[i]
      purple$segmentEndSupport[prev_i] <- purple$segmentEndSupport[i]
    }
    purple <- purple[-i, ]
  }
  cat("  [merge] converged after", iter, "iterations,", nrow(purple), "segments remain\n")
  as_tibble(purple)
}

get_qualifying_bulk_breakpoints <- function(knum, chain, opt) {
  purple_file <- file.path(opt$purple_root, knum, "purple", paste0("DNA-T-", knum, ".purple.cnv.somatic.tsv"))
  if (!file.exists(purple_file)) return(tibble(chr = character(0), pos = numeric(0)))
  purple <- read_tsv(purple_file, show_col_types = FALSE) %>% mutate(chr = paste0("chr", chromosome))
  cat("  Raw segments:", nrow(purple), "\n")
  purple_merged <- merge_weak_segments(purple, opt$min_evidence_count)

  purple_bp <- purple_merged %>%
    group_by(chr) %>% arrange(start, .by_group = TRUE) %>%
    mutate(seg_size = end - start + 1, prev_seg_size = lag(seg_size),
           prev_cn = lag(copyNumber), cn_diff = abs(copyNumber - prev_cn),
           is_first = row_number() == 1) %>%
    ungroup() %>%
    filter(!is_first, seg_size >= opt$min_bulk_segment_size,
           prev_seg_size >= opt$min_bulk_segment_size, cn_diff >= opt$min_cn_diff,
           !segmentStartSupport %in% c("CENTROMERE", "TELOMERE")) %>%
    transmute(chr, pos = start)
  cat("  Qualifying breakpoints (post-merge):", nrow(purple_bp), "\n")
  if (nrow(purple_bp) == 0) return(tibble(chr = character(0), pos = numeric(0)))
  lifted <- liftover_points(purple_bp$chr, purple_bp$pos, chain)
  cat("  Post-liftover:", nrow(lifted), "\n")
  lifted %>% transmute(chr, pos)
}

match_bp <- function(query, subject, tol) {
  query$matched <- sapply(query$pos, function(p) any(abs(subject$pos[subject$chr == query$chr[which(query$pos==p)[1]]] - p) <= tol))
  query
}
# safer chromosome-aware matcher
match_bp_by_chr <- function(query, subject, tol) {
  query$matched <- FALSE
  for (cn in unique(query$chr)) {
    qi <- which(query$chr == cn); si <- subject$pos[subject$chr == cn]
    if (length(si) == 0) next
    query$matched[qi] <- sapply(query$pos[qi], function(p) any(abs(si - p) <= tol))
  }
  query
}

permutation_test <- function(bulk_bp, sc_bp, tol, chr_lengths, n_perm) {
  set.seed(42)
  observed <- mean(match_bp_by_chr(bulk_bp, sc_bp, tol)$matched)
  null_rates <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    sc_shuf <- sc_bp
    for (cn in unique(sc_bp$chr)) {
      idx <- which(sc_bp$chr == cn)
      len <- chr_lengths$length[chr_lengths$chr == cn]
      if (length(len) == 0 || length(idx) == 0) next
      sc_shuf$pos[idx] <- sample.int(len, length(idx), replace = TRUE)
    }
    null_rates[p] <- mean(match_bp_by_chr(bulk_bp, sc_shuf, tol)$matched)
  }
  list(observed = observed, null_mean = mean(null_rates), null_sd = sd(null_rates),
       p_value = (sum(null_rates >= observed) + 1) / (n_perm + 1))
}

chr_lengths <- tibble(
  chr = paste0("chr", c(1:22, "X", "Y")),
  length = c(248956422,242193529,198295559,190214555,181538259,170805979,
             159345973,145138636,138394717,133797422,135086622,133275309,
             114364328,107043718,101991189,90338345,83257441,80373285,
             58617616,64444167,46709983,50818468,156040895,57227415)
)

run_one <- function(src_id) {
  knum <- sample_key$knum[sample_key$src == src_id]
  cat("\nn", src_id, " <-> ", knum,
      "\nn", sep="")

  cat("  Computing bulk breakpoints (merged + filtered)...\n")
  bulk_bp <- get_qualifying_bulk_breakpoints(knum, chain, opt)

  ascent_file <- file.path(opt$ascent_root, paste0(src_id, ".phased_on_manual",
                            if (src_id %in% bf_only) "_bf" else "", ".Rds"))
  if (!file.exists(ascent_file)) { cat("  MISSING scWGS Rds -- skipping\n"); return(invisible(NULL)) }
  d <- read_rds(ascent_file)
  slot <- tail(names(d$segments), 1)
  if (slot != "masked") cat("  ** NOTE: this patient has no `masked` slot -- using `", slot,
                             "` instead. Not directly comparable to patients on `masked`. **\n", sep="")
  sc_segs <- d$segments[[slot]] %>% mutate(chr = as.character(chr))
  sc_bp <- sc_segs %>% group_by(chr) %>% arrange(start.pos, .by_group=TRUE) %>%
    filter(row_number() > 1) %>% ungroup() %>% transmute(chr, pos = start.pos)
  cat("  scWGS breakpoints (slot:", slot, "):", nrow(sc_bp), "\n")

  if (nrow(bulk_bp) == 0 || nrow(sc_bp) == 0) {
    cat("  Insufficient data for comparison -- skipping stats\n") #not sure
    return(tibble(src = src_id, knum = knum, slot = slot, n_bulk_bp = nrow(bulk_bp),
                   n_sc_bp = nrow(sc_bp), n_matched = NA, recovery_pct = NA,
                   null_mean_pct = NA, p_value = NA))
  }

  bulk_matched <- match_bp_by_chr(bulk_bp, sc_bp, opt$tolerance)
  n_matched <- sum(bulk_matched$matched)
  recovery <- 100 * n_matched / nrow(bulk_bp)
  cat(sprintf("  Recovery: %d / %d (%.1f%%)\n", n_matched, nrow(bulk_bp), recovery))

  perm <- permutation_test(bulk_bp, sc_bp, opt$tolerance, chr_lengths, opt$n_perm)
  cat(sprintf("  Permutation: null=%.1f%%+/-%.1f%%  p=%.4f\n",
              100*perm$null_mean, 100*perm$null_sd, perm$p_value))

  write_tsv(bulk_matched, file.path(opt$outdir, paste0(src_id, "_", knum, "-final_bulk_matched.tsv")))

  tibble(src = src_id, knum = knum, slot = slot, n_bulk_bp = nrow(bulk_bp), n_sc_bp = nrow(sc_bp),
         n_matched = n_matched, recovery_pct = round(recovery, 1),
         null_mean_pct = round(100*perm$null_mean, 1), p_value = perm$p_value)
}

if (!is.null(opt$src)) {
  print(run_one(opt$src))
} else {
  results <- map_dfr(sample_key$src, run_one)
  print(results)
  write_tsv(results, file.path(opt$outdir, "final_cohort_summary.tsv"))
  cat("\nCohort summry written to", file.path(opt$outdir, "final_cohort_summary.tsv"), "\n")
}
