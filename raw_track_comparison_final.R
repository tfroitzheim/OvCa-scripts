#!/usr/bin/env Rscript

# Same raw-track plots as raw_track_comparison.R, but bulk breakpoints are
# now defined on the MERGED segment table (weak-evidence micro-segments
# folded into their larger neighbor first) -- matching the canonical
# final_breakpoint_comparison.R logic.


# this pulls the WGS data --> e.g. cobalt ratio, mirrored BAF, I guess it very much depend sin what format your data comes?


suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(patchwork); library(GenomicRanges); library(rtracklayer)
  library(optparse)
})

sample_key <- tribble(
  ~src,    ~knum,
  "src16", "K19085-21", "src20", "K1190-22", "src21", "K1554-22",
  "src23", "K1862-22",  "src24", "K1883-22", "src26", "K2208-22",
  "src27", "K2765-22",  "src29", "K3784-22", "src32", "K11788-22",
  "src35", "K742-23",   "src39", "K6634-23"
)
bf_only <- c("src16", "src32")
chr_list <- c(as.character(1:22), "X", "Y")

option_list <- list(
  make_option("--purple_root", type="character"),
  make_option("--ascent_root", type="character"),
  make_option("--chain", type="character"),
  make_option("--src", type="character"),
  make_option("--downsample_n", type="integer", default=4000),
  make_option("--min_bulk_segment_size", type="integer", default=100000),
  make_option("--min_cn_diff", type="double", default=0.3),
  make_option("--min_evidence_count", type="integer", default=50),
  make_option("--tolerance", type="integer", default=100000),
  make_option("--max_labeled_excluded", type="integer", default=15),
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

thin <- function(df, n) {
  if (nrow(df) <= n) return(df)
  df[round(seq(1, nrow(df), length.out = n)), ]
}

merge_weak_segments <- function(purple, min_evidence_count, max_iter = 5000) {
  purple <- purple %>% arrange(chr, start) %>% as.data.frame()
  iter <- 0
  repeat {
    iter <- iter + 1
    if (iter > max_iter) break
    weak <- which(purple$bafCount < min_evidence_count | purple$depthWindowCount < min_evidence_count)
    if (length(weak) == 0) break
    i <- weak[1]; chr_i <- purple$chr[i]
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
  as_tibble(purple)
}

load_all_bulk_breakpoints <- function(knum, chain, opt) {
  purple_file <- file.path(opt$purple_root, knum, "purple", paste0("DNA-T-", knum, ".purple.cnv.somatic.tsv"))
  if (!file.exists(purple_file)) return(tibble(chr = character(0), pos = numeric(0),
                                                 qualifies = logical(0), exclude_reason = character(0),
                                                 method = character(0)))
  purple_raw <- read_tsv(purple_file, show_col_types = FALSE) %>% mutate(chr = paste0("chr", chromosome))
  purple <- merge_weak_segments(purple_raw, opt$min_evidence_count) %>%
    group_by(chr) %>% arrange(start, .by_group = TRUE) %>%
    mutate(seg_size = end - start + 1, prev_seg_size = lag(seg_size),
           prev_cn = lag(copyNumber), cn_diff = abs(copyNumber - prev_cn),
           is_first = row_number() == 1) %>%
    ungroup() %>%
    filter(!is_first) %>%
    mutate(
      qualifies = seg_size >= opt$min_bulk_segment_size & prev_seg_size >= opt$min_bulk_segment_size &
                  cn_diff >= opt$min_cn_diff & !segmentStartSupport %in% c("CENTROMERE", "TELOMERE"),
      exclude_reason = case_when(
        segmentStartSupport %in% c("CENTROMERE", "TELOMERE") ~ "centromere/telomere",
        seg_size < opt$min_bulk_segment_size | prev_seg_size < opt$min_bulk_segment_size ~ "segment too small",
        cn_diff < opt$min_cn_diff ~ "no real CN change",
        TRUE ~ NA_character_
      )
    ) %>%
    transmute(chr, pos = start, qualifies, exclude_reason, method)

  if (nrow(purple) == 0) return(tibble(chr = character(0), pos = numeric(0),
                                        qualifies = logical(0), exclude_reason = character(0), method = character(0)))
  lifted <- liftover_points(purple$chr, purple$pos, chain)
  lifted %>% inner_join(purple %>% mutate(orig_idx = row_number()) %>%
                           select(orig_idx, qualifies, exclude_reason, method), by = "orig_idx") %>%
    select(-orig_idx)
}

build_chr_plot <- function(knum, src, chr, d, chain, downsample_n, bulk_bp_all, tol, max_labeled) {
  chr_hg38 <- paste0("chr", chr)

  cobalt_file <- file.path(opt$purple_root, knum, "cobalt", paste0("DNA-T-", knum, ".cobalt.ratio.tsv.gz"))
  cobalt <- read_tsv(cobalt_file, show_col_types = FALSE) %>%
    filter(chromosome == chr, tumorGCRatio > -1) %>% thin(downsample_n) %>% mutate(orig_idx = row_number())
  if (nrow(cobalt) == 0) return(NULL)
  lift_c <- liftover_points(paste0("chr", cobalt$chromosome), cobalt$position, chain)
  cobalt_hg38 <- lift_c %>% inner_join(cobalt %>% select(orig_idx, tumorGCRatio), by = "orig_idx") %>%
    transmute(pos_hg38 = pos, tumorGCRatio)

  amber_file <- file.path(opt$purple_root, knum, "amber", paste0("DNA-T-", knum, ".amber.baf.tsv.gz"))
  amber <- read_tsv(amber_file, show_col_types = FALSE) %>%
    filter(chromosome == chr) %>% thin(downsample_n) %>% mutate(orig_idx = row_number())
  amber_hg38 <- if (nrow(amber) > 0) {
    lift_a <- liftover_points(paste0("chr", amber$chromosome), amber$position, chain)
    lift_a %>% inner_join(amber %>% select(orig_idx, tumorModifiedBAF), by = "orig_idx") %>%
      transmute(pos_hg38 = pos, tumorModifiedBAF)
  } else tibble(pos_hg38 = numeric(0), tumorModifiedBAF = numeric(0))

  slot <- tail(names(d$segments), 1)
  sc_segs <- d$segments[[slot]]
  sc_segs_gr <- sc_segs %>% mutate(chr = as.character(chr), row_id = row_number()) %>%
    transmute(chr, start_bp = start.pos, end_bp = end.pos, row_id) %>% filter(chr == chr_hg38)
  clone_ids <- d$clones$clone_id
  seg_means <- d$clones$segment_means
  sc_df <- map_dfr(seq_along(clone_ids), function(i) {
    cn_vec <- seg_means[[i]]
    if (is.null(cn_vec) || length(cn_vec) != nrow(sc_segs)) return(NULL)
    sc_segs_gr %>% mutate(cn = cn_vec[row_id], clone = clone_ids[i])
  })

  sc_bp_pos <- if (nrow(sc_segs_gr) > 1) {
    sc_segs_gr %>% arrange(start_bp) %>% filter(row_number() > 1) %>% pull(start_bp)
  } else numeric(0)

  bp_chr <- bulk_bp_all %>% filter(chr == chr_hg38)
  bp_qual <- bp_chr %>% filter(qualifies)
  bp_excl <- bp_chr %>% filter(!qualifies)

  n_matched <- if (nrow(bp_qual) > 0 && length(sc_bp_pos) > 0) {
    sum(sapply(bp_qual$pos, function(p) any(abs(sc_bp_pos - p) <= tol)))
  } else 0
  stat_line <- sprintf("Bulk breakpoints matched: %d / %d qualifying  |  excluded: %d  |  scWGS slot: %s",
                        n_matched, nrow(bp_qual), nrow(bp_excl), slot)

  x_range <- range(c(cobalt_hg38$pos_hg38, amber_hg38$pos_hg38, sc_df$start_bp, sc_df$end_bp), na.rm = TRUE)
  if (any(!is.finite(x_range))) return(NULL)

  bp_layer_sc   <- if (length(sc_bp_pos) > 0)
    geom_vline(xintercept = sc_bp_pos, color = "#2a2a2a", linewidth = 0.3, linetype = "dotted") else NULL
  bp_layer_excl <- if (nrow(bp_excl) > 0)
    geom_vline(data = bp_excl, aes(xintercept = pos, color = method), linewidth = 0.3, linetype = "22", alpha = 0.6) else NULL
  bp_layer_qual <- if (nrow(bp_qual) > 0)
    geom_vline(xintercept = bp_qual$pos, color = "#791313", linewidth = 0.5) else NULL
  bp_labels <- if (nrow(bp_excl) > 0 && nrow(bp_excl) <= max_labeled)
    geom_text(data = bp_excl, aes(x = pos, y = Inf, label = method), angle = 90,
              hjust = 1.05, vjust = -0.3, size = 2.2, color = "grey40") else NULL

  method_colors <- c("BAF_WEIGHTED" = "#f0997b", "STRUCTURAL_VARIANT" = "#93b8ed",
                      "LONG_ARM" = "#c0dd97", "NON_DIPLOID" = "#ecc9f0")

  p_cobalt <- ggplot(cobalt_hg38, aes(pos_hg38, tumorGCRatio)) +
    bp_layer_sc + bp_layer_excl + bp_layer_qual + bp_labels +
    geom_point(size = 0.4, alpha = 0.3, color = "#2a78d6") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
    scale_color_manual(values = method_colors, name = "Excluded (method)", na.translate = FALSE) +
    coord_cartesian(xlim = x_range, ylim = c(0, 4)) +
    labs(title = paste0(src, " / ", knum, " -- chr", chr,
                         "  (dark red = qualifying bulk bp [merged], light = excluded, dotted = scWGS bp)"),
         subtitle = stat_line, y = "COBALT ratio", x = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "top", legend.key.size = unit(10, "pt"))

  p_amber <- ggplot(amber_hg38, aes(pos_hg38, tumorModifiedBAF)) +
    bp_layer_sc + bp_layer_excl + bp_layer_qual +
    geom_point(size = 0.4, alpha = 0.3, color = "#eb6834") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60") +
    scale_color_manual(values = method_colors, guide = "none") +
    coord_cartesian(xlim = x_range, ylim = c(0.4, 1)) +
    labs(y = "AMBER mirrored BAF", x = NULL) +
    theme_minimal(base_size = 11)

  p_ascent <- if (nrow(sc_df) > 0) {
    ggplot(sc_df, aes(x = start_bp, xend = end_bp, y = cn, yend = cn)) +
      bp_layer_sc + bp_layer_excl + bp_layer_qual +
      geom_segment(aes(color = clone), linewidth = 1) + coord_cartesian(xlim = x_range) +
      labs(y = "ASCENT clone CN ratio", x = paste0("chr", chr, " position (hg38)")) +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  } else {
    ggplot() + labs(y = "ASCENT clone CN ratio (no data)", x = paste0("chr", chr, " position (hg38)")) +
      theme_minimal(base_size = 11)
  }

  p_cobalt / p_amber / p_ascent + plot_layout(heights = c(1.2, 1, 1.3))
}

run_patient <- function(src) {
  knum <- sample_key$knum[sample_key$src == src]
  cat("\n=== ", src, " <-> ", knum, " ===\n", sep="")
  ascent_file <- file.path(opt$ascent_root, paste0(src, ".phased_on_manual",
                            if (src %in% bf_only) "_bf" else "", ".Rds"))
  if (!file.exists(ascent_file)) { cat("  MISSING scWGS Rds -- skipping patient\n"); return(invisible(NULL)) }
  d <- read_rds(ascent_file)

  cat("  Loading + merging + filtering bulk breakpoints...\n")
  bulk_bp_all <- load_all_bulk_breakpoints(knum, chain, opt)
  if (nrow(bulk_bp_all) == 0) {
    cat("  No bulk breakpoints available -- skipping patient\n")
    return(invisible(NULL))
  }
  cat("  Total lifted:", nrow(bulk_bp_all), " | qualifying:", sum(bulk_bp_all$qualifies),
      " | excluded:", sum(!bulk_bp_all$qualifies), "\n")

  outfile <- file.path(opt$outdir, paste0(src, "_", knum, "-raw_track_final.pdf"))
  pdf(outfile, width = 12, height = 8)
  n_pages <- 0
  for (chr in chr_list) {
    p <- tryCatch(build_chr_plot(knum, src, chr, d, chain, opt$downsample_n, bulk_bp_all,
                                  opt$tolerance, opt$max_labeled_excluded),
                  error = function(e) { cat("  chr", chr, ": ERROR --", conditionMessage(e), "\n"); NULL })
    if (!is.null(p)) { print(p); n_pages <- n_pages + 1 } else { cat("  chr", chr, ": skipped (no data)\n") }
  }
  dev.off()
  cat("  Wrote", n_pages, "pages to", outfile, "\n")
}

if (!is.null(opt$src)) {
  run_patient(opt$src)
} else {
  walk(sample_key$src, run_patient)
}
