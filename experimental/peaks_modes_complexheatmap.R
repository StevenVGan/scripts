#!/usr/bin/env Rscript
# UpSet plot and BED output using ComplexHeatmap (distinct / intersect / union modes).
#
# Usage: Rscript peaks_modes_complexheatmap.R [--names "N1,N2,..."] MODE OUTPUT_BASE BED1 BED2 [BED3 ...]
#
# MODE: distinct | intersect | union
#   distinct:  Mutually exclusive partitions (1=in, 0=not in). "all" = intersect of all.
#   intersect: 1=in, 0=ignored. "all" = intersect of all.
#   union:    1=in, 0=ignored, multiple 1s = OR. "all" = union of all.
#
# Output: OUTPUT_BASE.pdf (UpSet plot) and OUTPUT_BASE.bed (regions for "all" combination)

raw_args <- commandArgs(trailingOnly = TRUE)
set_names_param <- NULL
args <- character(0)
i <- 1
while (i <= length(raw_args)) {
  if (raw_args[i] == "--names" && i < length(raw_args)) {
    set_names_param <- strsplit(raw_args[i + 1], ",", fixed = TRUE)[[1]]
    i <- i + 2
    next
  }
  args <- c(args, raw_args[i])
  i <- i + 1
}

if (length(args) < 4) {
  stop("Usage: Rscript peaks_modes_complexheatmap.R [--names \"N1,N2,...\"] MODE OUTPUT_BASE BED1 BED2 [BED3 ...]")
}

mode <- match.arg(args[1], c("distinct", "intersect", "union"))
out_base <- args[2]
bed_files <- args[-(1:2)]
n <- length(bed_files)

if (!is.null(set_names_param)) {
  if (length(set_names_param) != n) {
    stop("--names must have ", n, " entries (one per BED file), got ", length(set_names_param))
  }
  set_names <- substr(trimws(set_names_param), 1, 25)
} else {
  set_names <- sub("\\.(bed|narrowPeak|broadPeak)$", "", basename(bed_files))
  set_names <- substr(set_names, 1, 25)
}

# Dependencies
for (pkg in c("ComplexHeatmap", "GenomicRanges", "IRanges")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " not installed. Install: BiocManager::install(c('ComplexHeatmap','GenomicRanges'))")
  }
}
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(GenomicRanges)
})

# Read BED files as GRanges
read_bed_gr <- function(path) {
  d <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  GRanges(seqnames = d[, 1], ranges = IRanges::IRanges(start = d[, 2], end = d[, 3]))
}

peak_gr_list <- lapply(bed_files, read_bed_gr)
names(peak_gr_list) <- set_names

# Use merged union regions as elements so size = peak count (not base pairs).
# GRanges input uses sum(width) = bp; list of vectors uses length = count.
all_gr <- suppressWarnings(Reduce(GenomicRanges::union, peak_gr_list))
merged <- GenomicRanges::reduce(all_gr)
merged_ids <- paste0(seqnames(merged), ":", start(merged), "-", end(merged))
names(merged) <- merged_ids

# For each set: which merged regions overlap it?
peak_list <- lapply(peak_gr_list, function(gr) {
  ov <- GenomicRanges::findOverlaps(merged, gr, type = "any")
  unique(merged_ids[queryHits(ov)])
})
names(peak_list) <- set_names

# Make combination matrix (size = count of regions)
m <- make_comb_mat(peak_list, mode = mode)

# Code for "all" combination (1 in every set)
all_code <- paste(rep("1", n), collapse = "")
if (!all_code %in% comb_name(m)) {
  stop("Combination ", all_code, " not found. Check input sets.")
}

# Extract region IDs for "all" combination, convert to BED
ids_all <- extract_comb(m, all_code)
if (length(ids_all) > 0) {
  # Parse "chr:start-end" back to BED
  parsed <- strsplit(as.character(ids_all), "[:-]")
  df <- do.call(rbind, lapply(parsed, function(x) data.frame(chr = x[1], start = as.integer(x[2]), end = as.integer(x[3]))))
  write.table(df, paste0(out_base, ".bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  message("BED (", mode, "): ", out_base, ".bed (", length(ids_all), " regions)")
} else {
  write.table(data.frame(character(), integer(), integer()), paste0(out_base, ".bed"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  message("BED (", mode, "): ", out_base, ".bed (0 regions)")
}

# Setting colors (match intersect_peaks_upsetR.R)
# This can also be done with hexadecimal
main_bar_col <- "violetred4"
sets_bar_col <- "turquoise4"
matrix_col <- "slateblue4"
shade_col <- adjustcolor("wheat4", alpha.f = 0.35)  # transparent shade

# UpSet plot: transparent shade, no bar borders, adjusted ratio
# shade_col with alpha for transparency; gpar(col=NA) removes bar borders
pdf(paste0(out_base, ".pdf"), width = 10, height = 10)
UpSet(m, column_title = paste0("UpSet plot (", mode, " mode)"),
      comb_col = matrix_col, bg_col = shade_col, bg_pt_col = shade_col,
      top_annotation = upset_top_annotation(m, gp = gpar(fill = main_bar_col, col = NA)),
      right_annotation = upset_right_annotation(m, gp = gpar(fill = sets_bar_col, col = NA)))
dev.off()
message("UpSet plot -> ", out_base, ".pdf")
