#!/usr/bin/env Rscript
# Venn diagrams (VennDiagram) from peak BEDs. Partition counts match peaks_ops UpSet
# (mutually exclusive regions via count_distinct logic from peaks_ops_upsetR.R).
#
# Usage: Rscript peak_vennDiagram.R [--names "A,B,..."] OUTPUT_BASE BED1 BED2 [BED3 BED4]
# Supports n = 2, 3, or 4 only.

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

if (length(args) < 3) {
  stop("Usage: Rscript peak_vennDiagram.R [--names \"N1,N2,...\"] OUTPUT_BASE BED1 BED2 [BED3 BED4]")
}

out_base <- args[1]
bed_files <- args[-1]
n <- length(bed_files)

if (n < 2 || n > 4) {
  message("VennDiagram supports 2–4 sets; got ", n, ". Skipping Venn plot.")
  quit(save = "no", status = 0)
}

if (!is.null(set_names_param)) {
  if (length(set_names_param) != n) {
    stop("--names must have ", n, " entries, got ", length(set_names_param))
  }
  set_names <- trimws(set_names_param)
} else {
  set_names <- sub("\\.bed$", "", basename(bed_files))
}

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  stop("VennDiagram not installed. Install: install.packages('VennDiagram')")
}

count_bed_lines <- function(path) {
  if (!file.exists(path)) return(0)
  d <- tryCatch(read.table(path, sep = "\t", header = FALSE, comment.char = ""),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(0)
  sum(nchar(trimws(as.character(d[, 1]))) > 0 & d[, 1] != "")
}

# Same as peaks_ops_upsetR.R count_distinct (partition sizes for UpSet bars)
count_distinct <- function(idx) {
  in_sets <- bed_files[idx]
  out_idx <- setdiff(seq_len(n), idx)
  tmp_intersect <- tempfile(fileext = ".bed")
  if (length(idx) == 1) {
    exit_code <- system(paste0("cp ", shQuote(in_sets[1]), " ", shQuote(tmp_intersect)))
  } else {
    cmd <- paste0("bedtools intersect -a ", shQuote(in_sets[1]), " -b ", shQuote(in_sets[2]))
    for (f in in_sets[-c(1, 2)]) {
      cmd <- paste0(cmd, " | bedtools intersect -a stdin -b ", shQuote(f))
    }
    exit_code <- system(paste0(cmd, " > ", shQuote(tmp_intersect)))
  }
  if (exit_code != 0) {
    unlink(tmp_intersect, force = TRUE)
    return(0)
  }
  if (length(out_idx) == 0) {
    cnt <- count_bed_lines(tmp_intersect)
    unlink(tmp_intersect, force = TRUE)
    return(cnt)
  }
  out_sets <- bed_files[out_idx]
  tmp_union <- tempfile(fileext = ".bed")
  tmp_result <- tempfile(fileext = ".bed")
  files_cat <- paste(shQuote(out_sets), collapse = " ")
  exit_code <- system(paste0("cat ", files_cat,
                               " | sort -k1,1 -k2,2n | bedtools merge -i stdin > ", shQuote(tmp_union)))
  if (exit_code != 0) {
    unlink(c(tmp_intersect, tmp_union, tmp_result), force = TRUE)
    return(0)
  }
  exit_code <- system(paste0("bedtools subtract -a ", shQuote(tmp_intersect),
                               " -b ", shQuote(tmp_union), " > ", shQuote(tmp_result)))
  cnt <- 0
  if (exit_code == 0) cnt <- count_bed_lines(tmp_result)
  unlink(c(tmp_intersect, tmp_union, tmp_result), force = TRUE)
  cnt
}

# All non-empty subsets as index vectors
all_subsets <- function(nn) {
  out <- list()
  for (k in seq_len(nn)) {
    out <- c(out, combn(nn, k, simplify = FALSE))
  }
  out
}

subs <- all_subsets(n)
part <- vapply(subs, count_distinct, numeric(1))
names(part) <- vapply(subs, function(idx) paste(idx, collapse = ","), character(1))

# Sum partition sizes where subset S contains all of required indices
sum_in_sets <- function(required) {
  tot <- 0
  for (j in seq_along(subs)) {
    s <- subs[[j]]
    if (all(required %in% s)) tot <- tot + part[j]
  }
  tot
}

area_i <- function(i) sum_in_sets(i)

# Non-negative clamp for VennDiagram
clamp_nonneg <- function(x) pmax(0, x)

if (n == 2) {
  ab <- clamp_nonneg(sum_in_sets(c(1, 2)))
  a1 <- clamp_nonneg(area_i(1))
  a2 <- clamp_nonneg(area_i(2))
  if (ab > min(a1, a2)) ab <- min(a1, a2)
  vg <- VennDiagram::draw.pairwise.venn(
    area1 = a1, area2 = a2, cross.area = ab,
    category = set_names,
    fill = c("#7FCDBB", "#2E86AB"),
    alpha = 0.6, cat.fontface = "bold", margin = 0.1
  )
  pdf(paste0(out_base, "_venn.pdf"), width = 7, height = 7)
  grid::grid.draw(vg)
  invisible(dev.off())
} else if (n == 3) {
  n123 <- clamp_nonneg(sum_in_sets(c(1, 2, 3)))
  n12 <- clamp_nonneg(sum_in_sets(c(1, 2)))
  n13 <- clamp_nonneg(sum_in_sets(c(1, 3)))
  n23 <- clamp_nonneg(sum_in_sets(c(2, 3)))
  a1 <- clamp_nonneg(area_i(1))
  a2 <- clamp_nonneg(area_i(2))
  a3 <- clamp_nonneg(area_i(3))
  vg <- VennDiagram::draw.triple.venn(
    area1 = a1, area2 = a2, area3 = a3,
    n12 = n12, n13 = n13, n23 = n23, n123 = n123,
    category = set_names,
    fill = c("#7FCDBB", "#2E86AB", "#E94F37"),
    alpha = 0.5, cat.fontface = "bold", margin = 0.08
  )
  pdf(paste0(out_base, "_venn.pdf"), width = 8, height = 8)
  grid::grid.draw(vg)
  invisible(dev.off())
} else {
  # n == 4; draw.quad.venn parameter order
  n1234 <- clamp_nonneg(sum_in_sets(c(1, 2, 3, 4)))
  n123 <- clamp_nonneg(sum_in_sets(c(1, 2, 3)))
  n124 <- clamp_nonneg(sum_in_sets(c(1, 2, 4)))
  n134 <- clamp_nonneg(sum_in_sets(c(1, 3, 4)))
  n234 <- clamp_nonneg(sum_in_sets(c(2, 3, 4)))
  n12 <- clamp_nonneg(sum_in_sets(c(1, 2)))
  n13 <- clamp_nonneg(sum_in_sets(c(1, 3)))
  n14 <- clamp_nonneg(sum_in_sets(c(1, 4)))
  n23 <- clamp_nonneg(sum_in_sets(c(2, 3)))
  n24 <- clamp_nonneg(sum_in_sets(c(2, 4)))
  n34 <- clamp_nonneg(sum_in_sets(c(3, 4)))
  a1 <- clamp_nonneg(area_i(1))
  a2 <- clamp_nonneg(area_i(2))
  a3 <- clamp_nonneg(area_i(3))
  a4 <- clamp_nonneg(area_i(4))
  vg <- VennDiagram::draw.quad.venn(
    area1 = a1, area2 = a2, area3 = a3, area4 = a4,
    n12 = n12, n13 = n13, n14 = n14, n23 = n23, n24 = n24, n34 = n34,
    n123 = n123, n124 = n124, n134 = n134, n234 = n234, n1234 = n1234,
    category = set_names,
    fill = c("#7FCDBB", "#2E86AB", "#E94F37", "#44AF69"),
    alpha = 0.45, cat.fontface = "bold", margin = 0.06
  )
  pdf(paste0(out_base, "_venn.pdf"), width = 9, height = 9)
  grid::grid.draw(vg)
  invisible(dev.off())
}
