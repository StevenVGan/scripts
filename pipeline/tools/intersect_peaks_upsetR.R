#!/usr/bin/env Rscript
# UpSet plot using UpSetR.
# Usage: Rscript intersect_peaks_upsetR.R [--names "N1,N2,..."] OUTPUT_BASE BED1 BED2 [BED3 ...]

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
  stop("Usage: Rscript intersect_peaks_upsetR.R [--names \"N1,N2,...\"] OUTPUT_BASE BED1 BED2 [BED3 ...]")
}

out_base <- args[1]
bed_files <- args[-1]
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

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  stop("UpSetR not installed. Install: install.packages('UpSetR', repos='https://cloud.r-project.org')")
}

# Count valid BED lines (avoids wc -l counting empty lines or headers)
count_bed_lines <- function(path) {
  if (!file.exists(path)) return(0)
  d <- tryCatch(read.table(path, sep = "\t", header = FALSE, comment.char = ""), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(0)
  sum(nchar(trimws(as.character(d[, 1]))) > 0 & d[, 1] != "")
}

# bedtools intersect -a A -b B -b C returns A ∩ (B ∪ C), not A ∩ B ∩ C.
# Chain intersects to require overlap with ALL files.
count_intersection <- function(idx) {
  if (length(idx) == 1) {
    count_bed_lines(bed_files[idx])
  } else {
    others <- bed_files[idx[-1]]
    cmd <- paste0("bedtools intersect -a ", shQuote(bed_files[idx[1]]), " -b ", shQuote(others[1]))
    for (f in others[-1]) {
      cmd <- paste0(cmd, " | bedtools intersect -a stdin -b ", shQuote(f))
    }
    tmp <- tempfile(fileext = ".bed")
    cnt <- 0
    exit_code <- system(paste0(cmd, " > ", shQuote(tmp)))
    if (exit_code == 0) {
      cnt <- count_bed_lines(tmp)
    } else {
      warning("bedtools intersect failed (exit ", exit_code, ") for combo: ",
              paste(set_names[idx], collapse = "&"))
    }
    unlink(tmp, force = TRUE)
    cnt
  }
}

# Generate all non-empty subsets (2^n - 1)
combos <- list()
for (k in seq_len(n)) {
  cmb <- combn(n, k, simplify = FALSE)
  combos <- c(combos, cmb)
}

expr <- numeric(length(combos))
names(expr) <- vapply(combos, function(idx) {
  paste(set_names[idx], collapse = "&")
}, character(1))

for (i in seq_along(combos)) {
  expr[i] <- count_intersection(combos[[i]])
}

# Skip if no multi-set overlap (all 2+ way intersections are 0)
multi_way <- vapply(combos, function(idx) length(idx) >= 2, logical(1))
if (!any(expr[multi_way] > 0)) {
  message("No overlapping regions across sets; skipping UpSet plot")
  quit(save = "no", status = 0)
}

# Remove zero counts
expr <- expr[expr > 0]
if (length(expr) == 0) {
  message("No overlapping regions; skipping UpSet plot")
  quit(save = "no", status = 0)
}

# Colors: one per intersection degree (2-way, 3-way, 4-way, ...)
degree_colors <- c(
  "1" = "#7FCDBB",   # 1-way (teal)
  "2" = "#2E86AB",   # 2-way (blue)
  "3" = "#E94F37",   # 3-way (red)
  "4" = "#44AF69",   # 4-way (green)
  "5" = "#F18F01",   # 5-way (orange)
  "6" = "#5C4D7D"    # 6+ way (purple)
)
sets_bar_col <- c("turquoise4")
matrix_col <- c("slateblue4")
shade_col <- c("wheat4")

# Build queries: one per bar, colored by intersection degree
queries <- list()
for (i in seq_along(expr)) {
  sets_in_combo <- strsplit(names(expr)[i], "&", fixed = TRUE)[[1]]
  deg <- length(sets_in_combo)
  col <- degree_colors[as.character(deg)]
  if (is.na(col)) col <- "gray50"
  queries[[i]] <- list(
    query = UpSetR::intersects,
    params = as.list(sets_in_combo),
    color = col,
    active = TRUE
  )
}

# Text scale: c(main title, main ticks, set title, set ticks, set names, numbers)
text_scale <- c(1.5, 1.25, 1.25, 1, 1.7, 1.5)
upset_data <- UpSetR::fromExpression(expr)

pdf(paste0(out_base, ".pdf"), width = 8, height = 6)
UpSetR::upset(upset_data,
              mainbar.y.label = "Peak intersections", sets.x.label = "Peak count",
              text.scale = text_scale,
              sets.bar.color = sets_bar_col,
              matrix.color = matrix_col, shade.color = shade_col,
              queries = queries, query.legend = "none")
invisible(dev.off())
