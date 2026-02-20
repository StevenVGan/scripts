#!/usr/bin/env Rscript
# UpSet plot with mode-specific counting (distinct, intersect, union).
# Usage: Rscript peaks_ops_upsetR.R --mode MODE [--names "N1,N2,..."] OUTPUT_BASE BED1 BED2 [BED3 ...]
#
# Modes:
#   distinct:  1=in, 0=NOT in. 1 1 0 = setdiff(intersect(A,B), C). Mutually exclusive partitions.
#   intersect: 1=in, 0=ignored. 1 1 0 = intersect(A,B). Bars can overlap.
#   union:     1=in, 0=ignored. 1 1 0 = union(A,B). Bars can overlap.

raw_args <- commandArgs(trailingOnly = TRUE)
set_names_param <- NULL
mode_param <- NULL
args <- character(0)
i <- 1
while (i <= length(raw_args)) {
  if (raw_args[i] == "--names" && i < length(raw_args)) {
    set_names_param <- strsplit(raw_args[i + 1], ",", fixed = TRUE)[[1]]
    i <- i + 2
    next
  }
  if (raw_args[i] == "--mode" && i < length(raw_args)) {
    mode_param <- raw_args[i + 1]
    i <- i + 2
    next
  }
  args <- c(args, raw_args[i])
  i <- i + 1
}

if (is.null(mode_param)) {
  stop("--mode is required. Use: distinct, intersect, or union")
}
mode <- match.arg(mode_param, c("distinct", "intersect", "union"))

if (length(args) < 3) {
  stop("Usage: Rscript peaks_ops_upsetR.R --mode MODE [--names \"N1,N2,...\"] OUTPUT_BASE BED1 BED2 [BED3 ...]")
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

# Count valid BED lines
count_bed_lines <- function(path) {
  if (!file.exists(path)) return(0)
  d <- tryCatch(read.table(path, sep = "\t", header = FALSE, comment.char = ""), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(0)
  sum(nchar(trimws(as.character(d[, 1]))) > 0 & d[, 1] != "")
}

# intersect mode: count = |intersect(sets)|
count_intersect <- function(idx) {
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

# union mode: count = |union(sets)|
count_union <- function(idx) {
  if (length(idx) == 1) {
    count_bed_lines(bed_files[idx])
  } else {
    tmp_merged <- tempfile(fileext = ".bed")
    files_to_cat <- paste(shQuote(bed_files[idx]), collapse = " ")
    exit_code <- system(paste0("cat ", files_to_cat, " | sort -k1,1 -k2,2n | bedtools merge -i stdin > ", shQuote(tmp_merged)))
    cnt <- 0
    if (exit_code == 0) {
      cnt <- count_bed_lines(tmp_merged)
    } else {
      warning("bedtools merge failed (exit ", exit_code, ") for combo: ",
              paste(set_names[idx], collapse = "&"))
    }
    unlink(tmp_merged, force = TRUE)
    cnt
  }
}

# distinct mode: count = |setdiff(intersect(in_sets), union(out_sets))|
count_distinct <- function(idx) {
  in_sets <- bed_files[idx]
  out_idx <- setdiff(seq_len(n), idx)
  # Step 1: intersect(in_sets)
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
  # Step 2: if out_sets non-empty, subtract union(out_sets)
  if (length(out_idx) == 0) {
    cnt <- count_bed_lines(tmp_intersect)
    unlink(tmp_intersect, force = TRUE)
    return(cnt)
  }
  out_sets <- bed_files[out_idx]
  tmp_union <- tempfile(fileext = ".bed")
  tmp_result <- tempfile(fileext = ".bed")
  files_cat <- paste(shQuote(out_sets), collapse = " ")
  exit_code <- system(paste0("cat ", files_cat, " | sort -k1,1 -k2,2n | bedtools merge -i stdin > ", shQuote(tmp_union)))
  if (exit_code != 0) {
    unlink(c(tmp_intersect, tmp_union, tmp_result), force = TRUE)
    return(0)
  }
  exit_code <- system(paste0("bedtools subtract -a ", shQuote(tmp_intersect), " -b ", shQuote(tmp_union), " > ", shQuote(tmp_result)))
  cnt <- 0
  if (exit_code == 0) {
    cnt <- count_bed_lines(tmp_result)
  }
  unlink(c(tmp_intersect, tmp_union, tmp_result), force = TRUE)
  cnt
}

# Select count function by mode (for bar heights)
count_combo <- switch(mode,
  intersect = count_intersect,
  union = count_union,
  distinct = count_distinct
)

# UpSetR fromExpression expects MUTUALLY EXCLUSIVE partition sizes. It computes
# set sizes (right bar) by summing all combo sizes that include that set.
# If we pass overlapping counts (intersect/union), the sum inflates the set size.
# Fix: always use partition (distinct) counts for fromExpression so set sizes are correct.
# Bar heights: distinct mode = partition (correct); intersect/union = partition (set sizes correct).
expr_for_upset <- count_distinct

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
  expr[i] <- expr_for_upset(combos[[i]])
}

# Skip if no multi-set overlap (all 2+ way are 0)
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

# Colors by degree
degree_colors <- c(
  "1" = "#7FCDBB", "2" = "#2E86AB", "3" = "#E94F37",
  "4" = "#44AF69", "5" = "#F18F01", "6" = "#5C4D7D"
)
sets_bar_col <- c("turquoise4")
matrix_col <- c("slateblue4")
shade_col <- c("wheat4")

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

text_scale <- c(1.5, 1.25, 1.25, 1, 1.7, 1.5)
upset_data <- UpSetR::fromExpression(expr)

pdf(paste0(out_base, ".pdf"), width = 8, height = 6)
UpSetR::upset(upset_data,
              mainbar.y.label = paste0("Peak ", mode, " (count)"),
              sets.x.label = "Peak count",
              text.scale = text_scale,
              sets.bar.color = sets_bar_col,
              matrix.color = matrix_col, shade.color = shade_col,
              queries = queries, query.legend = "none")
invisible(dev.off())
