# Shared helpers for the signals benchmark. data.table only, so every script
# runs natively rather than inside the (x86, emulated) signals image.

suppressMessages(library(data.table))
setDTthreads(0)

BIN <- 5e5

# Directory of the bench/ tree, from wherever the calling script lives.
bench_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath("bench")
}

# Chromosome arm, using the same rule as signals::coord_to_arm(): a position is
# on p while it is below the end of the last p band in signals' cytoband_map.
# The boundaries are extracted from that map (see data/arm_boundaries_hg19.csv)
# rather than hand-coded - an earlier hand-coded centromere for chr13 produced a
# wrong conclusion, and the scorer must agree with the method it scores.
.arm_bounds <- NULL
arm_of <- function(chr, start) {
  if (is.null(.arm_bounds)) {
    b <- fread(file.path(bench_dir(), "data", "arm_boundaries_hg19.csv"),
               colClasses = list(character = "chr"))
    .arm_bounds <<- setNames(b$p_end, b$chr)
  }
  ifelse(start < .arm_bounds[as.character(chr)], "p", "q")
}
unit_of <- function(chr, start) paste0(chr, arm_of(chr, start))

# Allele counts in the DLP haplotype_calling_count format, summed to
# (cell, chr, 500kb bin, block). a0/a1 are the two alleles as called, with no
# phasing applied.
read_block_counts <- function(path, chrs = NULL, mincounts = 2) {
  h <- fread(path, colClasses = list(character = c("chromosome", "cell_id")))
  setnames(h, "chromosome", "chr")
  if (!is.null(chrs)) h <- h[chr %in% chrs]
  h[, allele_id := paste0("a", allele_id)]
  w <- dcast(h, cell_id + chr + start + hap_label ~ allele_id,
             value.var = "readcount", fun.aggregate = sum, fill = 0L)
  rm(h); gc(verbose = FALSE)
  w[, start := floor(start / BIN) * BIN + 1]
  b <- w[, .(a0 = sum(a0), a1 = sum(a1)), by = .(cell_id, chr, start, hap_label)]
  b[, n := a0 + a1]
  b <- b[n >= mincounts]
  b[, unit := unit_of(chr, start)]
  b[]
}

# Per-bin total copy number from any file carrying cell_id, chr, start, state -
# an HMMcopy reads file or a signals output.
read_states <- function(path) {
  d <- fread(path, select = c("cell_id", "chr", "start", "state"),
             colClasses = list(character = c("chr", "cell_id")))
  d[, unit := unit_of(chr, start)]
  d[]
}

read_hscn <- function(path, cols = c("cell_id", "chr", "start", "state", "A", "B")) {
  d <- fread(path, select = cols, colClasses = list(character = c("chr", "cell_id")))
  setorder(d, cell_id, chr, start)
  d[, unit := unit_of(chr, start)]
  d[]
}

# Phase-free symmetric mixture: the allelic ratio q <= 0.5 a cell shows over a
# set of blocks without knowing which allele is which, and the log likelihood
# ratio against q = 0.5.
mix_ll <- function(a, n, q) sum(log(0.5 * dbinom(a, n, q) + 0.5 * dbinom(a, n, 1 - q)))
fit_q <- function(a, n) optimize(function(q) -mix_ll(a, n, q), c(0.02, 0.5))$minimum

# Short, stable label for a cell id, for printing.
short_id <- function(x) sub(".*-(R[0-9]+-C[0-9]+)$", "\\1", x)
