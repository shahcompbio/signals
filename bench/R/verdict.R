#!/usr/bin/env Rscript
# Phase-free pairwise verdict: for pairs of cells that are both imbalanced over
# a chromosome arm, did they gain/retain the SAME parental haplotype or
# OPPOSITE ones?
#
#   Rscript bench/R/verdict.R <allele_counts> <cn_source> <out.rds>
#       [--cells FILE] [--select hash|lr|strat|file] [--maxcells N]
#
# <cn_source>  any file with cell_id, chr, start, state. Candidates are drawn
#              from ITS cells, so it must be QC-filtered: use a signals output
#              that passes both hard gates (bench/R/score.R), or pass --cells.
# --cells      restrict candidates to the cell_id column of FILE (e.g. a
#              signals output or qc.csv).
# --select     how to choose cells on an arm with more than --maxcells
#              candidates. Arms at or under the cap take every candidate.
#                strat  round-robin over bands of q (0-.1-.2-.3-.4-.5), then
#                       by hash within a band (default). Cells in different
#                       imbalance states tend to be different subclones, which
#                       is where opposite-parent pairs - the ones relabelling
#                       collapses - come from. On S2 it caught 48 relabelling
#                       failures against 2-3 for every other rule, and tied on
#                       S1 (all 8 known failures)
#                hash   stable pseudo-random choice by a hash of (arm, cell);
#                       proportional to clone size, so it under-samples small
#                       subclones
#                lr     strongest phase-free imbalance first
#                file   first N in file order - the original behaviour, kept
#                       only to reproduce verdicts cached before this change
# --maxcells   cap per arm (default 30); cost grows with its square
#
# The statistic marginalises the unknown block phase, so it uses no phasing at
# all and its answer cannot move when the phasing changes:
#
#   logLR = sum_h log[ .5 B(a_ih;n_ih,q_i) B(a_jh;n_jh,q_j)
#                    + .5 B(a_ih;n_ih,1-q_i) B(a_jh;n_jh,1-q_j) ]
#         - sum_h log[ .5 B(a_ih;n_ih,q_i) B(a_jh;n_jh,1-q_j)
#                    + .5 B(a_ih;n_ih,1-q_i) B(a_jh;n_jh,q_j) ]
#
# calibrated by permuting the block pairing (keeping both marginals); the
# statistic is z, and |z| >= 5 counts as determined. Outside file mode, each
# pair's permutations are seeded from the pair's own identity, so a pair's z
# does not depend on which other arms or pairs were tested.
#
# WHY QC MATTERS. QC-failed cells look non-diploid on every arm, so they pass
# the candidate filter everywhere, and they look strongly imbalanced too
# (on S1: median q 0.21 at a tenth of the depth of real candidates, versus
# 0.36) - so ranking by evidence does not screen them out. Built from the raw
# HMMcopy reads, the S1 arms where relabelling was shown to fail had 114-120
# candidates instead of 10-11, and the old first-30 cut kept almost none of the
# real ones: no determined pair survived on any of those arms.
#
# Known limitations:
#   - needs TWO imbalanced cells, so it is blind where only one carries an
#     event, and saturates once a method gets every determined pair right
#   - assumes one orientation per cell per arm; a genuine mid-arm switch
#     dilutes z toward 0 (a strong |z| is trustworthy, a weak one is not)
#   - 40 permutations estimate the null's spread noisily: reseeding flips 1.4%
#     of verdicts on S1 (55 of 4039 pairs, median |z| change 0.37), all with
#     |z| < 6. This does not bias method comparisons - every method is scored
#     against the same cached panel - but which borderline pairs are in the
#     panel is partly chance. Raise NPERM if that matters

local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(f), "common.R"))
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list(cells = NULL, select = "strat", maxcells = "30")
pos <- character(0); i <- 1
while (i <= length(args)) {
  if (args[i] %in% c("--cells", "--select", "--maxcells")) {
    opt[[sub("^--", "", args[i])]] <- args[i + 1]; i <- i + 2
  } else { pos <- c(pos, args[i]); i <- i + 1 }
}
if (length(pos) != 3) stop("usage: verdict.R <allele_counts> <cn_source> <out.rds> [options]")
counts <- pos[1]; cnsrc <- pos[2]; out <- pos[3]
SELECT <- match.arg(opt$select, c("strat", "hash", "lr", "file"))
MAXCELLS <- as.integer(opt$maxcells)
NPERM <- 40
ZDET <- 5
MINBLOCKS <- 30
set.seed(1)   # only governs file mode; other modes seed per pair

logLR <- function(ai, ni, aj, nj, qi, qj)
  sum(log(0.5 * dbinom(ai, ni, qi) * dbinom(aj, nj, qj) +
          0.5 * dbinom(ai, ni, 1 - qi) * dbinom(aj, nj, 1 - qj))) -
  sum(log(0.5 * dbinom(ai, ni, qi) * dbinom(aj, nj, 1 - qj) +
          0.5 * dbinom(ai, ni, 1 - qi) * dbinom(aj, nj, qj)))

# deterministic, platform-independent string hash (for selection and seeding)
hkey <- function(s) vapply(s, function(x) {
  h <- 0
  for (ch in utf8ToInt(x)) h <- (h * 31 + ch) %% 2147483647
  h
}, numeric(1), USE.NAMES = FALSE)

message("reading allele counts: ", counts)
b <- read_block_counts(counts)
message("reading copy number: ", cnsrc)
cn <- read_states(cnsrc)

# a cell is a candidate on an arm if its total copy number there admits
# imbalance at all: mostly odd, or mostly non-diploid
cand <- cn[, .(frac_odd = mean(state %% 2 == 1, na.rm = TRUE),
               frac_nd = mean(state != 2, na.rm = TRUE)), by = .(cell_id, unit)][
             frac_odd > 0.5 | frac_nd > 0.5]
if (!is.null(opt$cells)) {
  keep <- unique(fread(opt$cells, select = "cell_id")$cell_id)
  n0 <- uniqueN(cand$cell_id)
  cand <- cand[cell_id %in% keep]
  message(sprintf("--cells: kept %d of %d candidate cells", uniqueN(cand$cell_id), n0))
}

# phase-free evidence for every candidate (cheap: ~1 s for 2000 cell-arms)
ev <- b[cand[, .(cell_id, unit)], on = .(cell_id, unit), nomatch = 0][, {
  o <- optimize(function(q) -mix_ll(a0, n, q), c(0.02, 0.5))
  .(q = o$minimum, lr = -o$objective - mix_ll(a0, n, 0.5), nb = .N)
}, by = .(cell_id, unit)]

# sanity check on the candidate pool: QC-failed cells tend to have no allele
# counts at all, so a high share of count-free candidates means an unfiltered
# copy number source. Only arms that carry allele data count - the acrocentric
# p arms (13p 14p 15p 21p 22p) have copy number but no haplotype blocks for any
# cell, and would otherwise trip this on every sample.
withdata <- cand[unit %in% unique(b$unit)]
nocount <- 1 - nrow(ev) / max(nrow(withdata), 1)
# QC-passing sources measure 0%, raw HMMcopy reads 27% on S1; 10% sits in the gap
if (nocount > 0.1)
  warning(sprintf(paste("%.0f%% of candidate (cell, arm) pairs have no allele counts - the",
                        "copy number source looks unfiltered. Use a gated signals output or --cells."),
                  100 * nocount), call. = FALSE)

pick <- function(u) {
  if (SELECT == "file") {
    cells <- intersect(cand[unit == u]$cell_id, unique(b[unit == u]$cell_id))
    return(head(cells, MAXCELLS))
  }
  e <- ev[unit == u & nb >= MINBLOCKS]
  if (nrow(e) > MAXCELLS) {
    e[, h := hkey(paste(u, cell_id))]
    e <- switch(SELECT,
      hash  = e[order(h)][seq_len(MAXCELLS)],
      lr    = e[order(-lr, cell_id)][seq_len(MAXCELLS)],
      strat = {
        e[, band := cut(q, c(0, .1, .2, .3, .4, .5), include.lowest = TRUE)]
        setorder(e, band, h)
        e[, r := seq_len(.N), by = band]
        e[order(r, band)][seq_len(MAXCELLS)]
      })
  }
  sort(e$cell_id)
}

units <- if (SELECT == "file") unique(cand$unit) else sort(unique(cand$unit))
res <- rbindlist(lapply(units, function(u) {
  cells <- pick(u)
  if (length(cells) < 2) return(NULL)
  bu <- b[unit == u & cell_id %in% cells]
  qs <- bu[, .(q = fit_q(a0, n)), by = cell_id]
  rbindlist(lapply(combn(cells, 2, simplify = FALSE), function(pr) {
    mi <- bu[cell_id == pr[1], .(start, hap_label, a0, n)]
    mj <- bu[cell_id == pr[2], .(start, hap_label, a0, n)]
    m <- merge(mi, mj, by = c("start", "hap_label"), suffixes = c(".i", ".j"))
    if (nrow(m) < MINBLOCKS) return(NULL)
    qi <- qs[cell_id == pr[1]]$q; qj <- qs[cell_id == pr[2]]$q
    obs <- logLR(m$a0.i, m$n.i, m$a0.j, m$n.j, qi, qj)
    if (SELECT != "file") set.seed(hkey(paste(u, pr[1], pr[2])))
    nul <- replicate(NPERM, {
      p <- sample.int(nrow(m))
      logLR(m$a0.i, m$n.i, m$a0.j[p], m$n.j[p], qi, qj)
    })
    data.table(unit = u, cell_i = pr[1], cell_j = pr[2], nblocks = nrow(m),
               q_i = qi, q_j = qj, z = (obs - mean(nul)) / sd(nul))
  }))
}), fill = TRUE)

res[, verdict := fifelse(z >= ZDET, "same parent",
                 fifelse(z <= -ZDET, "opposite parents", "undetermined"))]
attr(res, "inputs") <- list(counts = normalizePath(counts), cn_source = normalizePath(cnsrc),
                            cells = opt$cells, select = SELECT, maxcells = MAXCELLS,
                            nperm = NPERM, zdet = ZDET, created = format(Sys.time()))
saveRDS(res, out)
det <- res[verdict != "undetermined"]
cat(sprintf("select=%s maxcells=%d: %d pairs, %d determined across %d arms; %d opposite-parent\n",
            SELECT, MAXCELLS, nrow(res), nrow(det), uniqueN(det$unit),
            sum(det$verdict == "opposite parents")))
