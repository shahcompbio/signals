#!/usr/bin/env Rscript
# Phase-free pairwise verdict: for pairs of cells that are both imbalanced over
# a chromosome arm, did they gain/retain the SAME parental haplotype or
# OPPOSITE ones?
#
#   Rscript bench/R/verdict.R <allele_counts.csv.gz> <cn_source> <out.rds>
#
# <cn_source> is any file with cell_id, chr, start, state. Use the output of a
# signals run that passes both hard gates (bench/R/score.R): its cells are the
# QC-passing set and its state equals the HMMcopy input, so it is the same for
# every method. Do NOT use the raw HMMcopy reads - see the warning below.
#
# The statistic marginalises the unknown block phase, so it uses no phasing at
# all and its answer cannot move when the phasing changes:
#
#   logLR = sum_h log[ .5 B(a_ih;n_ih,q_i) B(a_jh;n_jh,q_j)
#                    + .5 B(a_ih;n_ih,1-q_i) B(a_jh;n_jh,1-q_j) ]
#         - sum_h log[ .5 B(a_ih;n_ih,q_i) B(a_jh;n_jh,1-q_j)
#                    + .5 B(a_ih;n_ih,1-q_i) B(a_jh;n_jh,q_j) ]
#
# logLR scales with block count and depth, so it is calibrated by permuting the
# block pairing (keeping both marginals); the reported statistic is z, and
# |z| >= 5 counts as determined.
#
# Known limitations:
#   - needs TWO imbalanced cells, so it is blind where only one cell carries an
#     event, and saturates once a method gets every determined pair right
#   - assumes one orientation per cell per arm; a genuine mid-arm switch dilutes
#     z toward 0 (a strong |z| is trustworthy, a weak one is not evidence)
#   - takes the first 30 candidate cells per arm, IN FILE ORDER. This makes the
#     set of pairs tested an accident of row order, and it is the weakest part
#     of the apparatus. Measured on S1: building from the raw HMMcopy reads
#     instead of a signals output (920 cells vs 726 QC-passing, interleaved)
#     changed the determined pairs from 1751 to 1167, with only 1087 in common,
#     and the opposite-parent pairs from 78 to 4. Not one of the five arms where
#     relabelling was shown to fail (1p 2q 5q 9q 18p) kept a determined pair -
#     so that defect would have been invisible. The statistic itself is stable:
#     zero verdicts flip among pairs both versions determine. See PLAN.md for
#     the fix (evidence-ranked, QC-restricted selection).

local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(f), "common.R"))
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: verdict.R <allele_counts> <cn_source> <out.rds>")
counts <- args[1]; cnsrc <- args[2]; out <- args[3]

NPERM <- 40
ZDET <- 5
MAXCELLS <- 30
MINBLOCKS <- 30
set.seed(1)

logLR <- function(ai, ni, aj, nj, qi, qj)
  sum(log(0.5 * dbinom(ai, ni, qi) * dbinom(aj, nj, qj) +
          0.5 * dbinom(ai, ni, 1 - qi) * dbinom(aj, nj, 1 - qj))) -
  sum(log(0.5 * dbinom(ai, ni, qi) * dbinom(aj, nj, 1 - qj) +
          0.5 * dbinom(ai, ni, 1 - qi) * dbinom(aj, nj, qj)))

message("reading allele counts: ", counts)
b <- read_block_counts(counts)
message("reading copy number: ", cnsrc)
cn <- read_states(cnsrc)

# a cell is a candidate on an arm if its total copy number there admits
# imbalance at all: mostly odd, or mostly non-diploid
cand <- cn[, .(frac_odd = mean(state %% 2 == 1, na.rm = TRUE),
               frac_nd = mean(state != 2, na.rm = TRUE)), by = .(cell_id, unit)][
             frac_odd > 0.5 | frac_nd > 0.5]

res <- rbindlist(lapply(unique(cand$unit), function(u) {
  cells <- cand[unit == u]$cell_id
  bu <- b[unit == u & cell_id %in% cells]
  cells <- intersect(cells, unique(bu$cell_id))
  if (length(cells) < 2) return(NULL)
  if (length(cells) > MAXCELLS) cells <- cells[seq_len(MAXCELLS)]
  bu <- bu[cell_id %in% cells]
  qs <- bu[, .(q = fit_q(a0, n)), by = cell_id]
  rbindlist(lapply(combn(cells, 2, simplify = FALSE), function(pr) {
    mi <- bu[cell_id == pr[1], .(start, hap_label, a0, n)]
    mj <- bu[cell_id == pr[2], .(start, hap_label, a0, n)]
    m <- merge(mi, mj, by = c("start", "hap_label"), suffixes = c(".i", ".j"))
    if (nrow(m) < MINBLOCKS) return(NULL)
    qi <- qs[cell_id == pr[1]]$q; qj <- qs[cell_id == pr[2]]$q
    obs <- logLR(m$a0.i, m$n.i, m$a0.j, m$n.j, qi, qj)
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
                            nperm = NPERM, zdet = ZDET, maxcells = MAXCELLS,
                            created = format(Sys.time()))
saveRDS(res, out)
cat("pairs:", nrow(res), "  determined:", sum(res$verdict != "undetermined"),
    "across", uniqueN(res[verdict != "undetermined"]$unit), "units\n")
print(res[verdict != "undetermined", .N, by = verdict])
