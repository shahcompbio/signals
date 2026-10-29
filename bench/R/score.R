#!/usr/bin/env Rscript
# Score signals outputs on the real-data panel.
#
#   Rscript bench/R/score.R --verdict v.rds [--input hmmcopy_reads.csv.gz]
#       [--out metrics.csv] [--detail dir/] [label=]hscn.csv.gz ...
#
# One row per run. Columns, grouped by role:
#
#   hard gates   gate_AB_ne_state    bins where A + B != state      (must be 0)
#                gate_state_changed  bins whose state differs from
#                                    the HMMcopy input (--input)     (must be 0)
#   correctness  pairs_*             agreement with the phase-free
#                                    pairwise verdict, overall and
#                                    split by same / opposite parent
#   diagnostic   purity_*            within a (cell, arm) with >= 30
#                                    imbalanced bins, the share of
#                                    bins voting the majority
#                                    orientation
#                segments, isolated_bins, bkpt_*, flat_*
#
# None of the diagnostics are correctness metrics. All are functions of the A/B
# labelling, so any method that relabels can improve them without being more
# right:
#   - segments: on S1 the most correct configuration had the most segments
#   - purity: relabelling mirrors everything to A >= B, making orientation
#     uniform by construction. On S2 the relabelled run scores 2.5% of cell-arms
#     below 0.9 purity against ~13.5% for every other run, while being the
#     least correct of them (opposite-parent pairs 3/5 vs 5/5).
# Compare diagnostics only between configurations that do not relabel.
#
# --detail writes the disagreeing pairs and per-arm purity for each run.

local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(f), "common.R"))
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list(verdict = NULL, input = NULL, out = NULL, detail = NULL)
runs <- character(0)
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a %in% c("--verdict", "--input", "--out", "--detail")) {
    opt[[sub("^--", "", a)]] <- args[i + 1]; i <- i + 2
  } else { runs <- c(runs, a); i <- i + 1 }
}
if (is.null(opt$verdict) || !length(runs))
  stop("usage: score.R --verdict v.rds [--input reads] [--out f] [--detail dir] run ...")

lab <- ifelse(grepl("=", runs), sub("=.*", "", runs), sub("\\.csv(\\.gz)?$", "", basename(runs)))
paths <- sub("^[^=]*=", "", runs)

V <- readRDS(opt$verdict)[verdict != "undetermined"]
message(sprintf("verdict: %d determined pairs across %d arms", nrow(V), uniqueN(V$unit)))
inp <- if (!is.null(opt$input)) {
  message("reading input states: ", opt$input)
  fread(opt$input, select = c("cell_id", "chr", "start", "state"),
        colClasses = list(character = c("chr", "cell_id")))
} else NULL
if (!is.null(opt$detail)) dir.create(opt$detail, showWarnings = FALSE, recursive = TRUE)

score_one <- function(label, path) {
  message("scoring ", label)
  d <- read_hscn(path)
  r <- list(run = label, bins = nrow(d))

  # --- hard gates -------------------------------------------------------------
  ok <- !is.na(d$A) & !is.na(d$B) & !is.na(d$state)
  r$gate_AB_ne_state <- sum(ok & d$A + d$B != d$state)
  r$gate_state_changed <- if (is.null(inp)) NA_integer_ else {
    m <- inp[d[, .(cell_id, chr, start, out_state = state)], on = .(cell_id, chr, start), nomatch = 0]
    sum(m$state != m$out_state, na.rm = TRUE)
  }

  # --- correctness: pairwise parent-of-origin ---------------------------------
  # a cell's orientation on an arm: sign of the majority of its imbalanced bins
  o <- d[ok & A != B, .(sgn = sign(sum(sign(A - B)))), by = .(cell_id, unit)]
  vi <- o[V[, .(cell_id = cell_i, unit)], on = .(cell_id, unit)]$sgn
  vj <- o[V[, .(cell_id = cell_j, unit)], on = .(cell_id, unit)]$sgn
  call <- fifelse(is.na(vi) | is.na(vj) | vi == 0 | vj == 0, NA_character_,
                  fifelse(vi == vj, "same parent", "opposite parents"))
  sc <- !is.na(call)
  agree <- sc & call == V$verdict
  r$pairs_scored <- sum(sc); r$pairs_agree <- sum(agree)
  for (v in c("same parent", "opposite parents")) {
    k <- if (v == "same parent") "same" else "opp"
    r[[paste0("pairs_", k, "_scored")]] <- sum(sc & V$verdict == v)
    r[[paste0("pairs_", k, "_agree")]] <- sum(agree & V$verdict == v)
  }

  # --- consistency: sgn purity ------------------------------------------------
  pu <- d[ok & A != B, .(nb = .N, nplus = sum(A > B), nminus = sum(A < B)), by = .(cell_id, unit)][nb >= 30]
  pu[, purity := pmax(nplus, nminus) / nb]
  r$purity_cellarms <- nrow(pu)
  r$purity_median <- if (nrow(pu)) median(pu$purity) else NA_real_
  r$purity_frac_below_0.9 <- if (nrow(pu)) mean(pu$purity < 0.9) else NA_real_

  # --- diagnostics: segments, isolated bins, transitions ----------------------
  d[, st := paste0(A, "|", B)]
  same <- d$cell_id == shift(d$cell_id) & d$chr == shift(d$chr)
  same[1] <- FALSE
  # a new segment starts at every cell/chromosome boundary and every A|B change
  r$segments <- sum(!same | d$st != shift(d$st), na.rm = TRUE)
  nxt_same <- shift(same, type = "lead"); nxt_same[is.na(nxt_same)] <- FALSE
  prv <- shift(d$st); nxt <- shift(d$st, type = "lead")
  r$isolated_bins <- sum(same & nxt_same & prv == nxt & d$st != prv, na.rm = TRUE)

  t <- data.table(dA = d$A - shift(d$A), dB = d$B - shift(d$B), ds = d$state - shift(d$state))
  t <- t[same & ok & shift(ok, fill = FALSE)]
  t[, nchanged := (dA != 0) + (dB != 0)]
  bk <- t[ds != 0]; fl <- t[ds == 0]
  r$bkpt_transitions <- nrow(bk)
  r$bkpt_two_allele <- sum(bk$nchanged == 2)
  r$bkpt_two_allele_pct <- if (nrow(bk)) 100 * mean(bk$nchanged == 2) else NA_real_
  r$flat_transitions <- nrow(fl)
  r$flat_changes <- sum(fl$nchanged > 0)

  if (!is.null(opt$detail)) {
    fails <- V[sc & !agree][, call := call[sc & !agree]]
    fwrite(fails[order(-abs(z))], file.path(opt$detail, paste0(label, ".disagreements.csv")))
    fwrite(pu[, .(cell_arms = .N, median_purity = median(purity),
                  frac_below_0.9 = mean(purity < 0.9)), by = unit][order(unit)],
           file.path(opt$detail, paste0(label, ".purity_by_arm.csv")))
  }
  as.data.table(r)
}

res <- rbindlist(Map(score_one, lab, paths))
res[, pairs_pct := round(100 * pairs_agree / pairs_scored, 2)]
res[, pairs_opp_pct := round(100 * pairs_opp_agree / pairs_opp_scored, 1)]
if (!is.null(opt$out)) fwrite(res, opt$out)

show <- res[, .(run, gate_AB_ne_state, gate_state_changed,
                pairs = paste0(pairs_agree, "/", pairs_scored),
                opp = paste0(pairs_opp_agree, "/", pairs_opp_scored),
                purity_med = round(purity_median, 3),
                purity_lt0.9 = round(purity_frac_below_0.9, 3),
                segments, isolated_bins,
                bkpt_2allele = paste0(bkpt_two_allele, " (", round(bkpt_two_allele_pct, 1), "%)"),
                flat_changes)]
print(show)
gate_fail <- res[gate_AB_ne_state > 0 | (!is.na(gate_state_changed) & gate_state_changed > 0)]$run
if (length(gate_fail)) message("\nHARD GATE FAILED: ", paste(gate_fail, collapse = ", "))
