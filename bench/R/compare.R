#!/usr/bin/env Rscript
# What moved between two signals outputs, and where.
#
#   Rscript bench/R/compare.R <baseline.csv.gz> <run.csv.gz> [--top N] [--out f.csv]
#
# Reports total and allele-specific changes, then the arms and cells that moved
# most. Total copy number should never differ between two runs on the same input
# - if it does, that is reported first, since it means the phasing has leaked
# into the copy number call.

local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source(file.path(dirname(f), "common.R"))
})

args <- commandArgs(trailingOnly = TRUE)
top <- 10; out <- NULL; pos <- character(0); i <- 1
while (i <= length(args)) {
  if (args[i] == "--top") { top <- as.integer(args[i + 1]); i <- i + 2 }
  else if (args[i] == "--out") { out <- args[i + 1]; i <- i + 2 }
  else { pos <- c(pos, args[i]); i <- i + 1 }
}
if (length(pos) != 2) stop("usage: compare.R <baseline> <run> [--top N] [--out f.csv]")

x <- read_hscn(pos[1]); y <- read_hscn(pos[2])
k <- c("cell_id", "chr", "start")
m <- merge(x[, c(k, "unit", "state", "A", "B"), with = FALSE],
           y[, c(k, "state", "A", "B"), with = FALSE],
           by = k, all = TRUE, suffixes = c(".base", ".run"))

only_base <- sum(is.na(m$state.run) & !is.na(m$state.base))
only_run <- sum(is.na(m$state.base) & !is.na(m$state.run))
m <- m[!is.na(state.base) & !is.na(state.run)]
m[, state_diff := state.base != state.run]
m[, as_diff := fifelse(is.na(A.base) | is.na(A.run), is.na(A.base) != is.na(A.run),
                       A.base != A.run | B.base != B.run)]

cat(sprintf("baseline: %s\nrun:      %s\n\n", pos[1], pos[2]))
cat(sprintf("bins compared  %d   (only in baseline %d, only in run %d)\n",
            nrow(m), only_base, only_run))
cat(sprintf("state changed  %d%s\n", sum(m$state_diff),
            if (sum(m$state_diff)) "   <- total copy number moved; phasing has leaked into it" else ""))
cat(sprintf("A|B changed    %d   (%.3f%%)  in %d cells\n",
            sum(m$as_diff), 100 * mean(m$as_diff), uniqueN(m[as_diff == TRUE]$cell_id)))

if (any(m$as_diff)) {
  # a swapped orientation (A|B -> B|A) is a phasing change; anything else
  # changed the allele-specific state itself
  m[, swap := as_diff & !is.na(A.base) & !is.na(A.run) & A.base == B.run & B.base == A.run]
  cat(sprintf("  of which pure orientation swaps %d, state changes %d\n",
              sum(m$swap), sum(m$as_diff & !m$swap)))

  by_arm <- m[, .(bins = .N, changed = sum(as_diff), swaps = sum(swap),
                  cells = uniqueN(cell_id[as_diff])), by = unit][changed > 0][order(-changed)]
  by_arm[, pct := round(100 * changed / bins, 2)]
  cat(sprintf("\narms with most change (top %d of %d):\n", min(top, nrow(by_arm)), nrow(by_arm)))
  print(head(by_arm, top))

  by_cell <- m[as_diff == TRUE, .(changed = .N, arms = paste(sort(unique(unit)), collapse = " ")),
               by = cell_id][order(-changed)]
  by_cell[, cell := short_id(cell_id)]
  cat(sprintf("\ncells with most change (top %d of %d):\n", min(top, nrow(by_cell)), nrow(by_cell)))
  print(head(by_cell[, .(cell, changed, arms)], top))
  if (!is.null(out)) fwrite(by_arm, out)
}
