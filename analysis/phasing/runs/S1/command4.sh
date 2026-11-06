#!/usr/bin/env bash
# Control run: seeded, but WITHOUT --rephase_by_imbalance and
# --relabel_unsupported_phase.
#
# The point is to separate three effects that are currently confounded. The
# pre-fix baseline (hscn.csv.gz) was unseeded and had a different fitted rho
# (0.02831 vs 0.02777), so diffing it against run 2 mixes the seed change, the
# re-phasing and the relabelling together. This run isolates them:
#
#   run 4 vs pre-fix baseline  -> effect of seeding alone
#   run 2 vs run 4             -> effect of re-phasing + relabelling alone
#
# Same seed as runs 2 and 3 so the first pass is identical to theirs.

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_4.png \
--heatmap heatmap_4.png \
--heatmapraw heatmapraw_4.png \
--csvfile hscn_4.csv.gz \
--qccsvfile qc_4.csv \
--Rdatafile signals_4.Rdata \
--mincells 2 \
--sex FEMALE
