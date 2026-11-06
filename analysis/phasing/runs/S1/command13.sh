#!/usr/bin/env bash
# As run 12, but the count-free-bin inheritance now runs AFTER the
# minimum-breakpoint DP instead of before it.
#
# Run 12 improved isolated bins (367 -> 297) but left chr13 at 2.43, unchanged
# from run 11, because the inheritance ran first and the DP then re-oriented the
# free 13q segment, stranding the pinned 13p bins it had just been matched to.
# Only 4 cell-chromosome pairs changed at all.
#
# (My chr13 diagnosis in run 11 was also wrong: I cut the arm at 19.5Mb when
# hg19 chr13's centromere is ~18Mb, so what looked like an internal 13p flip was
# really the 13p/13q junction. 13p is 36 bins, all count-free.)
#
# Verified offline against run 12's own output: chr13 2.43 -> 2.33, chr14
# unchanged at 1.73. That would make every chromosome neutral-or-better than the
# run 5 baseline.
#
#   run 13 vs run 12 -> the reordering alone
#   run 13 vs run 5  -> the whole chain against the baseline

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_13.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_13.png \
--heatmap heatmap_13.png \
--heatmapraw heatmapraw_13.png \
--csvfile hscn_13.csv.gz \
--qccsvfile qc_13.csv \
--Rdatafile signals_13.Rdata \
--mincells 2 \
--sex FEMALE
