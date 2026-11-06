#!/usr/bin/env bash
# Viterbi fix, OFF (legacy decode). The control arm for command5.sh: same
# image, same seed, same flags, only the decode differs.
#
#   run 5 vs run 6  -> the Viterbi fix in isolation
#   run 6 vs run 2  -> should be ~identical; validates that cpp_legacy really
#                      does reproduce signals 0.16.0
#
# REQUIRES a rebuilt image (signals >= 0.17.0) - --viterbiver does not exist in
# 0.16.0 and the run will fail with an unused-argument error.

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp_legacy \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_6.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_6.png \
--heatmap heatmap_6.png \
--heatmapraw heatmapraw_6.png \
--csvfile hscn_6.csv.gz \
--qccsvfile qc_6.csv \
--Rdatafile signals_6.Rdata \
--mincells 2 \
--sex FEMALE
