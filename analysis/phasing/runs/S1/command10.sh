#!/usr/bin/env bash
# Arm-level detection, oriented pooling, and unit-aware relabelling.
#
# Run 9 fixed both target cells (chr4q R59-C30: 12 -> 1 run; chr2q R57-C12
# restored to 2 runs) but was still worse than run 5 overall - 13,213 vs 12,150
# segments - and it was all chromosome 17: mean 4.62 -> 6.76, every other
# chromosome neutral or better.
#
# Cause was NOT phasing: chr17's block phases are byte-identical between runs 5
# and 9. It was relabel_unsupported_phase. Its rule was "a chromosome is free
# only if no other cell supports ANY unit on it", so 17p having 2 candidates
# blocked relabelling of 17q, which has 0 and needs it. chr17 was relabelled in
# 118,175 rows in run 5 and 0 in run 9.
#
# The minimum-breakpoint DP now takes a per-bin `free` flag: supported units
# stay pinned, unsupported ones are still cleaned up, and the cross-boundary
# optimisation is kept.
#
#   run 10 vs run 9  -> the relabelling fix alone
#   run 10 vs run 5  -> everything since, together

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_10.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_10.png \
--heatmap heatmap_10.png \
--heatmapraw heatmapraw_10.png \
--csvfile hscn_10.csv.gz \
--qccsvfile qc_10.csv \
--Rdatafile signals_10.Rdata \
--mincells 2 \
--sex FEMALE
