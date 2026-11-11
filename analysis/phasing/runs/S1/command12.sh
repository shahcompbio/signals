#!/usr/bin/env bash
# As run 11, plus: a run of bins with no haplotype counts inherits its A/B
# labels from the nearest informative bin.
#
# Run 11 is already the best result so far - 11,753 segments vs the run 5
# baseline's 12,150, 367 vs 399 isolated bins, 17 chromosomes neutral, 5 better,
# 1 marginally worse. This is polish on that last one.
#
# chr13 sat 0.09 above baseline because pinning the count-free 13p region
# preserved an internal flip the HMM had emitted there (36 bins 0|1 then 3 bins
# 1|0) which per-bin mirroring used to mask. Those bins have no data and their
# A/B came from fillmissing anyway, so they now inherit from the nearest
# informative bin - consistent internally and with the segment they adjoin.
#
# Verified offline against run 11's output: chr13 2.43 -> 2.33, chr14 unchanged
# at 1.73. Expect ~11,680 segments total.
#
#   run 12 vs run 11 -> should touch only chr13/14/15/22
#   run 12 vs run 5  -> the whole chain against the baseline

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_12.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_12.png \
--heatmap heatmap_12.png \
--heatmapraw heatmapraw_12.png \
--csvfile hscn_12.csv.gz \
--qccsvfile qc_12.csv \
--Rdatafile signals_12.Rdata \
--mincells 2 \
--sex FEMALE
