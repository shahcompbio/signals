#!/usr/bin/env bash
# As run 10, plus: units with no haplotype data at all are left alone.
#
# Run 10 beat the run 5 baseline overall (11,975 vs 12,150 segments, 367 vs 399
# isolated bins) and restored chr17 (6.76 -> 4.64). One chromosome regressed:
# chr14, 1.95 -> 2.14.
#
# Cause: 14p has 36 bins in the CN table and ZERO haplotype counts, so its
# support count is 0 - which the relabelling read as "unsupported, free to
# relabel" when it actually means "no data". Relabelling those count-free bins
# put a label break at the 14p/14q junction in every cell. The same applies to
# 13p, 15p and 22p.
#
# A unit is now only eligible for relabelling if it appears in the imbalance
# table at all, i.e. if it has haplotype data.
#
#   run 11 vs run 10 -> the no-data-unit fix alone (expect chr13/14/15/22 only)
#   run 11 vs run 5  -> the whole arm-level chain against the baseline

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_11.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_11.png \
--heatmap heatmap_11.png \
--heatmapraw heatmapraw_11.png \
--csvfile hscn_11.csv.gz \
--qccsvfile qc_11.csv \
--Rdatafile signals_11.Rdata \
--mincells 2 \
--sex FEMALE
