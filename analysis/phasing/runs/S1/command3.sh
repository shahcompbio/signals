#!/usr/bin/env bash
# Outputs go to *_3.* so the previous run is kept as a baseline to diff against.
# Drop the suffixes once you are happy with the result.

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_3.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_3.png \
--heatmap heatmap_3.png \
--heatmapraw heatmapraw_3.png \
--csvfile hscn_3.csv.gz \
--qccsvfile qc_3.csv \
--Rdatafile signals_3.Rdata \
--mincells 2 \
--sex FEMALE
