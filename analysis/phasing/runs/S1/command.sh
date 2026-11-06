#!/usr/bin/env bash
# Outputs go to *_relabel.* so the previous run is kept as a baseline to diff against.
# Drop the suffixes once you are happy with the result.

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_relabel.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_relabel.png \
--heatmap heatmap_relabel.png \
--heatmapraw heatmapraw_relabel.png \
--csvfile hscn_relabel.csv.gz \
--qccsvfile qc_relabel.csv \
--Rdatafile signals_relabel.Rdata \
--mincells 2 \
--sex FEMALE
