#!/usr/bin/env bash
# Arm-level detection, with candidate phasing cells put into a common
# orientation before their counts are pooled.
#
# Run 8 fixed chromosome 4 but regressed 2q for R57-C12 (2 runs -> 12). Cause:
# the three cells pooled to phase 2q carried opposite parents (pairwise block
# agreement 78%, 37%, 24%), so summing their raw counts made them cancel and the
# pooled majority became a patchwork. orient_cells() now flips cells that
# disagree with the highest-scoring one before pooling.
#
# Effect on 2q for R57-C12, at bin level: 134/154 bins either side of BAF 0.5
# (alternating) -> 267/25, mean BAF 0.656, i.e. a clean 1|2.
# Narrow: 4 of 21 flagged units are affected (1p, 2q, 14q, 16p). 4q has a single
# candidate and is untouched, so the chromosome 4 fix carries over.
#
#   run 9 vs run 8  -> the orientation fix alone
#   run 9 vs run 5  -> arm-level detection + orientation, together

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_9.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_9.png \
--heatmap heatmap_9.png \
--heatmapraw heatmapraw_9.png \
--csvfile hscn_9.csv.gz \
--qccsvfile qc_9.csv \
--Rdatafile signals_9.Rdata \
--mincells 2 \
--sex FEMALE
