#!/usr/bin/env bash
# Arm-level DETECTION with a chromosome-level first pass.
#
# Run 7 set --phasebyarm, which made the first pass cluster per arm as well.
# That was worse: an arm has far fewer bins for UMAP, so the first-pass cell
# choice degraded badly (best_used_rank 99-227 across many arms) and total
# segments went 12,150 -> 13,639 against run 5.
#
# Run 7 also never applied any correction: the line extracting flagged units
# read the wrong column, so the second pass was skipped entirely. Fixed.
#
# This run keeps the first pass exactly as run 5 and changes only the unit used
# to detect and correct bad phasing. Expected: 4q re-phased from R59-C30, 4p
# left alone, plus 1p and 18p - three units, not 21.
#
#   run 8 vs run 5  -> the effect of arm-level detection alone

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_8.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_8.png \
--heatmap heatmap_8.png \
--heatmapraw heatmapraw_8.png \
--csvfile hscn_8.csv.gz \
--qccsvfile qc_8.csv \
--Rdatafile signals_8.Rdata \
--mincells 2 \
--sex FEMALE
