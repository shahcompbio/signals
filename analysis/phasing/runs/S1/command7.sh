#!/usr/bin/env bash
# Arm-level phasing, ON. Identical to command5.sh apart from --phasebyarm and
# the output names, so run 7 vs run 5 isolates the arm change.
#
# Motivating case: chromosome 4. 17 cells carry p-arm LOH and correctly phase
# the p arm; those cells are balanced on the q arm, so q-arm blocks are phased
# off noise, and R59-C30 (total CN 3 throughout, hence necessarily imbalanced)
# alternates 1|2 / 2|1 across 12 runs on q while reading a clean 2|1 on p.
#
# --phasebyarm switches BOTH the phasing unit and the imbalance detection unit.
# Expected to re-phase 4q from R59-C30, leave 4p alone (its cells are already
# the best available), and pick up 1p and 18p as well.
#
# REQUIRES A REBUILT IMAGE, and this is not self-evident: run_signals.R is
# bind-mounted from the host, so --phasebyarm is accepted immediately, but
# arm-level phasing itself is a package change baked into the image. Against an
# older installed signals the flag is honoured by the script and ignored by the
# package, the cell list is keyed by chromosome while the phasing filters on
# arm, nothing matches, and the run dies later with
#   Error in fifelse(state_min < 0, 0, state_min):
#     'no' is of type NULL but 'yes' is double
# run_signals.R now checks for this up front and stops with a clear message.

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --phasebyarm \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_7.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_7.png \
--heatmap heatmap_7.png \
--heatmapraw heatmapraw_7.png \
--csvfile hscn_7.csv.gz \
--qccsvfile qc_7.csv \
--Rdatafile signals_7.Rdata \
--mincells 2 \
--sex FEMALE
