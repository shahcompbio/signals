#!/usr/bin/env bash
# Joint cross-cell optimisation (C2), WITHOUT the relabelling step.
#
# Why relabelling is off
# ----------------------
# Measured on the well-characterised sample against the phase-free pairwise
# verdict (126 well-determined pairs, |z| >= 5):
#
#   run 13  arm re-phasing + relabelling      121/126   96%
#   run 17  C2 + relabelling                  117/126   93%
#   run 14  C2, no relabelling                126/126  100%
#
# split by verdict:
#                         pairs   r13    r14    r17
#   same parent             105    95%   100%   100%
#   opposite parents         21   100%   100%    57%
#
# Relabelling mirrors undetermined bins to A >= B, which makes two cells that
# gained OPPOSITE parents look identical. That is precisely the failure seen on
# this sample's chr15/19, where the phase-free test says opposite (z = -8.6,
# -9.2) and the pipeline reported the same parent.
#
# Note that run 14 also had the MOST segments of the three. Segment count
# rewards aggressive relabelling and says nothing about whether the parent-of-
# origin call is right; it should not be used as the quality metric here.
#
# What to expect
# --------------
#  * chr15 and chr19: the two gain cells should come out on OPPOSITE sides.
#    Check with:  Rscript phase-confidence.R <allele_counts> hscn_c2.csv.gz 15q
#  * Chromosomes with no phasing information will look NOISIER on the heatmap
#    than before. That is honest rather than wrong - relabelling was hiding it.
#    confidence_c2.csv says which chromosomes those are: low frac_determined.
#  * ~4x runtime, since each joint iteration is a full HSCN inference.

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --joint_phasing \
  --joint_iters 4 \
--hmmcopyqc SHAH_H003831_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H003831_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts SHAH_H003831_T01_01_DLP01_allele_count.csv.gz \
--ncores 5 \
--rephase_report confidence_c2.csv \
--qcplot qcplot_c2.png \
--heatmap heatmap_c2.png \
--heatmapraw heatmapraw_c2.png \
--csvfile hscn_c2.csv.gz \
--qccsvfile qc_c2.csv \
--Rdatafile signals_c2.Rdata \
--mincells 2 \
--sex FEMALE
