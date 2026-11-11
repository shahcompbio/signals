#!/usr/bin/env bash
# Joint cross-cell optimisation (roadmap C2), third attempt.
#
# Run 14 skipped relabelling entirely - it lived inside the old branch - so the
# comparison against run 13 was never like-for-like, and chromosome 8, whose
# appearance depends wholly on that step, showed its raw flipping.
# Run 15 died on `unit_by` not being defined outside that branch.
#
# Both fixed. Relabelling now runs for either phasing path and, under C2, is
# gated on the actual per-block evidence rather than a per-cell score proxy.
# The gate is the FRACTION of a bin's blocks reaching |logLR| >= 5, not the
# summed |logLR|: a sum grows with block count, so run 15's summed gate freed
# only 70 of 6088 bins while 30 weak blocks summing to 15 passed it.
#
# On run 15's own confidence table this pins 16p, 16q, 13q and 17q (74-98% of
# blocks determined) and frees the other 37 units, where the median block
# carries |logLR| under 4. 8p sits at 0.00 - correctly saying there is no
# phasing information there at all.
#
# Joint cross-cell optimisation of the block phases (roadmap C2).
#
# Replaces --rephase_by_imbalance entirely. No score threshold, no cell
# selection, no orientation forcing: every cell votes for every block weighted
# by its own local imbalance, so a cell carrying the opposite parental haplotype
# gets a negative weight and REINFORCES the vote instead of cancelling.
#
# Why this run exists: on the noisier sample, chr15/19 have two cells with gains
# on opposite parents. Their per-cell imbalance scores (8.8, 7.4) fell just under
# the threshold of 10, so relabel_unsupported_phase declared the units
# unsupported and mirrored every cell to A >= B - erasing which parent was
# gained and making the two cells look identical. The phase-free pairwise test
# says opposite (z = -8.6, -9.2).
#
# Validated on that unit in isolation:
#   current pipeline      0.303 / 0.391  same side   WRONG
#   joint optimisation    0.355 / 0.644  opposite    correct
#   phase-free truth      ~0.31 / ~0.75
#
# Two design points worth knowing:
#  * the seed is the cell with the most imbalance EVIDENCE (mixture LR), not the
#    most extreme ratio - a thin cell hits q = 0.03 on noise and poisons every
#    iteration;
#  * the weight takes its sign from the HMM state but its magnitude from the
#    phase-free ratio, because at rho = 0.077 the HMM calls 326 of 369 cells
#    imbalanced and those spurious weights swamp the real ones.
#
# Each iteration is a full HSCN inference, so this is ~4x the runtime. It stops
# early once no block flips (3 iterations on both samples).
#
# Writes per-unit phase confidence (median |logLR| per block, and the fraction
# of blocks determined) to the --rephase_report path.

docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --joint_phasing \
  --joint_iters 4 \
  --phase_llr_min 5 \
  --phase_min_frac 0.5 \
  --relabel_unsupported_phase \
  --rephase_report confidence_16.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_16.png \
--heatmap heatmap_16.png \
--heatmapraw heatmapraw_16.png \
--csvfile hscn_16.csv.gz \
--qccsvfile qc_16.csv \
--Rdatafile signals_16.Rdata \
--mincells 2 \
--sex FEMALE
