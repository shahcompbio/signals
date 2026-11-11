# Joint cross-cell haplotype phasing for signals — design brief

Self-contained handoff. Open a fresh session with this file to pick the topic up.

## Context

`signals` (`R/callHSCN.R`) infers haplotype-specific copy number per cell. Two passes
of a per-cell, per-chromosome HMM, with a phasing step between them.

Phasing is a **hard label per haplotype block**, shared across all cells:
`phase_haplotypes_bychr` (`callHSCN.R:544`) pools raw allele counts over a
per-chromosome subset of cells and calls the majority allele "A". The subset comes
from UMAP + hdbscan on `state_BAF` (`get_cells_per_chr_local`, `callHSCN.R:438`),
picking the cluster with the highest fraction of allelically imbalanced bins.

Blocks are then summed within 500kb bins (`combinedata.R:312`) and the HMM assigns
B-allele copy number per bin given `BAF = alleleB / totalcounts`.

## Failure mode that motivated this

The cluster-selection step is blind to the events it is meant to find. It clusters on
`state_BAF`, which is derived from the **first-pass mirrored phasing**. On a chromosome
where that phasing is uninformative, an imbalanced cell has first-pass BAF ≈ 0.5 —
indistinguishable from diploid cells on the exact feature used to select it. Circular.

Observed on a 725-cell DLP+ sample (`test_data/`): chr6 was phased from three cells
ranked 226th, 561st and 711th of 725 by allelic imbalance. The three cells that
actually carried imbalance (two 1-copy LOH, one 2:1) ranked 1, 2, 3 and were not used.
Result: those cells flipped A/B every few bins — 129, 117 and 68 segments on a
chromosome that should have one.

## What has already been fixed (in `run_signals.R`, script-level, no package changes)

1. **Detection.** `cell_imbalance_Z()` — Tarone's Z on *unphased* block counts. An
   imbalanced cell is overdispersed relative to Binomial(n, 0.5) no matter how blocks
   are phased, so this needs no phasing, clustering or HMM. Z scales with sequencing
   depth, so cells are scored within per-chromosome depth deciles (robust z within
   stratum); compare `score`, never raw `Z`.
2. **Targeted re-phasing.** `flag_failed_phasing()` + `rebuild_phasing()` — re-derive
   block phases for flagged chromosomes from the top-scoring cells, with an all-cell
   fallback for blocks those cells do not cover. (Without the fallback,
   `phase_haplotypes_bychr` omits uncovered blocks and `format_haplotypes`'s right
   join at `combinedata.R:127` then drops them from *every* cell.) Supplied to
   `callHaplotypeSpecificCN(phased_haplotypes = ...)`, which bypasses the clustering.
   chr6 result: 129 → 1, 117 → 2, 68 → 3 segments; cohort mean 3.14 → 2.65.
3. **Relabelling invariant.** `relabel_min_breakpoints()` + `relabel_unsupported_phase()`
   — where no *other* cell on a chromosome carries detectable imbalance, the A/B label
   is arbitrary, so choose the labelling implying the fewest allele-specific
   breakpoints. Emission-likelihood-neutral (verified: identical to machine precision
   across 725 cells), because signals' emissions are symmetric — `p(B) + p(state−B) = 1`
   exactly, so `dbetabinom(b, n, p) == dbetabinom(n−b, n, 1−p)`.

## The open problem

All of the above still treats phasing as: pick a subset of cells, take a hard majority
vote, done. Two structural weaknesses remain.

**(a) The vote is unweighted and subset-restricted.** Pooling near-balanced cells
destroys the phasing (empirically: top-3 cells excellent, top-50 nearly as bad as the
original). So the method needs a hard threshold on how many cells participate, and the
best threshold differs per cell — on chr6, top-3 was best for the LOH cells while
top-10 was best for the 2:1 cell.

**(b) The transition cost is not allele-aware.** `HaplotypeHMM` (`callHSCN.R:36`) uses a
flat matrix over B with a single `selftransitionprob`, so `2|1 → 1|0` (both alleles
change) and `2|1 → 0|1` (only A changes) score identically. The parsimonious criterion
is the count of allele-specific breakpoints, `1[ΔA≠0] + 1[ΔB≠0]`, which needs a
position-dependent transition matrix — the C++ Viterbi (`src/code.cpp:11`) takes one
fixed matrix.

## Proposed objective

Promote the block phases to free variables and optimise jointly with the per-cell states.

Variables: `φ_h ∈ {0,1}` per haplotype block (shared across cells); per-cell per-bin
allele-specific states `(A_ci, B_ci)` with `A + B` fixed to the HMMcopy total.

Maximise

    Σ_c Σ_i  log P(counts_ci | A_ci, B_ci, φ)  −  λ Σ_c Σ_i  C((A,B)_{c,i−1}, (A,B)_{c,i})

with `C = 1[ΔA≠0] + 1[ΔB≠0]`.

Block-coordinate ascent:

- **Given φ**, the problem decouples over cells and is the existing HMM (with the
  allele-aware transition cost). Already implemented apart from the cost.
- **Given the per-cell states**, the optimal update for each block has a closed form —
  a log-odds weighted vote over *all* cells:

      LLR_h = Σ_c (a0_ch − a1_ch) · log( q_cb / (1 − q_cb) )

  where `q_cb = A/(A+B)` for cell `c` in the bin containing block `h`, clamped by
  `loherror`, and `a0/a1` are that cell's raw allele counts in block `h`. Set
  `φ_h = 1[LLR_h > 0]`.

This is the generalisation of the current rule: the existing method is the special case
where selected cells get weight 1 and everyone else weight 0. Balanced cells have
`q = 0.5` and contribute exactly zero, so no clustering, no `minPts`, no threshold —
the weighting is self-selecting, and `|LLR_h|` is a free per-block confidence.

## Known pitfalls (measured, not hypothetical)

- **Initialisation matters more than anything else.** Seeded from the global mirrored
  phasing, the EM sits in a degenerate fixed point and never recovers (chr6: stayed at
  133/123/68 segments, BAF ≈ 0.46). Seeded from the imbalance-score-selected cells it
  converges correctly and then including all 725 cells is as good as using only the
  top 3. Seed from the score-selected phasing.
- **Do not replace the clustering wholesale.** On chr1, where the clustering works
  (it picked cells ranked 1, 5–12, 17), nothing tried beat it — signals as-run gave
  mean 3.38 segments/cell vs 3.40–3.71 for every alternative. Treat this as a fallback
  triggered by detection, or verify the joint version dominates before switching.
- **Double-dipping.** A cell in the vote influences its own call. Diluted across a large
  cohort, but with one imbalanced cell it self-confirms exactly. Consider leave-one-out
  for the called cell, and report the labelling-invariant imbalance score as independent
  evidence.
- **Sub-chromosomal events.** chr2 improves in aggregate but individual cells barely
  move — the imbalance there looks arm-level, and one whole-chromosome phase cannot
  capture events that differ across arms. Needs `phasebyarm` or segment-level phasing.
- **Identifiability floor.** With one imbalanced cell and a constant-orientation
  constraint, there are H blocks, H free binary phases and zero residual degrees of
  freedom — every cell yields a perfect consistent fit, so nothing is testable. Phasing
  carries information only once two or more cells constrain the same blocks.
- **Unseeded clustering.** `umap_clustering` is called without a seed in
  `get_cells_per_chr_local`, so which chromosomes fail varies between runs (chr6 ranked
  226th in one run, 2nd in the next on identical input). Seed it before benchmarking
  anything.

## Validation setup

- Data: `test_data/` — 725 cells, DLP+, `hscn.csv.gz`, `haplotype_calling_count.csv.gz`,
  `signals.Rdata` (readRDS, not load — it is an RDS despite the extension;
  `readRDS(...)$hscn` has `data`, `phasing`, `haplotype_phasing`, `haplotype_stats`).
- Ground-truth-ish cases: chr6 cells `R58-C15` and `R47-C24` (1-copy LOH, same allele),
  `R59-C30` (2:1). chr1 as a negative control where the existing method works.
- Useful metric: A/B segments per cell per chromosome (should be ~1 per CN segment),
  plus mean BAF against the state's expectation — but note the latter is
  labelling-dependent and is not evidence about the CN call.
- `hap_label` is **not** unique within a chromosome (268 of 6477 blocks on chr6 span two
  bins). Key on `(chr, start, end, hap_label)`, as signals does.
- Beware: signals' fitted `rho` (≈0.028) is a *bin-level* overdispersion. A parametric
  block-level null built from it is too tight — diploid control cells fail it. Use a
  coverage-matched empirical null across cells instead.
