# signals benchmarking: plan

A way to test and benchmark each past and future change to signals — so that a
change can be ranked against the alternatives, not just checked for crashes.

---

## Why this is needed now

Two findings from taking stock of what exists.

**1. The existing simulator cannot test phasing.** `simulate_data_cohort()` emits
`alleleA` / `alleleB` directly — counts already in the true A/B frame — with
`hap_label := 1:.N`, one haplotype block per bin. There is no hidden block phase
to recover, so every simulated dataset enters the pipeline already phased.
`test-hscn.R` asserts >99% `state_AS` accuracy on that data, which says nothing
about phasing. It also simulates whole clones of 10-20 cells, never the regime
where one or two cells carry an event — which is where every remaining problem
lives.

**2. The real-data metric has saturated.** The phase-free pairwise test (do two
imbalanced cells favour the same parental haplotype?) scores runs 14, 19 and 20
all at 1748/1748 on S1. It can no longer tell good methods apart. It is also
structurally blind to the hard case: it needs *two* imbalanced cells, so S1
chromosome 8 (2 of 726 cells) offers at most one informative pair.

So real data can now catch *regressions* but cannot *rank* improvements.
Ranking needs ground truth, and the simulator does not provide it where it
matters.

---

## Three tiers

| tier | answers | status |
|---|---|---|
| **1. Unit / invariant tests** | is the mechanism correct? | exists — ~930 testthat assertions, brute-force Viterbi, `A + B == state` |
| **2. Semi-synthetic benchmark** | which method is *better*? | **missing** |
| **3. Real-data panel** (S1, S2) | did anything regress on real data? | exists ad hoc; formalise |

### Tier 1 — unit and invariant tests

Keep as is. Runs via `dev/run-tests.sh` in the signals image. These check that a
mechanism does what it claims (Viterbi returns the MAP path; the allele-aware
matrix reduces to the flat one at constant copy number), not that it helps.

### Tier 2 — semi-synthetic benchmark

Use real datasets as the template and draw allele counts from a known truth:

- **real geometry** — actual haplotype block positions, blocks spanning two bins
  (the B7 case), real per-cell per-block depth including blocks a cell does not
  cover (the B3 case)
- **known truth** — per-cell `A`/`B` and the true parental haplotype of every
  block, emitted as `allele0` / `allele1` behind a **random hidden per-block
  phase**, so the pipeline has to recover it
- **real noise** — beta-binomial at each sample's fitted overdispersion (S1 0.028,
  S2 0.077)

Taking depth and block structure from real data sidesteps arguments about
whether a simulator is realistic.

Inject events on a grid aimed at failure modes seen in practice:

| axis | levels | motivated by |
|---|---|---|
| cells carrying the event | 1, 2, 5, 20, 200 | chr8: 2 of 726 cells |
| extent | arm, whole chromosome | chr4: p-arm LOH, whole-chromosome gain |
| parent across cells | same, opposite | relabelling collapsed opposite-parent pairs |
| event type | LOH, `2|1`, `3|1` | |
| depth | real, 0.5x, 0.25x | |

Metrics that need truth, so are unavailable on real data:

- parent-of-origin accuracy per (cell, unit)
- **within-segment switch rate** — the `1|2 <-> 2|1` noise this work started from
- one- vs two-allele moves at copy number breakpoints
- **calibration** — does `--rephase_report` flag low confidence where accuracy is
  actually low? At present nothing shows whether the confidence numbers mean
  anything.

### Tier 3 — real-data panel

S1 and S2, fixed seed, a named set of method configurations, one scorer.

| metric | role |
|---|---|
| `A + B == state` in every bin | **hard gate** |
| `state` identical to the HMMcopy input | **hard gate** |
| pairwise parent-of-origin correctness, split same / opposite | report; saturated |
| `sgn` purity — within a (cell, unit), does the cell vote one orientation? | diagnostic; gameable by relabelling (see below) |
| segments, isolated bins, breakpoint allele moves | diagnostic only |
| runtime | report |
| diff against a stored baseline | what moved, and where |

Segment count is diagnostic only: it is a function of the A/B labelling and on
S1 the most correct configuration had the *most* segments.

**Purity turned out to have the same flaw.** It was expected to be the
truth-free metric that still discriminates once the pairwise test saturates.
Validating the scorer showed otherwise: relabelling mirrors everything to
`A >= B`, making orientation uniform by construction, so on S2 the relabelled
run scores 2.5% of cell-arms below 0.9 purity against ~13.5% for every other
run — while being the least correct of them. Like segments, it is only
comparable between configurations that do not relabel.

This sharpens the case for tier 2. On real data, every metric that still
discriminates is a function of the labelling, and the one that is not has
saturated. Truth has to come from simulation.

---

## Steps

### Step 1 — commit the apparatus and the real-data panel *(done)*

- `bench/R/common.R` — readers, arm assignment from signals' own cytoband map
- `bench/R/verdict.R` — phase-free pairwise verdict, cached per dataset
- `bench/R/score.R` — every tier-3 metric into one tidy table
- `bench/R/compare.R` — diff a run against a baseline, per unit
- `bench/run.sh` — run one (configuration, dataset) into a results store keyed by
  git SHA
- `bench/configs.tsv`, `bench/datasets.tsv` — named configurations and datasets
- **validation:** the consolidated scorer must reproduce every number already
  reported for S1 runs 13, 14, 17, 18, 19, 20 and S2 runs 3, 4

### Step 1b — fix pair selection in the verdict *(done)*

Validating step 1 showed the set of pairs the verdict tests was an accident of
input. Rebuilt from the raw HMMcopy reads rather than a signals output, the S1
panel lost its determined pairs on all five arms where relabelling was known to
fail.

The cause was mainly QC-failed cells, not the file-order cut alone: they look
non-diploid and strongly imbalanced on every arm, flooding those arms with
114-120 candidates instead of 10-11, so evidence ranking could not have screened
them out. Fixed by restricting candidates to QC-passing cells and sampling
across bands of `q` on arms over the cap (`strat`). Chosen on the one defect
known on real data: `strat` catches all 8 known relabelling failures on S1, as
every rule does, and 48 on S2 against 2-3 for the alternatives.

Re-scoring history reverses no conclusion, but S2's relabelled run turns out to
get 48 of 52 opposite-parent pairs wrong where the old panel showed 2 of 5. The
old panel had too few opposite-parent pairs to see it; `score.R` now warns when a
panel has fewer than 20.

### Step 2 — semi-synthetic generator

`bench/R/simulate.R`: template a real dataset's geometry and depth, inject a
known event into chosen cells, hide the block phase, emit files `run_signals.R`
reads unchanged. Start with the cells-carrying-the-event axis (1, 2, 5, 20) on
one arm — that alone shows whether C2, `allele_aware`, B12 or anything later
helps in the regime that still fails.

Truth-based scorer alongside.

### Step 3 — grid and report

Methods x regimes x replicate seeds, as LSF array jobs on iris. One report
comparing methods per regime.

### Step 4 — revisit existing tests

`test-complexphasing.R` asserts on mean segment count, which is the metric this
work concluded rewards the wrong behaviour. Reconsider once tier 2 exists.

---

## Design decisions

- **Results are keyed by git SHA** of the signals tree and the script, so a
  number can always be traced to the code that produced it.
- **The pairwise verdict is cached** per dataset. It uses no phasing information,
  so it does not change between methods — and on S1 it takes minutes.
- **Arms come from signals' `cytoband_map`.** An earlier version of the scripts
  hand-coded centromere midpoints; one of those (19.5Mb for chr13, whose
  centromere is ~18Mb) produced a wrong conclusion. The scorer must use the same
  boundaries as the method it is scoring.
- **Outputs are never committed.** Each `hscn.csv.gz` is ~100MB.
- **Scripts depend on `data.table` only**, so they run natively rather than under
  x86 emulation.

---

## Open questions

- **Run on iris?** Each real-data run takes about an hour locally, mostly x86
  emulation on Apple silicon. A grid of methods x regimes x seeds is impractical
  locally and natural as LSF array jobs. If yes, `run.sh` should target
  singularity from the start.
- **Branch.** Infrastructure has so far lived on `test-infrastructure`, separate
  from `phasingimprovements`. But benchmarking evaluates the phasing work, so it
  is most useful when available on that branch.
- **Where real data lives.** Currently copied into `test_data/` and
  `test_data2/`. On iris it could be read in place from the Isabl data lake.

---

## What exists, for reference

| item | where | state |
|---|---|---|
| phase-free pairwise verdict | old job dir | survived; cached for S1, S2 |
| scorer, fails, segments | old job dir | survived |
| breakpoint allele moves, depth, `sgn` purity | `/tmp` | **lost** — reconstruct |
| per-unit confidence report | `run_signals.R --rephase_report` | in the script |
| unit tests | `tests/testthat/` | committed |
| test runner, dev image | `dev/` on `test-infrastructure` | partly committed |
| 20 S1 runs, 4 S2 runs | `test_data/`, `test_data2/` | on disk, gitignored |
