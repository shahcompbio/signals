# signals benchmark

Tools for scoring changes to signals on real data, so that a change can be
compared with the alternatives rather than just checked for crashes. See
[`PLAN.md`](PLAN.md) for the wider plan, including the semi-synthetic tier this
real-data panel cannot replace.

Everything here depends on `data.table` only and runs natively. Only `run.sh`
needs the signals container image.

## Quick start

```bash
# one-off: say where the data is
cp bench/datasets.example.tsv bench/datasets.tsv      # then edit

# one-off per dataset: the phase-free verdict (minutes; cached). Candidates are
# drawn from the copy number source's cells, so it must be QC-filtered: a signals
# output that passes both gates, or the HMMcopy reads plus --cells
Rscript bench/R/verdict.R <allelecounts> test_data/hscn_19.csv.gz bench/cache/S1.verdict.rds

# run a configuration on a dataset (needs an image: dev/build-dev-image.sh)
bench/run.sh c2 S1 --dry-run        # show the command and provenance
bench/run.sh c2 S1                  # ~1h locally under x86 emulation

# score any outputs, including runs made before this harness existed
Rscript bench/R/score.R --verdict bench/cache/S1.verdict.rds \
    --input test_data/<sample>_hmmcopy_reads.csv.gz \
    r19=test_data/hscn_19.csv.gz r20=test_data/hscn_20.csv.gz

# what moved between two runs, and where
Rscript bench/R/compare.R test_data/hscn_19.csv.gz test_data/hscn_20.csv.gz
```

## Files

| file | purpose |
|---|---|
| `R/verdict.R` | phase-free pairwise verdict: do two imbalanced cells share a parent? |
| `R/score.R` | every real-data metric for one or more runs, one row per run |
| `R/compare.R` | per-arm and per-cell diff of a run against a baseline |
| `R/common.R` | readers; arm assignment matching `signals::coord_to_arm()` |
| `run.sh` | run one (configuration, dataset), keyed by signals and script sha |
| `configs.tsv` | named method configurations |
| `datasets.example.tsv` | template for the local, gitignored `datasets.tsv` |
| `history.tsv` | runs made before this harness, with their provenance |
| `data/arm_boundaries_hg19.csv` | p-arm ends from signals' `cytoband_map` |

## Reading the scores

Metrics have different roles, and treating a diagnostic as a correctness metric
has produced wrong conclusions before.

| role | metric | notes |
|---|---|---|
| **hard gate** | `gate_AB_ne_state` | bins where `A + B != state`; must be 0 |
| **hard gate** | `gate_state_changed` | bins whose `state` differs from the HMMcopy input; must be 0. Needs `--input` |
| correctness | `pairs_*`, `pairs_opp_*` | agreement with the phase-free verdict. **Saturated** - the recommended configurations all score 100% |
| diagnostic | `purity_*` | share of a cell's imbalanced bins voting its majority orientation |
| diagnostic | `segments`, `isolated_bins` | |
| diagnostic | `bkpt_two_allele*`, `flat_changes` | allele moves at copy number breakpoints, and changes within constant copy number |

**Every diagnostic is a function of the A/B labelling**, so a method that
relabels can improve it without being more right. Two measured examples:

- on S1 the most correct configuration has the most segments
- on S2 the relabelled run has 2.5% of cell-arms below 0.9 purity against ~13.5%
  for every other run, while being the least correct

Compare diagnostics only between configurations that do not relabel. The
opposite-parent split matters more than the overall figure: relabelling fails
specifically by collapsing opposite-parent pairs, which the overall percentage
hides (99.5% overall, 90% on opposite-parent pairs).

**Where a diff lands matters too.** `compare.R` on runs 19 and 20 (adding
`--allele_aware`) shows the largest change on 14p - 5.5% of bins, against at
most 1.7% elsewhere. 14p is an acrocentric arm with almost no haplotype data, so
its calls there come almost entirely from the transition prior. A change that
moves calls mainly where there is no data to be right about is not evidence of
being more right.

## How the verdict chooses pairs

The verdict can only test pairs among the cells it selects, so selection
decides what the panel can detect. Two rules matter.

**Candidates must be QC-passing cells.** QC-failed cells look non-diploid on
every arm, so they pass the candidate filter everywhere, and they also look
strongly imbalanced (on S1: median `q` 0.21 at a tenth of the depth of real
candidates, against 0.36). Built from the raw HMMcopy reads, the S1 arms where
relabelling fails had 114-120 candidates instead of 10-11; the original
first-30-in-file-order cut kept almost none of the real ones and no determined
pair survived on any of those arms. Ranking by evidence does not screen them out
either. `verdict.R` warns when a source looks unfiltered: 0% of candidates lack
allele counts for a QC-passing source, 27% for the raw reads.

**On arms with more than 30 candidates, sample across imbalance states.** The
default, `strat`, takes cells round-robin across bands of `q`, which picks up
different subclones - and that is where opposite-parent pairs, the ones
relabelling collapses, come from. Judged by the one defect known on real data
(runs that differ only by relabelling), on how many known failures each rule
catches:

| rule | S1 | S2 |
|---|---|---|
| first 30 in file order (original) | 8 | 2 |
| `hash` (stable pseudo-random) | 8 | 3 |
| `lr` (strongest imbalance first) | 8 | 3 |
| **`strat`** (default) | **8** | **48** |

The count of opposite-parent pairs is not the criterion: `lr` had the most on
S1 but the weakest separation, since many of its extra pairs are ones
relabelling does not break.

**Effect on history.** Re-scoring all ten historical runs against the new panel
reverses no conclusion - the order of configurations is unchanged on both
samples - but it shows relabelling did far more damage on S2 than the old panel
could see:

| | old panel | new panel |
|---|---|---|
| S1 r17/r18 (relabelled), opposite-parent pairs wrong | 8 of 78 | 8 of 103 |
| S1 r14/r19/r20 | 0 | 0 |
| **S2 r1 (relabelled), opposite-parent pairs wrong** | **2 of 5** | **48 of 52** |
| S2 c2/r3/r4 | 0 | 0 |

The old caches are kept as `cache/S*.verdict.v0.rds`.

**Permutation noise.** 40 permutations estimate the null's spread noisily:
reseeding flips 1.4% of verdicts on S1, all with `|z| < 6`. Each pair is now
seeded from its own identity, so a panel is exactly reproducible and every
method is scored against the same pairs - but which borderline pairs are in the
panel is partly chance.

## Provenance

Never edit a script while it is running. R reads a script lazily, one expression
at a time from a byte offset, so an edit mid-run makes it resume at the wrong
place - this happened twice while building this panel. `run.sh` therefore runs a
snapshot of `run_signals.R` taken at launch and keeps it beside the results.


`run.sh` keys each output by the signals commit (`-dirty` if `R/`, `src/`,
`DESCRIPTION` or `NAMESPACE` have uncommitted changes) **and** a hash of the
script, and records the image id, flags and wall time in `meta.tsv`. Both are
needed: the same flags have meant different things at different script
versions. S1 run 14 passed `--relabel_unsupported_phase`, but at that version
relabelling only ran under `--rephase_by_imbalance`, so it was inert - visible
in the data as zero bins with `A + B > state`, which is `relabel_core`'s
signature.

## Data

The repository is public, so nothing identifying is committed:

- `datasets.tsv` (sample names) - gitignored; see `datasets.example.tsv`
- `cache/` (verdicts hold per-cell statistics) - gitignored
- `results/` (each run is ~250MB) - gitignored

## Validation

The scorer and verdict were checked against every number previously reported for
S1 runs 13, 14, 17, 18, 19, 20 and S2 runs r1, c2, r3, r4: pairwise and
opposite-parent agreement, `A + B != state`, segments, isolated bins and
breakpoint moves all reproduce exactly. `verdict.R --select file` reproduces both
original verdicts exactly (S1 4,039 and S2 2,626 pairs, maximum `|z|` difference
0). The `strat` panels are exactly reproducible across runs, and with `--cells`
are identical whether built from a signals output or the raw HMMcopy reads. Arm
assignment agrees with `signals::coord_to_arm()` on all 4,375 DLP bins.
