# signals phasing: roadmap

Status of the phasing work, validated on two DLP+ samples:

* **S1** `SHAH_H002216_T01_01_DLP01` — 725 cells, rho ~0.028
* **S2** `SHAH_H003831_T01_01_DLP01` — 369 cells, rho ~0.077 (the noisier one)

Branch: `phasingimprovements`. Package changes are in `R/`; a large amount of
phasing logic currently lives only in `test_data/run_signals.R`, which is
gitignored — see *Where the code lives*.

Legend: *(pkg)* changes the package · *(script)* lives in `run_signals.R` only.

**This document is itself gitignored** (it lives under `test_data/`). It is the
only record of the analysis below; moving it into the repo — `dev/` or `docs/` —
would be worth doing before the branch is merged.

Last updated after S1 run 18 (`test_data/command18.sh`) and S2 run 3
(`test_data2/command3.sh`).

---

## Root cause

`get_cells_per_chr_local` (`R/callHSCN.R`) picks the cells that phase each
chromosome by clustering on `state_BAF`, which derives from the **first-pass
mirrored phasing**. Where that phasing is uninformative an imbalanced cell has
first-pass BAF ~0.5 and is indistinguishable from a diploid cell on the very
feature used to select it. Everything below is either a consequence of that or
something found while fixing it.

---

## How to judge a change

**Do not use segment count.** It rewards aggressive relabelling by construction
and says nothing about whether a call is right. Measured on S1 against the
phase-free pairwise verdict, the two metrics point in opposite directions:

| S1 configuration | run | overall | **opposite-parent** | same-parent |
|---|---|---|---|---|
| arm re-phasing + relabelling | 13 | 94% | 99% | 94% |
| joint optimisation + relabelling | 17 | 99.5% | 90% | 100% |
| joint optimisation + relabelling, post-bugfix | 18 | 99.5% | 90% | 100% |
| **joint optimisation, no relabelling** | 14 | **100%** | **100%** | **100%** |

Measured on 1,751 well-determined pairs across 20 units (|z| >= 5); an earlier
version of this table used 126 pairs and the same ordering held. S2 agrees:
98% with relabelling (60% on opposite-parent pairs), 100% without.

Segment counts point the other way — the most correct configuration has the
most segments — which is why they are not used here. Use the phase-free
pairwise test (`test_data/phase-confidence.R`) as the quality metric: it asks
whether two cells gained the same parental haplotype, marginalising the block
phase, so its answer cannot move when the phasing changes.

Two further aggregation traps, both of which produced wrong conclusions here:

* **Never summarise per chromosome when phasing per arm.** A genuine phase
  switch at the centromere averages to BAF ~0.5 across the chromosome. On S2
  chr19 this looked like failed phasing; per arm both cells were cleanly
  separated (0.31/0.67) and correct.
* **Check arm boundaries against the cytoband map.** A hand-picked centromere
  coordinate (19.5Mb for chr13, whose centromere is ~18Mb) misattributed a
  p/q junction as an internal flip.

---

## Done — package

### D1. Seed the stochastic steps *(pkg, 0.16.0, committed)*
Three unseeded steps: the subsampling in `min_cells` that sets `minPts`, the
UMAP embedding that picks phasing cells, and the subsampling in `fitBB` (which
chooses binomial vs beta-binomial). Unseeded, **22,265 LOH calls were
run-to-run lottery** on S1. `seed` threaded through `callHaplotypeSpecificCN`,
`callAlleleSpecificCN`, `proportion_imbalance`, `get_cells_per_chr_local/global`,
`min_cells`, `fitBB`. Two runs with `--seed 42` are byte-identical. Also forces
`n_sgd_threads = 0` in `uwot::umap`, since its SGD is only reproducible
single-threaded.

### D2. Viterbi backtrace was not returning the MAP path *(pkg, 0.17.0, committed)*
Both `viterbi()` (C++) and `viterbiR()` seeded the final position with the
*predecessor* of the best final state, then took the per-column argmax of the
forward matrix instead of following the backpointers — so they largely ignored
the transition penalty and emitted isolated single-bin state changes. Verified
against exhaustive enumeration of all state paths. Legacy behaviour available as
`viterbiver = "cpp_legacy"`, verified byte-identical to 0.16.0.

S1 effect: isolated bins 655 → 399, segments 16,040 → 12,150, ~8,000 spurious
LOH calls removed, total CN unchanged. In isolation (data, phasing and rho held
fixed, chr2/6/13): 590 cell·chr pairs lost segments, 8 gained, and the legacy
path never scored higher.

### A1. `phasebyarm = TRUE` produced an empty phasing table *(pkg, 0.17.0, in tree)*
`phase_haplotypes_bychr` filters `chrarm == names(chrlist)[i]`, so the cell list
must be keyed by arm. Neither selection path did: `proportion_imbalance`
accepted `phasebyarm` and forwarded it to neither, and `get_cells_per_chr_local`
accepted it and ignored it. `"4p" == "4"` is never true, so arm phasing produced
**zero** phased blocks and `format_haplotypes` then dropped every block.

Also made `plot_clusters_used_for_phasing` arm-aware: strip the arm for
`chrfilt`, pass the arm's span as `chrstart`/`chrend` so a panel labelled
"chr 19p" covers 19p rather than all of chr19, expose `ncol`, and return the grid
shape so the caller can size the device per row. Rounded the Mb tick labels in
`plottinglist` (independent of arm phasing; any non-integral `chrstart` produced
labels like `100.000001`).

Tests: `tests/testthat/test-phasebyarm.R`.

---

### A3. Five phasing defects fixed together *(pkg, 0.17.0, in tree)*
Fixed in one pass as B4, B8, B5, B3, B7. Verified in the image by deparsing from
the namespace — note that grepping the installed lazy-load DB as text reports
every function missing and is not a valid check.

* **B4** — `get_cells_per_chr_local` now drops hdbscan's noise cluster (`"0"`)
  before ranking, as `get_cells_per_chr_global` already did. Kept if *every*
  cell landed in it. (In `dbscan`, noise is cluster `0`; `-1` is the
  scikit-learn convention. `R/clustering.R:182` remaps `0` through `"ZZ"` to the
  string `"0"`.)
* **B8** — `fillmissing` no longer repairs `A + B > state` by blanking `state`
  and refilling it from a neighbour. The alleles are clamped instead, so total
  copy number is exactly the `CNbins` input.
* **B5** — new `min_propA` argument, **default `0` (off)**. The best `propA` per
  unit and an `unphaseable` list are recorded unconditionally as attributes.
  Default chosen deliberately: five regressions in this work came from changing
  behaviour, so the mechanism shipped inert. See B5b — the evidence for turning
  it on arrived with run 18.
* **B3** — blocks with no counts in the selected cells fall back to the all-cell
  majority instead of vanishing. New internal `phase_with_fallback()`, applied
  to all three branches of `phase_haplotypes_bychr` and to `phase_haplotypes`.
* **B7** — counts are pooled per `(chr, hap_label)`, the phase decided once, and
  mapped back onto every bin the block spans. Under
  `global_phasing_for_balanced = TRUE` a block is assigned to the balanced or
  imbalanced side by a majority of its bins rather than split across both cell
  sets.

Measured effect (S1 run 17 -> 18, S2 `hscn_c2` -> `hscn_3`):

| | S1 | S2 |
|---|---|---|
| bins changed | 7 of 4.4M | 1,271 (`A`), 1,043 (`state`) |
| `A + B != state` before | 2,482 | 269 |
| `A + B != state` after | 2,484 | **0** |
| correctness | unchanged (99.5%) | unchanged (100%) |

Two things to read from this:

* **The invariant is now exact where the package owns the output** — S2 goes to
  zero. S1 does not, because `relabel_core` in the script runs *afterwards* and
  re-breaks it (fixed, see D5).
* **The fixes are near-no-ops on a clean sample.** `fillmissing` was already
  backfilling B3's dropped bins from neighbours, so B3 shows up as the *quality*
  of those bins, not as a change in bin count; B4 and B7 only bite when a noise
  cluster wins or a block straddles a bin boundary. The correctness ceiling on
  these samples is set by the relabelling policy, not by these five.

Tests: `tests/testthat/test-phasing-bugfixes.R`, 24 assertions, including an
old-vs-new comparison driven through `format_haplotypes` (the join that caused
B3's cohort-wide loss) and a check that per-bin keying really would have
disagreed for B7.

## Done — script only

### D3. Phase-free imbalance detection *(script)*
Two statistics, both independent of the phasing so neither can be circular:

* `cell_imbalance_Z` — Tarone's Z on unphased block counts, scored within
  per-unit depth deciles. **Superseded**: it is a *relative* score, so where
  imbalance is common it compares imbalanced cells against imbalanced cells and
  reports nothing unusual. S2 15q scored 8.8 while carrying cells with mixture
  LR 1800–2360.
* `phase_free_imbalance` — symmetric two-component mixture ML giving the allelic
  ratio `q` and the evidence `lr` against balance. **Absolute**, so it does not
  have that failure. Rank cells by `lr`, never by `q`: a thin cell reaches
  q = 0.03 on noise (125 blocks, 661 counts).

### A2. Per-arm detection *(script)*
Generic `unit` (`add_unit`) threaded through scoring, selection, flagging and
rebuilding, exposed as `--rephase_unit {chr,arm}`. Deliberately **independent of
`--phasebyarm`**: arm-level *clustering* is worse (an arm has far fewer bins for
UMAP, giving first-pass `best_used_rank` of 99–227 versus 1–2), while arm-level
*detection* is what catches chromosomes like 4.

S1 chr4: 4p gap 0.0, not flagged, its 17 p-arm LOH cells already the best
available; 4q gap 12.3, flagged, re-phased from the one cell carrying q-arm
imbalance (R59-C30, 12 runs → 1).

### C2. Joint cross-cell optimisation *(script)* — the current best method
`LLR_h = Σ_c (a0 − a1) · w_cb`, a vote over every cell with no threshold and no
cell selection. A balanced cell contributes exactly zero; a cell carrying the
opposite parent has `w < 0`, so its counts flip sign and it *reinforces* the vote
rather than cancelling. `|LLR_h|` is the per-block confidence.

Three implementation points that were each necessary:

* **Seed on evidence, not on the point estimate** (`seed_phasing` uses the
  top-`lr` cell per unit).
* **Weight sign from the HMM state, magnitude from the phase-free ratio.** At
  rho = 0.077 the HMM's hard Viterbi call labels 326 of 369 cells imbalanced;
  those spurious weights swamp the real ones.
* **It cannot be seeded from relabelled states.** Relabelling forces `A ≥ B`, so
  every weight is positive and the update is structurally unable to discover
  opposite parents.

Correctness: S1 100% (1,748/1,748), S2 100% (1,630/1,631) when run without
relabelling — best of every configuration tried on both samples, on both the
original 126/66-pair measurement and the later 1,751/1,631-pair one. S2 chr15
and chr19 both come out correct, where the previous configuration had chr15
wrong.

### D4. Minimum-breakpoint relabelling *(script)* — **superseded, do not promote**
Where nothing pins the phase, the A/B label is arbitrary, so pick the labelling
with the fewest allele-specific breakpoints (`1[ΔA≠0] + 1[ΔB≠0]`, exact by DP
over two orientations per CN segment, emission-likelihood-neutral).

The *criterion* is sound. The *policy* is not: relabelling mirrors undetermined
bins to `A ≥ B`, which makes two cells that gained opposite parents look
identical. Measured cost on S1: 100% correctness without it, 96% with it, 93%
when its gate was mis-tuned; on opposite-parent pairs specifically it fell to
57%. It also produced five separate regressions (partial mirroring across a
scattered free set, conflating "no support" with "no data", a summed rather than
per-block evidence gate, and two wiring errors).

**B2 is the principled replacement**: put the same preference inside the
likelihood, where evidence can outweigh it, instead of applying it as a blanket
post-hoc mirror. `relabel_min_breakpoints` itself is worth keeping as a utility.

### D5. `relabel_core` broke total copy number *(script, fixed)*
The count-free-bin alignment copied raw `A`/`B` from the nearest informative
bin, importing that bin's total: a `state 2` bin received `A=4, B=1` from a
state-5 neighbour. 2,484 bins on S1 run 18, which is why B8's invariant did not
hold end-to-end even after the package was fixed — the script runs last.

The intent was only to align *orientation*, which is a swap, not a copy. Now
swaps the bin's own `A`/`B` when its major allele disagrees with the neighbour's,
preserving `A + B == state` by construction:

```
OLD copy-values : A+B = 5 5 5 3 3 3  vs state 5 5 2 2 3 3  -> broken in 2 bins
NEW swap-only   : A+B = 5 5 2 2 3 3  vs state 5 5 2 2 3 3  -> broken in 0 bins
```

Worth noting as a process failure: this was introduced earlier in the same
session that fixed B8, and the existing 56 assertions did not catch it because
none of them asserted the invariant. 4 assertions added. Separately,
`test-run-signals-functions.R` had its `quit()` mid-file, so anything appended
after line 357 silently never ran — the summary is now last.

---

## Open — decisions

### Drop the relabelling policy?
**Recommendation: drop it.** The evidence is now much stronger than when this
question was first written. On S1 run 18, all 8 remaining disagreements are
opposite-parent pairs reported as same-parent — on 1p, 2q, 5q, 9q and 18p, at
|z| up to 22.6 over 1,345 blocks. Run 14, identical except for relabelling,
gets **all 8 right**. The failure is not noise at the margin; it is the policy
doing exactly what it is designed to do (minimise transitions) in the cases
where the truth has a transition.

The honest replacement for its cosmetic benefit is reporting confidence (B6), so
an undetermined region is *labelled* undetermined rather than smoothed into
looking certain, with B2 carrying the preference inside the likelihood.

### Promote C2 into the package?
In the script each iteration re-runs the whole inference, because it cannot
reach the per-bin states — ~4x runtime. Interleaved with the HMM it reuses the
forward pass, ~1.3x. It also belongs there structurally. Argument against:
validated on two samples, so a third would raise confidence before it becomes
hard to change.

### Move the two phase-free statistics into the package?
`phase_free_imbalance` and the pairwise same-parent test are the most
independently useful pieces and the least entangled with anything else. Nothing
in signals can currently answer "is this cell imbalanced here?" or "did these two
cells gain the same parent?" without first trusting the phasing.

---

## Open — package defects

### B1. `fix_assignments()` is probably redundant *(pkg)*
It smooths singleton bins, which were largely a D2 artefact; 399 remain
genome-wide. Its own comment reads "This is horrible and needs to be
re-written!". Test with it disabled against runs already on disk; if it is a
no-op, delete it.

### B2. Allele-aware transition cost *(pkg)*
Flat matrix over B with one `selftransitionprob`, so `2|1 → 1|0` (both alleles
change) and `2|1 → 0|1` (only A changes) score identically. Wanted
`cost = 1[ΔA≠0] + 1[ΔB≠0]`, needing a position-dependent transition matrix —
`src/code.cpp` takes one fixed matrix, so this touches compiled code. Now the
principled replacement for D4's policy rather than a segment-smoothing fix.



### B5b. `min_propA` will not fix chr8 — and chr8 is not unphaseable *(pkg)*
**This entry previously claimed chr8 has no allelic information and that
`min_propA` was the fix. Both were wrong.** Measured on run 19:

* Only **2 of 726** cells are imbalanced on chr8, so the pairwise test could
  only ever evaluate one informative pair. "0 of 45 pairs determined" counted 45
  pairs drawn from 10 cells passing a loose CN filter, most uninformative by
  construction — not 45 independent negatives.
* Both cells carry real evidence: `R29-C28` on 8p fits `q = 0.050` with
  `lr = 79.9` (genuine LOH), `R59-C20` fits `q = 0.387`, `lr = 171.5`.

`min_propA` cannot help, for two independent reasons:

1. **It is structurally bypassed.** Once the script supplies `chr_cell_list` or
   `phased_haplotypes`, `callHaplotypeSpecificCN` takes the
   "Using user provided cell list" branch and never calls
   `proportion_imbalance()`, where `min_propA` lives. With `--joint_phasing` it
   would only affect the first-pass seed, which C2 then overwrites.
2. **Its fallback is the wrong action here.** Below the floor a unit pools
   *every* cell, which on chr8 would dilute 2 informative cells among 726
   balanced ones. Pooling is defensible when nothing is imbalanced; chr8 is the
   opposite case — few but strong.

The real problem is **too little evidence**, not none. Counts from voting cells:

| unit | voting cells | blocks | counts |
|---|---|---|---|
| 13q | 675 | 3,810 | 3,608,059 |
| 8q | 29 | 2,495 | 194,642 |
| 8p | 10 | 1,237 | 41,470 |

8p has ~1% of 13q's evidence. Worse, the cell with the strong imbalance
(`q = 0.05`) contributes only 266 counts over 76 blocks at a median of 3 per
block, while the cell with plenty of counts has weak imbalance. So ~1 nat per
block is available on 8p, i.e. roughly 2.7:1 odds — far below determination.

**No phasing method will make chr8 look clean on this data.** The honest
outcome is B6: report the region as low-confidence rather than emit a confident
arbitrary phasing. `confidence_19.csv` already measures this correctly (8p
`median_abs_llr` 0.00, `frac_determined` 0.066); it just is not surfaced in the
output.

### B12. C2 votes with the current, possibly wrong, per-bin orientation *(script)*
Found while testing the above. `joint_phase_update` sets each cell's vote
direction from `sgn = sign(A - B)` **per bin**, taken from the current state
assignment, so a cell whose phasing already flips votes against itself:

| unit | voting cells | median sgn purity | cells below 0.9 |
|---|---|---|---|
| 13q | 675 | 1.000 | 1.9% |
| 8p | 2 | 0.853 | 50% |
| 8q | 3 | 0.618 | 100% |

`R29-C28` on 8p has 65 bins A-major and 27 B-major, so 29% of its bins vote
backwards. This is a fixed point: flipping states produce flipping votes produce
flipping phases. The imbalanced cells are *not* swamped by the balanced majority
(a balanced cell has `sgn = 0` and contributes exactly nothing; on 8p the two
imbalanced cells carry 79% of all weight) — they cancel themselves.

Forcing one orientation per (cell, unit) by majority changes **42% of 8q blocks
and 13% of 8p**, and raises median per-block `|LLR|` by ~50%, while touching
only 1% of 13q — so it is close to a no-op where the phasing is already
determined and material where it is not.

Worth fixing on its own merits, but note it does **not** rescue chr8: median
`|LLR|` goes 0.8 -> 1.2, still undetermined. Fixing it needs care, since a
per-unit majority is too coarse for a cell with a genuine mid-arm switch — the
orientation should come from the cell's own counts rather than from its current
labels.


### B6. Phasing confidence in the output *(pkg)* — half done
Exists as `hscn$phase_confidence` in the script and as
`test_data/phase-confidence.R`; should be a first-class output field plus a
per-cell phasable/not-phasable annotation in `qc_per_cell`.



### B9. No leave-one-out for phasing cells *(pkg)*
Cells used for phasing influence their own calls. Simulated bias: phasing cells
reported minor BAF 0.189 where truth was 0.25, while held-out cells with the
identical true state reported 0.357.

### B10. `$phasing` has two incompatible shapes *(pkg)*
`plot_clusters_used_for_phasing` expects a cell list per unit;
`plotHeatmapQC` (`heatmap_plot.R:2502`) expects `$cl$clustering` and `$prop`.
Not currently in the pipeline's path, but a trap.

### B11. `plotting.R:253` tests `chrstart` twice *(pkg)*
`(!is.null(chrstart) | !is.null(chrstart))` — the second was presumably meant to
be `chrend`, so passing only `chrend` silently does nothing.

---

## Open — projects

### C1. Segment-level phasing *(pkg)*
Arms are a cheap approximation that catches whole-arm events and misses focal
ones; S1 chr2's imbalance is sub-chromosomal but not specifically arm-level. The
principled version segments on the copy-number breakpoints already present.
Note that C2 partly removes the need: per-cell-per-bin weights are already
local, so the unit only matters for the seed and for reporting.

### C3. Use the reference-panel inter-block phase *(pkg)*
signals treats blocks as independent, discarding SHAPEIT's between-block phase.
Used as a prior it would help precisely where **no** cell is imbalanced — the one
regime in which every cell-based method, C2 included, is helpless. Highest
ceiling on this list and the least explored.

---

## Where the code lives

`test_data/` is gitignored, so all of the following exists only there:

| function | status |
|---|---|
| `phase_free_imbalance`, and the pairwise test in `phase-confidence.R` | move to pkg (highest value) |
| `seed_phasing`, `joint_phase_update` | move to pkg as a loop inside `callHaplotypeSpecificCN` |
| `relabel_min_breakpoints` | move as a utility; keep the criterion, drop the policy |
| `cell_imbalance_Z`, `select_phasing_cells`, `flag_failed_phasing`, `rebuild_phasing`, `orient_cells` | superseded by C2; do not move |
| `relabel_unsupported_phase`, `relabel_core` | costs correctness (D4) and has broken the CN invariant once (D5); do not move |
| `add_unit`, `expand_phasing_to_units` | internals only if the unit concept moves |
| `phasebyarm_supported` | script/package skew probe; meaningless in the package |

60 assertions in `test_data/test-run-signals-functions.R` cover these — note
its `quit()` must stay at the end of the file. `test_data2/` is now in
`.gitignore`.

Evaluation scripts used for the tables above are *not* checked in either: the
phase-free verdict (permutation-calibrated pairwise z) and the scorer that
grades an `hscn` output against it. They are the measurement apparatus for every
correctness claim in this document and should move with the statistics in B6.

---

## Cross-cutting caveats

- Validated on **two** samples. Thresholds (`--rephase_min_score`,
  `--phase_llr_min`, `--phase_min_frac`) were calibrated on S1 and behave
  differently on S2, which has ~3x the overdispersion.
- Tarone's Z is not a p-value — genome-wide overdispersion inflates the baseline
  for every cell. Only the ranking within a depth stratum is meaningful, and
  that ranking fails when imbalance is common.
- `rho` from `fitBB` is a **bin-level** overdispersion. A parametric block-level
  null built from it is too tight and fails its own diploid controls; use a
  coverage-matched empirical null.
- `average_distance` and mean-BAF-versus-expected are functions of the A/B
  labelling, so they are phasing diagnostics only — never evidence about a copy
  number call. Use the labelling-invariant mixture ML for that.
- The per-chromosome clustering was unseeded until D1, so any comparison against
  a pre-0.16.0 run is confounded: which chromosomes fail varied between runs
  (S1 chr6 ranked 226th in one run and 2nd in the next on identical input).

---

## Suggested next step

Drop relabelling and enable the floor, on both samples:

```
--seed 42 --viterbiver cpp --rephase_unit arm --joint_phasing --joint_iters 4 \
  --phase_llr_min 5 --phase_min_frac 0.5 --min_propA 0.05
```

Blocked on adding a `--min_propA` flag to `run_signals.R`. Expected: S1 goes
99.5% -> 100% (the 8 opposite-parent pairs), and chr8 is reported unphaseable
instead of phased.
