# Numerical repair acceptance: computational social science

**Final integration addendum:** the lead completed the frozen-source replay:
14,637/14,637 checks pass across 84 datasets, with zero flags and zero unexpected
errors. The full package suite passes 3,701 assertions and the independent
comparator suite passes 976. Input/check/tolerance comparison and artifact
integrity verification pass. The [integration record](numerical_repair.md)
records closure and the explicit KDE normalization-contract change; the
review-stage requirements below document the acceptance procedure.

Review date: 2026-10-18. Julia 1.13.2, two Julia threads. This acceptance review
concerns numerical measurement reliability and reproducibility; human-perception
validity remains outside scope. No package source, tests, historical artifacts or
other documentation was edited by this reviewer.

## Acceptance assessment and evidence boundary

**The focused acceptance gate passes on final source `a60109d…`.** All five
previously selected defect comparisons have zero discrepancies; all eight
observations in the independent mixed-underflow oracle now have exactly matching
Float64 percentiles, including the previously incorrect reference row 4. The
positive-reference control remains appropriately sensitive. This closes the
blocker found on intermediate source `c7c96e…`; it does not replace the complete
integration suite/audit gate or establish universal transcendental exactness.

The first repair iteration passed 3,616/3,616 package checks and the documentation
build but retained two real-KDE unit-invariance flags in the complete audit.
The revised ordering contract was introduced to address those flags. Its initial
fallback introduced the independently detected mixed-underflow counterexample,
which is now repaired and retained as a control. Final
counts and fingerprints belong to the integration report. Passing a local replay
or a broad finite audit cannot justify overlooking a mathematical contradiction.

## Independent replay of downstream effects

I reran the six public-API comparisons from the previous
[computational social science follow-up](panel_followup/computational_social_science.md)
on the frozen repaired source. The probe was copied to
`/tmp/platypus-css-repair-acceptance/replay.jl` with only its input/output path
constants changed; original investigation artifacts were not overwritten.
The results and reproducible probe are now preserved in
[`validation/consistency/repair_panel/computational_social_science/`](../../validation/consistency/repair_panel/computational_social_science/),
including `effects.csv`, `entity_scores.csv` and `environment.txt`.
The first replay was preserved under
[`repair_iterations/first_full/panel_css/`](../../validation/consistency/repair_iterations/first_full/panel_css/)
before refreshing the current six comparisons. The refreshed source fingerprint is
`a60109d4a70036f7b7474d7a410e105b17f182f869b218332a45e5a0acb05f8c`.
The density specialist's separate hash
`67f545b13c28e0716b6115d357ce92fa10e96aafba5a673801e101a3a9a29bca`
identifies `src/density.jl`, not the concatenated package-source fingerprint.

The operational cutoff remains `atypicality >= 0.75`, purely illustrative and
uncalibrated. These are actual headline changes, not maxima over arbitrary pair
scores or several summary fields.

| Comparison | Previous maximum headline change | Repaired maximum | Repaired order reversals | Repaired cutoff switches |
|---|---:|---:|---:|---:|
| Five-row LOF rank reflection | 0.600 | 0 | 0 | 0 |
| Eight-row graph feature reversal | 0.125 | 0 | 0 | 0 |
| Gaussian product grid plus zero-weight copies | 0.460 | 0 | 0 | 0 |
| Exact triangle-depth oracle versus package | 0.333 | 0 | 0 | 0 |
| Exact constant-KDE oracle versus package | 0.333 | 0 | 0 | 0 |

All five comparisons also have zero changed entities and zero split or newly
created computed ties. In particular, adding zero-mass rows no longer creates
six selected entities from a mathematically flat Gaussian reference, and LOF
reflection no longer swaps the identity of the selected entity.

For the positive-membership control, appending the same two rows with positive
mass changes all 25 original headline scores, with maximum change 0.481481 and
seven cutoff switches. This remains legitimate reference sensitivity. Its
presence alongside restored invariants is evidence that the repair distinguishes
excluded-query arithmetic from a changed empirical population; it does not make
scores invariant to every dataset change.

## What is preserved, and what deliberately changed

The LOF definition, kth-neighbor target, kernel-contribution leave-one-out KDE,
Gaussian and empirical-beta copula definitions, and graph connectivity target
remain recognizably the same statistical objects. The numerical mechanisms now
preserve selected identities and intended geometry rather than defining ties by
rounded distances. The four pair summaries and the graph's separate output
contract retain their meanings. The graph still has no invented pairwise
novelty/conventionality summaries.

There are three changes that must not be described as invisible implementation
details:

- **Depth policy changed.** Exact algebraic covariance rank and strict affine
  support replace spectral truncation and its inconsistent residual tolerance.
  Very thin full-rank references remain full rank. Numeric `depth_rtol` now
  raises an actionable error. This is a disclosed estimator-policy change for
  near-singular cases; earlier analyses using truncation require refitting and
  reconsideration, not an assumption of identical results.
- **Neighbor cost changed.** Requested KD-tree and brute backends currently both
  execute the conservative filtered exact scan. Metadata distinguishes request
  and execution. This supports correctness but removes the earlier acceleration
  claim; candidate work is quadratic even without a dense retained distance
  matrix. Previous timing artifacts do not describe the repaired implementation.
- **Normalization input changed for real-coordinate KDE.** Physical-unit
  Float64 `raw_scores` are retained, but their rounded ties need not be the ties
  used by the empirical CDF. Public `normalization_scores` supplies the ordering
  labels, and `normalization_basis` identifies their policy. These labels are
  ordinal: their distances and magnitudes have no density interpretation and
  cannot be compared across pairs or fits. Export them explicitly with
  `scoretable(result; normalization=true)` when reproducing normalization.

## Revised CDF contract, detected blocker and final resolution

The independent 512-bit analysis of the two residual audit flags shows that
genuinely different densities can have identical physical-unit Float64 log
outputs after a common unit shift. Therefore simultaneously requiring exact
ranking of those rounded physical logs and preservation of the genuine order is
impossible on that fixture. Exposing a separate order representation is a
defensible contract correction, rather than an excuse to merge genuine
near-ties. The audit must explicitly verify the weighted CDF of
`normalization_scores`; it must no longer claim that this is always the CDF of
`raw_scores`. Physical-unit raw-score checks and the independent density-order
oracles remain necessary to connect those labels to the intended estimator.

The first implementation on source `c7c96e…` compared positive represented
dimensionless masses exactly but put **every all-underflow row below every positive-mass row**.
That separation is invalid with nonuniform weights. A tiny weight can reduce a
positive represented exponential's contribution below a much larger Gaussian
mixture whose direct exponential terms all rounded to zero.

Public-API fixture: x coordinates `[0,100,38,139]`, y coordinates all zero,
`transform=:none`, `method=:kde`, `bandwidth=(1,1)`, and reference weights
`[1e-20,1,0,1e-10]`. The third row is an excluded query. A separately evaluated
2048-bit Gaussian-mixture oracle, using the represented inputs and explicit
bandwidth, gives the following comparison against failing source `c7c96e…`:

| Quantity | Heavy reference, row 2 | Query, row 3 |
|---|---:|---:|
| Oracle log dimensionless kernel mean | approximately -760.5000000001 | approximately -768.05170186 |
| Package physical raw log density | approximately -761.41893853 | approximately -768.97064039 |
| Package ordinal normalization label | 2 | 4 |

The raw logs and oracle agree that the query density is lower, by a substantial
log difference of about 7.55. That intermediate normalization order reversed
them. The query's weighted percentile should be approximately `1e-20`; that source returned
**1.0**, reversing its headline atypicality by almost one. This is neither a
physical-log collision nor a subtle unresolved transcendental near-tie.

The wider-weight control `[1e-300,1,0,1e-200]` also failed there. Its query's oracle
log-kernel mean is approximately -1412.776, well below reference row 2 at -760.5,
yet the query again received percentile 1. The saved oracle also reports a
minute ordering difference between reference rows 2 and 4; resolving that
separate tiny distinction is **not** needed to establish the primary query
failure and should not be conflated with it.

The historical failure is documented in the
[transcript-derived record](../../validation/consistency/repair_panel/computational_social_science/mixed_underflow_c7_transcript.txt).
The original full-precision `c7c96e…` CSV/environment files were overwritten by
integration replays before this reviewer could archive them. That record is
explicitly reconstructed from recorded tool outputs, not presented as a recovered
original CSV. The directory named `mixed_underflow_intermediate` preserves the
actual files found before the final replay: they already carry `a60109d…` and
eight passing oracle values; its README states this rather than inventing an
earlier failed state. The separate `first_full/panel_css` archive remains intact.

Final independent evidence is in the reproducible
[`mixed_underflow_probe.jl`](../../validation/consistency/repair_panel/computational_social_science/mixed_underflow_probe.jl),
[observations](../../validation/consistency/repair_panel/computational_social_science/mixed_underflow_final/observations.csv)
and fingerprinted environment file. The oracle does not call the package kernel
or ordering implementation. It asserts only its independently established
mathematical premises and records the package discrepancies.

The final repair uses a common scalar representation: an exponential offset
and an exact rational mass. Ordinary positive-mass sums retain their existing
distinctions. When a kernel underflows, its unweighted exponent is shifted and
the original weights remain outside the exponential/log calculation, preserving
small conditional-mass differences. Different offsets are compared through
adaptive directed local-MPFR bounds on the log mass ratio. Thus, an underflow
flag is no longer an order key. These comparisons certify the ordering of the
represented scalar keys; they do not make all kernel evaluations exact real
exponentials. No caller-global precision state is changed.

I independently reran both scripts after the specialist confirmed the source
freeze. Every one of the eight saved percentiles agrees exactly with the
2048-bit oracle after Float64 export. The wide-weight query returns `1e-300`
instead of 1, and the moderate-weight query returns approximately `1e-20`
instead of 1. Both row-4 percentiles are also correct. All five original defect
comparisons remain at zero discrepancy, and positive-weight appending still
changes 25 headline values with seven illustrative cutoff switches. No additional
confirmed defect emerged from this final bounded replay.

The separate geometry contracts are appropriate. Source-backed unweighted ranks
retain the intended integer midrank lattice for neighbor membership. Raw and
z-score neighbors use the represented transformed binary coordinates. Source-backed
raw/z-score depth uses original coordinates to preserve its affine identity.
Arbitrary decimal multiplication can change represented input distances, so
failure of an imagined real-number similarity is not automatically a package
defect. The revised tests appropriately preserve those decimal examples and
check each against its own independent rational geometry, while using exactly
preserved binary scaling for strict invariance controls.

## Review of numerical and provenance acceptance controls

I inspected the mathematics, statistics, HPC, Julia and communication repair
reports, the integrated transformation/scoring/provenance paths, and the new
`numerical_depth`, `numerical_density`, `numerical_neighbors`,
`numerical_contracts`, `numerical_precision` and `numerical_limits` tests.
Their specialist-run results are evidence reported by those reviewers; my
independent execution in this pass is the six-comparison replay above.

The controls address the important distinction between preserving exact equality
and erasing genuine distinctions: the perturbed quadrilateral, weak nonzero
Gaussian dependence, beta-polynomial perturbation and near-but-unequal neighbor
distances remain separate cases. Constant-support, Cartesian-product and
triangle identities provide expectations independent of the raw values emitted
by the implementation. Exact neighborhood oracles certify topology separately
from subsequent floating affinities and stationary mass.

Reference-mass provenance is materially improved. Original represented weights
reach transformations, moments and percentile masses; exported normalized
weights are not silently substituted for them. Metadata binds both to ordered
entity identifiers, permitting an explicit ID join after table sorting. It also
records effective backend and depth rank/support policy. Detached snapshots and
the outcome result's measurement metadata prevent accidental aliasing from
rewriting the analysis record. Full unique pair coverage is now validated, so
duplicating one pair cannot silently change the aggregation's feature weighting.

Metadata remains a record, not proof: it cannot verify a user-supplied dataset
version, cutoff or missing-data decision. Package version `0.2.0` alone does not
distinguish the pre-repair and repaired working trees. Preserve source revision
or fingerprint, dependency environment, input identity and analysis script
alongside the table. Replaying with rounded exported weights alone is not
necessarily bit-identical to replaying original supplied mass; use the recorded
`source_reference_weights` when exact reproduction matters.

## Remaining limits and final integration gate

The mixed-underflow ordering blocker above is resolved on the final fingerprinted
source. Exact predicates
and rational moments do not make arbitrary transcendental log densities, LOF
ratios or affinities exact. Distinct mathematical depths can round to the same
Float64 output. The documentation correctly retains those limits and avoids
claiming an arbitrary epsilon has solved every score-ordering problem.

Before declaring integration complete:

1. Retain the mixed-underflow oracle and finish the broad audit into the new directory, preserve its
   fingerprint and compare raw, pair and summary discrepancies with the original
   artifacts. A changed or rejected expectation needs an explicit mathematical
   explanation; no silent tolerance relaxation. Record the changed public CDF
   input explicitly rather than describing the final oracle as unchanged.
2. Include the independent mathematical fixtures, genuine-close-score controls,
   provenance tests and limits/precision checks in the complete passing suite.
   Both requested backends now use the same algorithm, so their agreement alone
   adds little independent correctness evidence; retain the rational oracles.
3. Distinguish integration correctness from the earlier recovery experiments.
   Their observed AUROCs, schema sensitivity and timings describe their recorded
   implementation/runtime. Recompute them before presenting those numbers as
   measurements of the repaired source. Correcting numerical invariance does
   not establish a universal anomaly ranking or remove the reported bridge-group
   and irrelevant-feature limitations.
4. Keep scaling claims bounded by the corrected algorithm. Exact arithmetic's
   bit lengths, allocations, ties and the conservative scan deserve workloads
   beyond the small correctness audit before recommending large empirical jobs.

The endogenous-reference design is unchanged. Dated cohorts, feature selection,
positive membership changes and outcome alignment remain analyst decisions.
Restored invariance to equivalent representations makes those substantive
comparisons more interpretable; it does not make independently fitted cohort
percentiles an absolute scale, create historical prediction semantics, or supply
causal identification or uncertainty estimates.

Reproduce the focused replay without altering historical panel observations:

```sh
julia --threads=2 --startup-file=no --project=. \
  validation/consistency/repair_panel/computational_social_science/probe.jl
julia --threads=2 --startup-file=no --project=. \
  validation/consistency/repair_panel/computational_social_science/mixed_underflow_probe.jl \
  /tmp/platypus-css-mixed-underflow-replay
```

The mixed-underflow probe accepts an optional output directory; use a fresh one
to avoid overwriting prior evidence.
