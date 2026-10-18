# Computational social science panel: measurement reliability after the audit

Reviewed 2026-10-18. This is a bounded root-cause review and diagnostic replay on
Julia 1.13.2, two Julia threads. No estimator, package test or original audit
output was edited. Human-perception validation is outside scope.

## Assessment

The independent numerical counterexamples change the status of the earlier
validation evidence. Smooth-data recovery and backend agreement demonstrated
behavior on those datasets; they did not establish reliable measurement under
equivalent encodings. The current audit and specialist oracles establish several
actual defects, including arbitrary score distinctions among mathematically tied
entities and changed neighborhoods under distance-preserving operations.

These defects must be separated from legitimate changes in an endogenous
measure when its reference collection, feature schema, neighborhood size or
statistical target changes. Replacing LOF with a graph or a copula estimator
cannot fix an invariance violation in the estimator being replaced. Conversely,
poor recovery of a planted cluster is not proof of an arithmetic defect.

The 630 flags are dependent stress-test observations, not a failure rate or a
method ranking. The parent's 525/105 split at raw-relative tolerance `1e-10` is
useful triage, not a count of exact ties: genuine close scores, changes in
represented weight ratios, numerical conditioning and different root causes
require separate diagnosis. The existing suite's 911/911 pass does not refute
an independent mathematical counterexample.

## What the metamorphic audit establishes

The harness in [`run.jl`](../../../validation/consistency/run.jl) generally
chooses appropriate comparisons. It aligns entity rows after permutation and
matches feature pairs by names. It carries raw kNN distance scaling and KDE
active-dimension log-Jacobian adjustments separately from percentile comparison.
It records supported model rejections separately. Neighborhood-size changes are
descriptive sensitivity, not failures. That separation should be preserved.

| Comparison | Interpretation and necessary qualifications |
|---|---|
| Reorder rows, exchange feature coordinates, change descriptive IDs/labels | Same observations and intended geometry. Alignment is essential. Numerical reductions, neighbor predicates or normalization can violate the intended identity. |
| Strictly increasing recoding before rank transformation | Same rank order and ties; valid when the recoding preserves represented information. The harness's rank-cubed construction does so at its small sample sizes. |
| Reflect coordinates | Preserves Euclidean distances and the relevant symmetric estimator definitions. Rounding transformed ranks can expose the implementation's failure to preserve that mathematical geometry. |
| Append zero-weight query rows | Positive reference measure should be unchanged for the supported weighted estimators, including preprocessing and normalization. The audit mostly appends existing rows; new off-support queries need separate coverage. |
| Multiply all membership weights by a positive constant | Same mathematical weights after normalization. Floating multiplication can slightly change represented ratios, so generic deviations need attribution; exact constant-support identities hold regardless and provide stronger controls. |
| Change `k`, add positive-weight observations, add/remove features, change transform or estimator | Changes the empirical reference or target. Ranking changes are sensitivity unless an additional mathematical identity establishes otherwise. |

Several passing checks are deliberately conditional. A mid-CDF oracle computed
from package raw values verifies the ranking formula **given those values**.
Stationarity verifies mass against the graph actually constructed. Tree/brute
agreement shares the distance kernel and boundary predicate. None independently
certifies the mathematical raw scores, intended topology or the ordering of
mathematical ties. Specialist rational/polynomial oracles and interventions
correctly fill some of those gaps.

The artifact verification checks deterministic input round-trips, identity
alignment, table consistency and fingerprints. This is valuable provenance,
but it is not statistical validation. Successful script exit means the
investigation completed: the audit intentionally records rather than throws on
violations.

## Focused downstream consequences

The new public-API
[`probe.jl`](../../../validation/consistency/panel/computational_social_science/probe.jl)
replays confirmed small fixtures and measures **headline atypicality** separately
from the original audit's maximum across pair scores and several summaries.
Its source fingerprint matches the original audit. It counts strict reversals
of computed entity order, tie splits and switches at the explicitly illustrative
cutoff `atypicality >= 0.75`. This is an operational illustration, not a
calibrated anomaly rule, estimated false-positive rate or application prevalence.

| Replay | Entities whose headline changes | Maximum change | Strict order reversals | Cutoff switches |
|---|---:|---:|---:|---:|
| Five-row LOF rank reflection | 5/5 | 0.600 | 3 | 2 |
| Eight-row graph feature reversal | 4/8 | 0.125 | 2 | 0 |
| Gaussian product grid plus zero-weight copies | 25/25 | 0.460 | 0 | 6 |
| Exact depth-triangle oracle versus package | 2/3 | 0.333 | 0 | 1 |
| Exact constant-KDE oracle versus package | 3/3 | 0.333 | 0 | 1 |

LOF selects one entity above the illustrative cutoff in each representation,
but it is a different entity: an unchanged selected **count** can conceal a
changed selected **identity**. In the graph fixture, two order reversals do not
change the selected top group at this cutoff, so not every numerical failure
changes every possible downstream decision. On the independent Gaussian grid,
an exactly flat reference should give no distinctions; adding excluded rows
creates six selected entities where there were none. Zero strict reversals there
is not reassuring because the original scores were all tied: 289 computed pair
ties split. The exact triangle and KDE controls similarly create spurious
distinctions, anchored by the specialist algebraic oracles rather than inference
from a small floating discrepancy.

As a deliberate sensitivity control, the probe appends the same Gaussian-grid
copies with **positive** weights. Headline values then change by up to 0.481 and
seven entities cross the illustrative cutoff. This is a changed reference,
not an invariance requirement; similar-sized changes alone do not identify a
bug. The control shows why the reference contract must precede interpretation.

Recounting the preserved audit shows that **603 of 630** flags change at least
one recorded summary beyond `1e-10`; 27 do not. Of those 27, 25 are raw-only
violations and two change pair scores without exceeding that tolerance in the
summaries. For pairwise results the summary maximum includes conventionality,
novelty and dispersion as well as atypicality; it must not be called a headline
ranking change. Existing counterexample rows maximize pair-score discrepancy,
not necessarily headline or raw discrepancy. Full downstream claims require
aligned per-entity comparison, as in this focused replay.

Saved evidence: [effects](../../../validation/consistency/panel/computational_social_science/effects.csv),
[entity scores](../../../validation/consistency/panel/computational_social_science/entity_scores.csv),
[audit counts](../../../validation/consistency/panel/computational_social_science/audit_effect_counts.csv),
[environment](../../../validation/consistency/panel/computational_social_science/environment.txt).

## Blind spots in the earlier recovery experiments and the new audit

1. **Recovery was not an independent correctness oracle.** The earlier six-method
   and graph comparisons used three seeds at one small sample size and defaults.
   They fairly retained poor cases and reported limited scope, but their random
   continuous clouds rarely exercised exact geometric ties, weighted constant
   support, excluded-row preprocessing or large common offsets. No errors were
   observed in the 18 earlier row-permutation checks; that remains a local fact,
   not a general reproducibility result.
2. **The constructed labels favor particular targets.** Reversing two of four
   correlated coordinates affects four of six pairs. Adding four noise features
   expands the summary to 28 pairs while retaining those four intended
   dependence violations. The observed dilution is a schema/aggregation effect,
   not evidence that numerical optimization failed. Bridge observations form a
   small coherent cluster, which local contrast and self-tuned connectivity may
   legitimately rank highly. Neither label design supplies a universal anomaly
   definition.
3. **Higher-order evidence needs its existing qualification.** Identical
   bivariate margins do not identify the full joint distribution. That does not
   prove that every finite-sample aggregation of pair scores is powerless. The
   graph's better recovery in the constructed modulo-constraint case is useful
   descriptive evidence, not a general theorem or validation of its numerical
   neighbor predicates.
4. **The new audit is a stress collection, not factorial Monte Carlo.** Seed and
   nominal size vary together; structured families override size, and several
   families share rank geometries. Method restrictions produce different
   denominators. Flag proportions therefore cannot estimate prevalence, isolate
   sample-size effects or support a league table of methods.
5. **The invariance set is still incomplete.** Only a selected row permutation,
   reversed feature order and simultaneous reflection are tested per case. Add
   individual-coordinate reflections and additional permutations. Raw-space
   graph orthogonal transformations and well-conditioned full-rank depth affine
   transformations are useful where their contracts apply; do not demand such
   invariance from axis-bandwidth KDE or a rank-transform rotation. The new
   singular-support shear counterexample is a particularly important omission.
6. **Weight and preprocessing semantics need targeted coverage.** Existing
   zero-weight additions are duplicate queries. Extreme/off-support queries,
   product weights and weak-but-genuine correlations need independent controls.
   Splitting membership weights over duplicate rows is not generally invariant
   for KDE's leave-one-kernel-out rule and effective sample size. Full refitting
   after deleting a row also differs from holding ranks/bandwidths fixed during
   kernel exclusion. Tests must not impose a different estimator accidentally.
7. **End-to-end coefficient and longitudinal stability was not measured.**
   Unstable scores can change outcome covariates and selected comparison groups;
   this review does not claim a measured coefficient bias or causal consequence.
   Numerical error, sampling uncertainty and deliberate cohort change must be
   assessed separately. Ordinary resampling does not repair a deterministic
   numerical contract violation.

## Priorities and acceptance protocol

**High priority: correct the numerical measurement before interpreting rank
stability.** Preserve proven equalities, coordinate symmetry, supported reference
exclusion and exact neighbor membership. Resolve the singular-support policy
and stable-centering defects identified by the mathematics/HPC panels. The
CDF-of-computed-values test remains useful but is insufficient. A blanket score
tolerance can erase genuine distinctions, as the perturbed-quadrilateral oracle
demonstrates. Do not relabel all 525 small-raw flags as ties or relax the audit
until it turns green.

**High priority: establish a two-layer repair gate.** First, each independently
confirmed fixture must satisfy its own mathematical expectation in raw values,
topology/support and normalized output. Include genuine-near-tie, weak-dependence
and just-off-support controls so an overmerging repair cannot pass by flattening
everything. Second, replay the preserved metamorphic audit into a new output
directory, preserving old artifacts and fingerprints. Review remaining flags by
numerical premise and root cause; any changed tolerance/estimator contract must
be explicit. A change in model semantics is not an invisible bug fix.

**Medium priority: evaluate consequences after correctness.** For repaired
invariant comparisons report aligned headline/tail changes, strict reversals,
tie blocks, and selected identities at declared operational cutoffs. Where a
top-k boundary is tied, report admissible selections rather than row-order
tie-breaking. If uncertainty-aware ranks are introduced, document unresolved
orderings separately from mathematical equalities. Preserve both genuine tiny
distinctions and exact ties in acceptance examples.

**Medium priority: rerun recovery as a separate experiment.** Repeat the existing
seeded recovery cases unchanged before expanding to independently crossed
sample sizes, feature counts, anomaly fractions and clustering. Compare raw
targets and their summaries, retain all failures and avoid tuning on the same
generated labels used to report performance. Reassess the original table after
repairs; its recorded values are not automatically invalid, but they describe
the earlier implementation and cannot certify the corrected one.

For feature/cohort studies, first verify invariance to mere representation and
excluded rows; then deliberately vary included features, positive membership
and time windows. Save full provenance and separate these substantive
sensitivities from numerical failures. Dated cohorts still use different
reference distributions: a score change is not an absolute change in an
entity's underlying properties.

## Procedure and reproduction

Read the full simulation report, harness, generators, artifact verifier and
statistics/mathematics/HPC follow-ups. Inspected the original recovery protocol
and artifacts. Replayed four preserved counterexamples plus the exact KDE/depth
expectations through the public API; added positive-membership reference change
as a non-invariance control. Recounted summary effects directly from saved
`checks.csv`. No new external theoretical claim required a literature search;
the numerical premises are derived or verified in the linked specialist notes.

```sh
JULIA_DEPOT_PATH=/tmp/platypus-review-julia-depot:/home/simon/.julia \
  julia --threads=2 --startup-file=no --project=. \
  validation/consistency/panel/computational_social_science/probe.jl
```

The probe records six comparisons and their entity-level values, and exits
successfully when the investigation completes. It is not a test that claims
the unresolved package defects have passed.
