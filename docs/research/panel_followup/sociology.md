# Sociology panel: reference definitions and operational consistency

Reviewed 2026-10-18 on Julia 1.13.2. Scope: the internal scientific contract,
representation invariance, reference membership, and reproducible substantive
workflows. Human-perception resonance is excluded. No estimator, package test,
existing audit output, or other panel report was changed.

## Assessment

The package's declared measurement scope is substantially coherent. It specifies
a supplied collection, feature schema, transformation, method and reference
weights; it does not claim to discover a unique social population. Labels do not
define reference groups. Its documentation correctly distinguishes local
isolation, centrality, dependence density and graph connectivity, and explains
the direction of novelty and atypicality. These strengths should survive repair.

The confirmed numerical failures violate that contract within a fixed numerical
target. Reordering columns or appending rows with no reference mass cannot
legitimately change the collection's substantive definition. Differences caused
by those operations should not be presented as sensitivity to a research design.
In contrast, replacing a reference cohort, adding measurement features, or
switching estimators changes the target and need not preserve the scores.

The evidence supports prioritizing numerical correctness over a new estimator
family. No finding here establishes that another sophisticated estimator would
recover the same target more faithfully. A method substitution must carry a new
target specification and separate evaluation, not serve as a way to make the
current failing properties disappear.

## 1. High: final-score distinctions lack an adequate numerical contract

The statistics, mathematics and HPC panelists provide independent oracles and
interventions, rather than solely cross-method comparisons:

- Equal simplex depths and constant weighted KDE densities become unequal
  percentiles; the statistics panel also confirms equal beta-copula endpoints.
- Column order and reflection can alter exact kth-distance neighborhoods and
  graph edges. This is a raw geometry defect, not just a display/ranking issue.
- A zero-weight augmentation changes Gaussian-product-grid rankings despite
  exact factorization of the positive reference.
- Depth support classification differs after an exact shear, while large
  representable translations degrade otherwise identical centered geometry.

These results matter to the package's substantive workflow because exported
scores and outcome predictors use normalized summaries. An arbitrary ordering
among equal mathematical scores can enter a median, a tail quantile, or an
outcome design matrix. That propagation follows from the source definitions;
this panel has not estimated changes in real empirical conclusions or run a new
outcome-effect study. A high-level caveat that scores are descriptive does not
remove the obligation to compute the stated descriptor consistently.

Acceptance should require raw identities, neighborhood/support decisions, and
final normalized outputs independently. An empirical-CDF test conditional on
already rounded raw scores cannot establish the intended ordering. A universal
epsilon tie rule would introduce another unacknowledged measurement decision:
the mathematics panel's genuine-near-tie control shows why small differences
cannot automatically be discarded. Any approximation policy must identify its
domain, determinism, error rationale and effect on reported ties.

## 2. Separate invariance obligations from changes to the measurement target

| Operation | Appropriate expectation and boundary |
|---|---|
| Reorder entities, retaining ID/weight/feature alignment | Same entity scores; identity alignment is essential. |
| Reorder or rename features, retaining the same columns and corresponding bandwidths | Same named pair scores and summaries; full-feature graph unchanged. |
| Change descriptive labels/IDs without changing alignment | No change in numerical scores. |
| Strictly increasing recoding under rank transformation | Same rank geometry if represented order and ties are preserved. No corresponding general requirement for raw or z-score distances. |
| Reflect a coordinate under the stated symmetric scoring rules | Same aligned scores; this tests the rule's geometry, not whether a variable's substantive meaning changed. |
| Add zero-weight queries to supported weighted estimators | Existing reference fits and existing entity scores unchanged; queries receive no reference mass. |
| Multiply all membership weights by a positive common constant | Same mathematical target. Distinguish harmless representation changes in weight ratios from exact identity fixtures that must still hold. |
| Change units jointly | Apply the method's correct raw-scale transformation: distances and densities are not numerically invariant in their original units. Preserve normalized ordering when the method is equivariant and parameters transform consistently. |
| Add positive-weight observations, alter relative membership, or change cohorts | Scores may change legitimately: marginal transforms, fits and percentiles use a different reference. |
| Add, duplicate, remove or mix measurement features | Scores may change legitimately: the pair collection or full-feature metric changes. General feature mixing is not an invariance of a pairwise summary. |
| Change k, bandwidth, regularization, support policy or estimator | Treat as specification sensitivity. A parameter change is not a numerical implementation cross-check. |
| Split one entity into multiple positive-weight records | Not an automatic invariance: entity-level leave-one-out and effective sample size can change even if empirical mass is preserved. |

Full-rank Mahalanobis affine invariance does not establish unrestricted affine
invariance for a thresholded singular-support approximation. The newly confirmed
axis/oblique support inconsistency needs an explicit support policy, as described
by the mathematics panel, rather than an overbroad blanket invariance promise.

## 3. Medium: record granularity needs an explicit place in the reference contract

The guides correctly distinguish membership weights from survey weights and
replicate counts. They should also state that a row is the unit of kernel
exclusion and effective-sample-size calculation. Repeated records can represent
different entities or multiple observations of one entity; the package cannot
infer which. A researcher must settle that resolution before constructing the
reference, and record any deduplication/aggregation rule.

The small [contract probe](../../../validation/consistency/panel/sociology/probe.jl)
provides a positive control for this distinction. With raw coordinates
`(0,0), (1,1), (2,0)`, weights `[2,1,1]`, and fixed KDE bandwidth 1, splitting
the first row into two identical rows with weights `[1,1]` preserves the
weighted empirical distribution. However, the focal log density changes from
`-3.217762559451067` to `-2.5288833906330743`, and its typicality changes from
0.5 to 0.75. Both log densities agree with independent analytic Gaussian-kernel
calculations. The original row excludes its whole kernel contribution; either
split row leaves its duplicate in the reference. Effective sample size changes
from 8/3 to 4, although the fixed bandwidth isolates the exclusion effect here.

This is **not a new numerical bug**. It prevents an invalid replication-invariance
test and clarifies the estimator's dependence on the declared entity unit.
Acceptance: retain a documented control distinguishing zero-mass augmentation
from positive-mass record splitting, explain KDE kernel-only leave-one-out, and
include entity resolution in example provenance.

## 4. Medium: exported membership metadata does not identify its rows

`analysis_metadata` stores `reference_weights` positionally and
`summarytable` stores entity IDs, but the metadata itself does not store the
ordered entity IDs that bind those weights to records. This is confirmed by the
probe on IDs `[10,20,30]`: metadata contains weights `[0.5,0.25,0.25]` and has no
`entity_ids` key. The graph configuration has the same positional arrangement,
but its uniform-weight restriction avoids an individual-weight ambiguity in
current graph fits. The material example concerns supported nonuniform pairwise
fits and their exported or nested measurement metadata.

Saving the original table and metadata together preserves the association only
if their original row order remains known. Sorting or joining the exported table
can lose that association, especially with equal scores. The existing advice to
join by `entity_id` cannot repair a weight vector that has no recorded ID map.
This is a provenance gap, not proof that the current fitting functions misalign
weights. No numerical error is claimed from this observation.

Repair direction: include ordered entity IDs alongside the normalized weight
vector, or export membership as an explicit ID/weight table. Document whether
the bundle records query rows as zero mass. Acceptance: export, permute/join the
score table, and reconstruct the exact ID-to-weight mapping; include nonuniform
weights and zero queries. A source/input checksum is complementary and does not
replace this map.

## 5. Comparability and feature multiplicity remain modeling choices

The shared interval from zero to one does not make different estimators or
independently fitted cohorts measurements on one absolute scale. For each fit,
a percentile is relative to that method's reference-score distribution; a cohort
score change can occur with unchanged entity features. The guides already say
this clearly. Repairing numerical inconsistencies must not weaken that boundary.

The all-pairs median treats each pair as an aggregation input. Adding a related
feature adds more pair inputs and can change which aspects of the feature schema
drive the summary. This is a declared schema dependence, not a failure of feature
permutation invariance. The OECD/JRC handbook recommends examining indicator
selection, dependence, normalization, weighting and aggregation when assessing
composite-score robustness. Its discussion of correlated indicators motivates
this review; its linear-weight formulas are not formulas for Platypus's nonlinear
median aggregation. [OECD/JRC handbook, sections 1.6–1.7](https://www.oecd.org/content/dam/oecd/en/publications/reports/2008/08/handbook-on-constructing-composite-indicators-methodology-and-user-guide_g1gh9301/9789264043466-en.pdf).

Acceptance: report numerical invariance checks separately from feature/cohort/
method sensitivity; preserve feature schema and pair names with results; avoid
automatic feature deduplication or a new pair-weighting scheme in a numerical
repair. Either would change the declared aggregation target.

## Procedure, sources and limits

1. Read the preserved simulation report, statistics/mathematics/HPC follow-up
   reports, measurement-scope and workflow guides, and source contracts in
   `types.jl`, `pairs.jl`, `estimation.jl`, `graph.jl`, and `outcome.jl`.
2. Trace the link from reference weights through transformation, fitting,
   percentile normalization, summaries and exported metadata. Separate fixed-
   target representation changes from legitimate target changes.
3. Consult the primary OECD/JRC handbook, specifically its weighting and
   robustness sections. Consult Emmott et al.'s author-hosted arXiv record, whose
   benchmark dimensions include anomaly frequency, clustering and feature
   relevance. That supports varied stress-test design, not an estimate of this
   package's real-world failure probability. [Emmott et al., current arXiv
   record](https://arxiv.org/abs/1503.01158).
4. Run the new contract probe on Julia 1.13.2. Its independent kernel expectations
   passed. Save the raw/normalized observations and metadata inspection in
   `validation/consistency/panel/sociology/`, without overwriting prior evidence.
5. Review the panel's exact counterexamples as established findings; this
   panelist did not independently rerun every numerical oracle or attribute
   every one of the original 630 flags. No new numerical defect is claimed.

Reproduce the contract controls with:

```sh
julia --startup-file=no --project=. validation/consistency/panel/sociology/probe.jl
```

Priority is to repair the independently confirmed representation and numerical
identity failures, specify unresolved ordering/support semantics, and strengthen
reference provenance. Then repeat the fixed-target audit before interpreting
broader specification sensitivity. Nothing in this review evaluates whether
human observers perceive the resulting entities as atypical.
