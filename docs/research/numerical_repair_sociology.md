# Numerical repair acceptance: reference identity and measurement contracts

**Final integration addendum:** the lead reran the saved 33-check acceptance
script on final source `a60109d4a70036f7b7474d7a410e105b17f182f869b218332a45e5a0acb05f8c`
with Julia 1.13.2/two threads: 33/33 passed. The final full package run includes
all 130 contract checks and passes 3,701 assertions. Updated execution/provenance
records are saved with the probe; prior `c7c96e` records are retained under
`repair_iterations/mixed_underflow_interim/panel_sociology`. The historical
review stages below remain identifiable. The final broad audit passes all
14,637 checks with zero flags and zero unexpected errors; its closure and
explicit KDE CDF contract change are in the [integration record](numerical_repair.md).

Reviewed 2026-10-18 on Julia 1.13.2. Scope is internal consistency, reference
identity, provenance and preservation of declared measurement targets.
Human-perception alignment is excluded. Source was treated as frozen; this
review changes only its report and acceptance artifacts. Earlier artifacts
remain preserved separately.

## Decision

**Accept the repaired identity/reference binding and pair-collection contracts
within the checks performed here.** No new source-level blocker was found in
this scope. The repair preserves the supplied-collection estimand and the
separation of local contrast, density, centrality, copula density and graph
connectivity. Exact algebraic depth rank is an explicitly changed numerical
support policy, rather than a claim that the old spectral policy is unchanged.

At this review stage, full numerical acceptance remained the integration lead's
decision after the 84-dataset audit and complete package checks. During this review,
`validation/consistency/repaired_results/` was still accumulating dataset files;
its final check/summary tables were not available. This report does not infer
that audit's outcome from a pilot or from the narrower tests below.

The review identified one stale public sentence in the measurement-scope table:
depth was still described as having a spectral singular-support convention.
The integration lead corrected it to exact algebraic rank and strict affine
support; I reread and confirmed that correction. Numeric `depth_rtol` is rejected
as documented in the source and estimation guide. The lead also reports the
complete package suite at **3,616/3,616** and a clean documentation build; these
are first-iteration integration results, separate from the checks I executed below. The final
build after all report updates and the full audit remain owned by integration.

## Final extension: explicit KDE normalization ordering

The acceptance extension reviewed frozen production source SHA256
`c7c96e892949130fa15a5b859a4f7c8eb1413e5e0e5529b26c7ad450c2c690c5`.
This digest was recomputed by concatenating sorted `src/` file contents with
newline separators, matching the integration convention. I read the updated
main repair and statistical reports, the actual density/summary implementation,
new public method contracts, and the extended structural tests.

**Accept the explicit normalization representation within this review's scope.**
Real-coordinate KDE still estimates the same declared kernel density relative
to the same row-level reference and weights. Its floating physical log-density
export can lose a small real ordering when a unit change adds a large common
constant. The new `normalization_scores` exposes a separate ordering derived
from dimensionless kernel masses before that lossy adjustment. Consequently,
the implementation preserves the density-ranking target while intentionally
changing the old contract that ranked the stored physical Float64 logs.

For real-coordinate KDE, the new values are dense ordinal labels, not densities,
probabilities or distances between entities. Their gaps and magnitudes have no
measurement interpretation across pairs or fits. They encode the computed order
and equality only. Equal physical raw logs can legitimately receive different
percentiles when these labels resolve a genuine contrast; consumers rebuilding
scores must use `normalization_scores` with the ID-bound source reference mass,
not rerank `raw_scores`. `scoretable(result; normalization=true)` and recorded
`normalization_basis=:scale_free_kernel_order` make this inspectable. Other
pairwise methods, including rank-space KDE, retain `:raw_score` normalization.

A zero-mass query can introduce a new ordinal level and renumber existing labels
without changing the positive reference or its percentiles. Requiring literal
label stability across such fits would be an invalid extra invariance: preserve
the relative ordering and reference-weighted mid-CDF, not arbitrary label
magnitudes. The extended contract tests explicitly retain the zero-query
percentile control. No category, audience, historical-newness or inferential
interpretation is introduced by the new representation.

The underlying exponential evaluations and underflow fallback remain finite
precision. Exact summation of their represented masses does not establish every
transcendental ordering. This boundary is stated in the source and public guide;
it does not justify treating ordinal keys as exact physical density estimates.
The separate depth support-policy change and the first full audit's residual
flags also remain explicitly identified, rather than being retrospectively
described as an unchanged fully passing contract.

The saved 33-check acceptance script was rerun unchanged on this frozen source:
**33/33 passed**. The extended structural/normalization suite passed
**130/130** on the same source; its result is recorded in the linked
[`numerical_contracts.log`](../../validation/consistency/repair_panel/sociology/numerical_contracts.log).
The exact commands, source digest and check counts are also retained in
[`final_provenance.json`](../../validation/consistency/repair_panel/sociology/final_provenance.json).
This review introduces no new source-level blocker. Final broad-suite/audit
counts belong to the integration report and must refer to this final source,
not the earlier 3,616-check iteration.

## Identity and reference binding

`src/estimation.jl`, `_summarize`, now records ordered `entity_ids` alongside
`reference_weights` and `source_reference_weights`. The latter preserves the
original represented weights used for exact reference geometry and percentile
mass; the former is their normalized floating export. This resolves the prior
positional-metadata gap: a sorted or joined score table can reconstruct each
entity's membership, including a zero-mass query, by ID.

`analysis_metadata` returns a deep copy. Results copy IDs, labels and pair names;
metadata mutation through the accessor does not rewrite fitted provenance.
`fit_outcome` nests the complete measurement metadata, including the identity
map. `src/graph.jl` separately records its ordered entity IDs and uniform
reference weights. These are analysis snapshots, not a certification that a
user-supplied ID or dataset declaration is truthful. Direct mutation of public
result fields and concurrent mutation of fitting inputs remain outside this
snapshot guarantee.

I independently reran `test/numerical_contracts.jl`: **91/91 checks passed**.
It covers nonuniform mass joined through unsorted IDs, detached metadata,
nested outcome identity, graph identity, requested/effective backend reporting,
complete pair coverage and legitimate pair reordering, among other numerical
and indexing controls.

I also ran **33/33 additional checks** using IDs `[30,10,40,20,99]`,
weights `[1,2,3,4,0]`, and coordinates
`[0 3; 1 0; 2 2; 4 1; 1.5 1.5]`. For KDE, depth and the Gaussian copula:

- Sort the score table, join metadata by ID, then sort the joined table. The
  recovered original masses are `[2,4,1,3,0]` and normalized masses are
  `[0.2,0.4,0.1,0.3,0]`; query ID 99 retains zero mass.
- Mutate the returned metadata's ID and source-weight vectors; fresh metadata
  still contains the original identity and mass.
- Change IDs and descriptive labels while preserving row/coordinate/weight
  alignment; raw and normalized scores are unchanged.
- Check pair names, complementary atypicality and the default lower-tail
  quantile's ordering relative to conventionality.

The first inline run incorrectly assumed a DataFrames join would preserve
left-table row order, producing six assertion failures in the probe. Sorting
**after** the join corrected that test assumption; all 33 checks then passed
without changing package code. This illustrates why identity-based joins must
not be replaced by a positional comparison.

## PairSet integrity and preserved targets

`src/types.jl`, `_validate_pairset`, requires the combinatorial number of pairs,
unique valid `i < j` index tuples, and matching feature names. Together these
conditions ensure every unordered pair appears once. Missing or duplicated
projections can no longer silently change the inputs to median/quantile/IQR
aggregation. Arbitrary order of the complete pair collection is still valid;
the passing permutation controls demonstrate that the repair does not impose an
irrelevant canonical column order.

Known transformations require source provenance, and source-backed pairs require
a known transformation. Stored pair coordinates must agree with the declared
source transformation. Legacy manually supplied coordinates remain explicitly
`:unknown`; they cannot establish the raw-feature provenance needed by weighted
preprocessing or copula fitting. Validation cannot infer a user's intended
meaning of a feature or repair an intentionally misassigned identifier, but it
now checks the structural relationships the API promises.

The specialist repairs preserve the target distinctions recorded in the manual:

| Repair | Contract consequence |
| --- | --- |
| Exact kth-boundary decisions and conservative scan | Retains tie-inclusive neighborhoods; actual backend is recorded even when `:kdtree` was requested. Quadratic candidate work is disclosed. |
| Exact source-rank and reference-mass arithmetic | Preserves represented reference geometry through preprocessing and percentiles; avoids introducing an epsilon-based definition of score ties. |
| KDE and copula symmetry/cancellation repairs | Keeps the existing density definitions and restricted beta family; no new estimator is substituted merely to eliminate a failed identity. |
| Exact algebraic depth and strict support | Deliberately replaces tolerance-based spectral truncation; thin full-rank inputs remain full rank and the obsolete numeric tolerance is rejected. |
| Pair coverage and export identity | Preserves equal participation of the declared pair set and retains the mapping from reference membership to entity identity. |

The normalized output is a descriptive weighted mid-CDF of explicit
`normalization_scores`. Except for real-coordinate KDE these equal computed raw
scores; the exception and its changed contract are discussed above. Summaries
remain the median, selected quantile, IQR and complementary
median. `novelty` still increases toward typicality, so lower values mean greater
lower-tail unusualness. IQR remains spread across pairs rather than a standard
error. Neither exact arithmetic in selected operations nor a shared percentile
scale adds null calibration, historical newness, or cross-method equivalence.

## Record granularity remains a specification choice

The previous sociology control is retained conceptually: splitting one
positive-mass entity into identical positive-mass rows changes kernel-only
leave-one-out exclusion and effective sample size. This is different from
appending a zero-mass query and is not an invariance the repair should impose.

The additional inline run rechecked the analytic Gaussian-kernel control with
coordinates `(0,0), (1,1), (2,0)`, masses `[2,1,1]`, and fixed bandwidth one.
Splitting the first row into two rows of mass one changes its expected density
from `[K(1,1)+K(2,0)]/2` to `[K(0,0)+K(1,1)+K(2,0)]/3`, where
`K(x,y)=exp(-(x*x+y*y)/2)/(2*pi)`. Both numerical log densities agree with these
independent formulas, and the resulting focal values differ. These are the
remaining three checks in the 33-check inline set.

The initial reading found that the guides explained kernel-only leave-one-out
without stating this record-splitting distinction explicitly. The integration
lead added the entity-resolution and row-splitting explanation to the
measurement-scope guide; I confirmed it is present. This completes the earlier
documentation recommendation without changing the estimator. Both documentation
issues raised in this acceptance review are therefore resolved.

## Procedure and evidence limits

I read the mathematics, statistics, Julia, HPC and communication numerical-repair
reports, the earlier sociology follow-up, current types/preprocessing/estimation/
graph/outcome source, and the scope/workflow/estimation guides. The primary
algorithm evidence is recorded by those specialist reports; this acceptance
pass introduces no new literature claim or algorithm.

Reproduce the existing structural/provenance suite from the repository root:

```sh
JULIA_DEPOT_PATH=/tmp/platypus-review-julia-depot:/home/simon/.julia \
julia --startup-file=no --project=. -e \
'using Test, Platypus, DataFrames; include("test/numerical_contracts.jl")'
```

The 33 additional checks were initially executed from standard input. The
corrected checks are now preserved in
[`acceptance.jl`](../../validation/consistency/repair_panel/sociology/acceptance.jl),
with the saved rerun in
[`execution.log`](../../validation/consistency/repair_panel/sociology/execution.log).
Run that same script from the repository root:

```sh
JULIA_DEPOT_PATH=/tmp/platypus-review-julia-depot:/home/simon/.julia \
julia --startup-file=no --project=. \
validation/consistency/repair_panel/sociology/acceptance.jl
```

The first saved script and its execution log were copied without replacement to
[`first_full/panel_sociology`](../../validation/consistency/repair_iterations/first_full/panel_sociology)
before recording the final rerun in the active acceptance directory. No earlier
audit probe or CSV was overwritten. The specialist numerical suites
were inspected as evidence but not all rerun by this panelist. The integration
report must supply the final frozen-source fingerprint, complete suite totals,
full audit comparison and any remaining exceptions before declaring broad
numerical acceptance.
