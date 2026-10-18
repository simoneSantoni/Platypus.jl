# Estimation API

For software attribution and BibTeX, see [Citing Platypus.jl](citation.md).

This reference describes the 0.2 development API. Use Julia 1.13.2 for the
default development environment. Read the [Estimation guide](../guide/estimation.md)
for method selection, assumptions, reference weighting, and computational cost.

| Entry point | Input | Result | Scoring target |
|---|---|---|---|
| [`fit_typicality`](@ref) | `PairSet` | `TypicalityResult` | Pairwise neighborhood contrast, density, or depth |
| [`fit_copula`](@ref) | Source-backed rank `PairSet` | `TypicalityResult` | Gaussian or empirical beta copula density |
| [`fit_graph_typicality`](@ref) | `EntitySet` | `GraphTypicalityResult` | Connectivity across all selected features |
| [`fit_outcome`](@ref) | `TypicalityResult`, outcome vector | `OutcomeResult` | Descriptive least-squares association |

Build pair inputs with [`build_pairs`](@ref). All scoring entry points require
at least two entities; pairwise fits also require at least two features and two
positive reference weights. Results retain input row order and entity IDs.

## Numerical contracts

Neighbor membership uses filtered exact squared-distance comparisons with
rational fallback. Both `backend=:kdtree` and `backend=:brute` currently execute
the conservative scan and report `effective_backend=:brute`; candidate work is
quadratic. Source-backed unweighted rank fits preserve doubled integer midranks;
raw/z-score fits use the represented transformed binary coordinates.

Depth uses exact rational 2×2 covariance, exact algebraic rank, and strict
affine support. Omit `depth_rtol`: numeric truncation arguments now raise.
For source-backed `:none` and `:zscore` depth, the solver uses raw coordinates
because invertible standardization leaves the target unchanged. Rank depth uses
exact reference-weighted ranks. Per-pair metadata records `effective_rank`,
`rank_policy`, `support_policy`, `solver`, and `rank_truncation=false`.

Density and copula computations preserve selected algebraic and reflection
identities; kernels, normal quantiles, and final outputs still use floating
arithmetic. No approximate-tie tolerance is applied. Exact decision predicates
are not a promise of exact representation of every real-valued estimator.

## Pairwise scoring

### Arguments and defaults

`fit_typicality(ps; kwargs...)` accepts the following keywords:

| Keyword | Default | Accepted values and behavior |
|---|---|---|
| `method` | `:lof` | `:knn`, `:lof`, `:kde`, or `:depth` |
| `k` | `0` | For kNN/LOF: `0` selects `clamp(round(Int, sqrt(N)), 1, N-1)`; otherwise an integer in `1:N-1`. KDE/depth require `0`. Boundary ties can include more than `k` neighbors. |
| `weights` | `nothing` | Uniform reference by default, or a `ReferenceWeights` in input row order. Nonuniform weights are supported only by KDE/depth. |
| `novelty_quantile` | `0.1` | Finite real number in `[0,1]`; quantile across normalized pair scores |
| `normalize` | `true` | Must remain `true`; use raw-score exports to inspect unnormalized values |
| `bandwidth` | `nothing` | KDE only: automatic choice, a finite positive scalar, or `(h1, h2)`. Rank-space bandwidths must be at most one. Values are in transformed coordinate units. |
| `backend` | `:kdtree` | `:kdtree` or `:brute`; both currently use the exact candidate scan for neighbor methods |
| `parallel` | `false` | Set `true` to score independent pairs using available Julia threads |
| `depth_rtol` | `nothing` | Retired; every non-`nothing` value raises. Depth uses exact algebraic rank. |
| `provenance` | `(;)` | Caller-supplied `NamedTuple`, copied into result metadata |

Weights affect supported marginal transformations, reference fitting, and the
final mid-CDF. Zero-weight rows remain queries. A positive weight that becomes
zero during Float64 conversion or normalization is rejected. KDE omits the
focal row's kernel contribution and renormalizes the remaining reference mass;
the marginal transformation and bandwidth stay fixed at full-reference values.
Depth includes the focal row when it has positive weight.

### Score scales and summaries

[`fit_typicality`](@ref) returns normalized pair scores and their median,
lower quantile, IQR, and complement-of-median summaries. It retains both
`raw_scores` and the explicit `normalization_scores` used for ordering.
[`fit_copula`](@ref) fits Gaussian or empirical beta copula densities and uses
the same summary contract.

Let ``r_{ep}`` be the value in `normalization_scores` and ``w_i`` reference
weights. The common normalization is

```math
s_{ep} =
\frac{\sum_i w_i\,\mathbf{1}(r_{ip}<r_{ep})+
\tfrac12\sum_i w_i\,\mathbf{1}(r_{ip}=r_{ep})}{\sum_i w_i}.
```

The summaries are

```math
C_e=\operatorname{median}_p s_{ep},\qquad
V_e=Q_q(s_{e\bullet}),\qquad
D_e=\operatorname{IQR}_p s_{ep},\qquad
A_e=1-C_e.
```

Here `novelty` is ``V_e``: **smaller means more unusual**. Normalization
is a descriptive reference percentile, not null-model calibration. Neither
`dispersion` nor the lower quantile measures sampling uncertainty.

For Gaussian KDE on real coordinates (`:none` or `:zscore`), `raw_scores`
contains physical log densities while `normalization_scores` contains Float64
ordinal labels of dimensionless leave-one-out kernel masses. The labels encode
order and equality only, not density magnitude. Exact accumulation of represented
exponential terms preserves distinctions that a physical-unit additive log
constant can erase. Consequently, identical raw Float64 logs can have different
percentiles. Other methods, including rank-space KDE, set
`normalization_scores` equal to `raw_scores`.

Metadata records `normalization_basis=:scale_free_kernel_order` or `:raw_score`.
Real-coordinate KDE compares direct and shifted kernel masses using a common
representation
``\exp(o_e)m_e``: ordinary rows have zero offset and direct rational mass; if
any direct kernel underflows, the row is recomputed after shifting its
unweighted kernel exponents, retaining original reference-weight ratios as
exact rationals outside the exponential. A rational square of a represented
half-exponential extends the range when the shifted exponential underflows.
Differing offsets are compared with adaptive directed MPFR bounds
on the log mass ratio, using explicit local precision and rounding. Equal
offsets use rational comparisons. No global precision state is changed.
These decisions certify the represented expressions, not the exact original
Gaussian mixture. Even the extended shifted terms can underflow;
all-negative-infinite kernel exponents yield represented mass zero, tied with
other represented zeros.

```@docs
fit_typicality
fit_copula
TypicalityResult
```

### Copula keywords

`fit_copula` shares `weights`, `novelty_quantile`, `normalize`, and `provenance`
with pairwise scoring. It has no `k`, `bandwidth`, `backend`, or `parallel`
keyword. Pairs must come from `build_pairs(es; transform=:rank)`.

| Keyword | Default | Accepted values and behavior |
|---|---|---|
| `family` | `:gaussian` | `:gaussian` or `:beta` |
| `regularization` | `0.0` | Gaussian correlation shrinkage in `[0,1)`; singular fits require a positive value large enough to resolve singularity. Beta requires zero. |
| `query_clip` | `nothing` | Gaussian only: a finite value in `(0,0.5)` permits clipping zero-weight endpoint queries. Positive-reference endpoints always raise. |

Gaussian fits support nonuniform reference weights but reject constant
normal-score margins. Beta fits require uniform positive weights and untied
margins; fitting is in-sample. These are model restrictions, not missing-value
or numerical-error fallbacks.

### Pairwise result fields

Write `N = result.n_entities` and `P = result.n_pairs`. Matrix rows align with
`entity_ids`; columns align with `pair_names`.

| Field | Type or shape | Meaning |
|---|---|---|
| `pair_scores` | `Matrix{Float64}`, `N × P` | Weighted mid-percentiles in `[0,1]`; higher means more typical |
| `raw_scores` | `Matrix{Float64}`, `N × P` | Negative kth distance, reciprocal LOF, log density, or depth, depending on method |
| `normalization_scores` | `Matrix{Float64}`, `N × P` | Ordering values used by the mid-CDF; real-coordinate KDE uses ordinal kernel-order labels |
| `conventionality`, `novelty`, `dispersion`, `atypicality` | Each `Vector{Float64}`, length `N` | Summaries defined above; a single pair has zero dispersion |
| `pair_names` | `Vector{Tuple{Symbol,Symbol}}`, length `P` | Ordered feature names for each matrix column |
| `entity_ids`, `labels` | Vectors of length `N` | Integer IDs and optional string labels |
| `method`, `k`, `novelty_quantile`, `normalized` | `Symbol`, `Int`, `Float64`, `Bool` | Effective fit settings; copula methods are `:gaussian_copula`/`:beta_copula`, non-neighbor `k` is zero, and `normalized` is true |
| `configuration` | `NamedTuple` | Recorded settings and provenance; use `analysis_metadata` for a deep copy |

The three score matrices occupy `24*N*P` bytes of Float64 storage, excluding
summaries and metadata. They are detached from inputs but remain mutable Julia
arrays. Editing them does not recompute summaries or configuration; refit to
obtain a coherent new result. KDE `-Inf` raw scores can represent zero density.

### Fit and export example

```@example api_estimation
using Platypus

es = EntitySet([0 0 1; 1 2 4; 2 1 0; 3 4 2; 4 3 5; 5 6 3],
               [:x, :y, :z]; ids=[11, 22, 33, 44, 55, 66])
ps = build_pairs(es; transform=:none)
weights = ReferenceWeights([1, 2, 1, 2, 1, 0])
result = fit_typicality(ps; method=:kde, bandwidth=(0.8, 0.8), weights,
                        provenance=(dataset="api-example",))

normalized = scoretable(result)
physical_logs = scoretable(result; raw=true)
ordering = scoretable(result; normalization=true)
metadata = analysis_metadata(result)
@assert size(result.pair_scores) == (6, 3) # hide
@assert size(normalized) == size(physical_logs) == size(ordering) == (18, 4) # hide
@assert metadata.normalization_basis == :scale_free_kernel_order # hide
@assert metadata.entity_ids == es.entity_ids # hide
summarytable(result)
```

Use a rank `PairSet` for copulas. The following fits use the uniform reference;
the earlier `weights` variable is not automatically reused:

```@example api_estimation
rank_pairs = build_pairs(es)
gaussian = fit_copula(rank_pairs; family=:gaussian, regularization=0.05)
beta = fit_copula(rank_pairs; family=:beta)
@assert gaussian.normalization_scores == gaussian.raw_scores # hide
@assert beta.normalization_scores == beta.raw_scores # hide
(gaussian.method, beta.method)
```

## Graph scoring

Graph scores use all features together and return a different result type.
The stationary masses are proportional to weighted degree. If the graph is
disconnected, this selects component mass proportional to component degree;
the stationary distribution is not unique. See the result metadata for the
component count and selection rule.

| Keyword | Default | Accepted values and behavior |
|---|---|---|
| `k` | `0` | Automatic `min(max(2, round(Int, sqrt(N))), N-1)`, or an integer in `1:N-1` |
| `transform` | `:rank` | `:rank`, `:zscore`, or `:none` |
| `weights` | `nothing` | Only uniform positive `ReferenceWeights` are supported |
| `backend` | `:kdtree` | `:kdtree` or `:brute`; effective backend is currently `:brute` |
| `max_edges` | `1_000_000` | Positive integer limiting stored undirected edges; exceeding it raises instead of dropping ties |
| `provenance` | `(;)` | Caller-supplied `NamedTuple` copied into metadata |

`GraphTypicalityResult` contains length-`N` vectors `typicality`, `atypicality`,
`stationary`, and `degree`; a symmetric `N × N` sparse Float64 `adjacency`
matrix; `entity_ids`, `labels`, `n_entities`, and `configuration`. The adjacency
has no stored self-loops. Isolated vertices receive zero stationary mass when
edges exist; an edgeless graph receives uniform mass. Access `configuration`
through `analysis_metadata` when retaining a detached copy.

```@docs
fit_graph_typicality
GraphTypicalityResult
```

```@example api_estimation
graph = fit_graph_typicality(es; k=2, max_edges=100)
@assert analysis_metadata(graph).effective_backend == :brute # hide
@assert size(graph.adjacency) == (6, 6) # hide
summarytable(graph)
```

## Outcome association

This is descriptive OLS. It does not create a historical reference, supply
independent-observation uncertainty, or support causal interpretation.
Outcome values must match the fitted result's row order.

| Keyword | Default | Accepted values and behavior |
|---|---|---|
| `predictors` | `(:conventionality, :novelty, :dispersion)` | Unique selection from these three fields and `:atypicality` |
| `intercept` | `true` | Include an intercept column |
| `interaction` | `false` | Add conventionality × (1 − novelty); both main effects must be selected |

The outcome vector must contain one finite real value per entity. There must
be more observations than coefficients and a full-rank design. Including both
conventionality and atypicality with an intercept is redundant and raises.
Single-pair results have identical conventionality/novelty and zero dispersion,
so the default predictor set is rank deficient for them; select fewer terms.

`OutcomeResult` stores `coefficients` and `terms` in design-column order,
`observed`, `fitted`, `residuals`, `entity_ids`, the `design` matrix, scalar `r2`,
and `configuration`. A constant outcome has `r2=NaN`. The configuration includes
a detached measurement configuration under `measurement`. Outcome results
support `analysis_metadata` and `plot_outcome_fit`, but not the score/summary
table accessors.

```@docs
fit_outcome
OutcomeResult
```

```@example api_estimation
outcome = fit_outcome(result, [1, 3, 2, 5, 4, 8]; predictors=(:atypicality,))
@assert outcome.entity_ids == result.entity_ids # hide
@assert outcome.terms == [:intercept, :atypicality] # hide
(terms=outcome.terms, coefficients=outcome.coefficients, r2=outcome.r2)
```

## Scores, tables, and configuration

The accessors return detached vectors. Tables return detached data, and
`analysis_metadata` returns a deep copy of the configuration. Save configuration
with exported tables; tables alone cannot reconstruct the analysis. Metadata
includes ordered `entity_ids` binding the `reference_weights` vector to rows.
Pairwise fits also retain `source_reference_weights`, the original represented
weights used before normalization. Preserve these mappings before sorting or
joining exports.

Neighbor metadata separates `requested_backend` from `effective_backend`;
`backend` records the effective scan. The current policy is
`distance_policy=:filtered_exact`. The recorded geometry distinguishes exact
reference midranks from source/transformed binary coordinates. Preserve the
dataset, package revision, and dependency environment separately.

For a `TypicalityResult`, `scoretable(result; raw=true)` selects method-specific
raw scores instead of normalized scores. `scoretable(result; normalization=true)`
exports the values used in the mid-CDF; the two options are mutually exclusive.
For real-coordinate KDE these are ordinal labels, not densities, and their
magnitudes cannot be compared across pairs/fits. For a `GraphTypicalityResult`,
`raw=true` selects stationary mass; `normalization=true` is a pairwise-result
option. Graph results support `atypicality`, `summarytable`,
`scoretable`, and `analysis_metadata`, but have no pairwise novelty or dispersion.

### Export schemas

| Call | Rows | Columns |
|---|---|---|
| `summarytable(pair_result)` | One per entity | `entity_id`, `label`, `conventionality`, `novelty`, `dispersion`, `atypicality` |
| `scoretable(pair_result; ...)` | One per entity–pair combination | `entity_id`, `feature_a`, `feature_b`, `typicality` |
| `summarytable(graph_result)` | One per entity | `entity_id`, `label`, `typicality`, `atypicality`, `stationary`, `degree` |
| `scoretable(graph_result; raw=false)` | One per entity | `entity_id`, `typicality` |

Pair score tables list all pairs for the first entity, then all pairs for the
next. The column name `typicality` stays the same for normalized, raw, and
normalization exports; retain the export mode alongside the table. Graph
`raw=true` places stationary mass in that same column. Join by `entity_id`
(and feature names for pair scores), not by row position after sorting.

```@docs
atypicality
conventionality
novelty
dispersion
summarytable
scoretable
analysis_metadata
```
