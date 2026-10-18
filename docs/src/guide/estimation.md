# Estimation

Choose a method for the quantity you intend to describe. The estimators do not
share a single meaning of typicality: neighborhood isolation, relative local
density, global centrality, fitted dependence, and graph connectivity differ.

The final 84-dataset numerical audit passed 14,637 checks with zero violations
and zero unexpected comparison errors. All 3,701 package assertions and 976
independent mathematical comparator checks passed. Real-coordinate KDE exposes
a separate normalization ordering because physical-unit log adjustments can
erase a genuine density distinction when stored as Float64. The audit checks
the weighted mid-CDF of this exposed ordering; this is an explicit output-contract
change. The numerical limits below still apply beyond the tested scenarios.

## Pairwise methods

```@example estimation
using Platypus
x = collect(range(0.05, 0.95; length=20))
y = x[[3, 8, 1, 17, 12, 6, 20, 10, 4, 15, 2, 18, 9, 5, 13, 19, 7, 16, 11, 14]]
es = EntitySet(hcat(x, y), [:x, :y])
pairs = build_pairs(es)
result = fit_typicality(pairs; method=:lof, k=4,
    provenance=(dataset="permuted-grid-v1",))
summarytable(result)[1:3, :]
```

| Method | Raw score, larger = more typical | Reference weights | Assumptions and controls |
| --- | --- | --- | --- |
| `:knn` | Negative kth-neighbor distance | Uniform positive only | Euclidean geometry; choose `k` |
| `:lof` | Reciprocal Local Outlier Factor | Uniform positive only | Relative local reachability; choose `k`; includes boundary ties |
| `:kde` | Log product kernel density | Nonuniform supported | Smoothing and support assumptions; choose `bandwidth` |
| `:depth` | Mahalanobis depth | Nonuniform supported | Exact rational moments and algebraic rank; strict affine support |
| Gaussian copula | Log fitted copula density | Nonuniform supported | Gaussian dependence; nonconstant margins; explicit singular-fit shrinkage |
| Empirical beta copula | Log empirical beta copula density | Uniform positive only | Untied continuous margins; in-sample empirical beta construction |

Call the first four through [`fit_typicality`](@ref); use [`fit_copula`](@ref)
for the last two. Copula methods require `transform=:rank`.

### Exact neighborhood decisions

`k=0` selects rounded ``\sqrt{N}``, clamped to `1:N-1`. Explicit `k` must
lie in that range. Membership uses squared-distance comparisons: outward-rounded
floating intervals certify clear orderings, and overlapping intervals trigger
exact rational arithmetic. Returned distances are rounded separately. LOF
includes every kth-distance boundary tie under the declared geometry, so a
neighborhood can contain more than `k` rows.

Both requested backends, `:kdtree` (the retained default keyword) and `:brute`,
currently use the conservative scan. No KD-tree candidate pruning is active,
because its earlier floating bounds were not certified for the exact predicate.
Metadata records `requested_backend` and `effective_backend=:brute`.
This policy has quadratic candidate work and linear scratch storage rather than
a dense distance matrix.

Source-backed unweighted rank fits use doubled integer midranks, preserving
rank-lattice geometry before unit-square rounding. Raw and z-score fits instead
interpret the represented transformed Float64 coordinates as exact binary
rationals. These are explicit geometry policies. An arbitrary decimal rescaling
can change represented input ties; it is not automatically an exact similarity.

Duplicates need a defined extension of LOF: a zero-reachability cluster receives
LOF one; a point with positive reachability and a zero-reachability neighbor can
receive inverse-LOF zero. No fixed distance epsilon changes the scale of the
data. Inspect sensitivity to `k`, especially with small clusters or many ties.

### Kernel density

For rank coordinates, KDE uses reflected Epanechnikov product kernels on
``[0,1]^2``. For `:none` and `:zscore` it uses Gaussian product kernels on
real coordinates. Reflection handles bounded support; it does not force the
estimated margins to be uniform and does not make the KDE a copula estimator.

`bandwidth` accepts a positive scalar or a two-element tuple. Rank-space
bandwidths must be at most one. Automatic bandwidths use weighted reference
variance and effective sample size with the bivariate ``N_{\rm eff}^{-1/6}``
rate, and an explicit kernel-scale adjustment. Centered moments use exact
rational arithmetic on represented coordinates and original reference weights;
square roots and kernel values remain floating calculations. This is a smoothing
heuristic, not a guarantee of an optimal bandwidth.

Each query excludes its own kernel contribution, but ranks and bandwidths remain
fixed at their full-reference values. This is not refitting the whole pipeline
leave-one-out. A constant reference coordinate is dropped; a query outside that
constant support receives zero density. Compact support can also produce zero
density, represented by `-Inf` in `raw_scores`. Normalized scores remain finite.
When both reference coordinates are constant, the on-support log density is
returned as exactly zero after validating remaining reference mass. Coordinate
exchange is evaluated symmetrically. Rank-space reflection uses exact reference
ranks and a canonical evaluation of reflected kernel distances.

#### Physical density and normalization ordering

For `transform=:none` or `:zscore`, the Gaussian KDE returns physical log
densities in `raw_scores`. It derives percentile order from dimensionless
kernel masses before subtracting the common log bandwidth/normalization constant.
For rows whose direct exponentials remain representable, with computed Float64
terms ``t_{ej}``, the mass is

```math
a_e = \frac{\sum_{j\ne e}w_j t_{ej}}{\sum_{j\ne e}w_j}.
```

The sums and reference-weight ratios use exact arithmetic on the represented
terms. Dense Float64 ordinal labels of these masses are stored in
`normalization_scores`. Labels preserve the computed ordering and ties; they
are not density values and their magnitude cannot be compared between pairs or
fits. Adding a large physical-unit log constant can make two raw Float64 log
densities equal while these labels remain distinct. Different percentiles for
such raw ties are therefore intentional and visible in the exported ordering.

Real-coordinate KDE compares direct and shifted kernel masses through one
common representation, ``\exp(o_e)m_e``. Unequal weights can make an underflowed
kernel contribution outweigh a positive direct contribution, so underflow is
not itself an ordering rule. Independent oracle controls cover these mixed
cases, including very unequal reference weights.

Ordinary rows have offset zero and exact direct mass. If **any** direct kernel
underflows, the whole row is recomputed by shifting its **unweighted kernel
exponents** by their maximum. Original reference-weight ratios remain exact
rational factors outside the exponential; weighted represented terms are
summed exactly and divided by exact remaining reference mass. When a shifted
exponential underflows, squaring its represented half-exponential as a rational
extends the range without a floating product. This includes cases where a
missing term can outweigh the row's positive direct contributions.

Equal offsets permit exact rational mass comparisons. Otherwise, directed MPFR
bounds compare the log mass ratio with the exact offset difference, doubling
local precision from 128 bits until the intervals separate. Exact special cases
are handled first. Precision and rounding are explicit per value; caller/global
precision state is not changed. The comparison certifies the represented
``\exp(o_e)m_e`` expressions, **not** exact original Gaussian densities.
The extended shifted terms may still underflow relative to the largest kernel
term; if all kernel exponents are negative infinity, the represented mass is
zero and ties with other represented zeros. No epsilon or approximate-tie bucket is used.
Rank-space KDE and every other pairwise method use raw values for normalization.

```@example estimation
real_kde = fit_typicality(build_pairs(es; transform=:none); method=:kde)
@assert analysis_metadata(real_kde).normalization_basis == :scale_free_kernel_order
ordering_table = scoretable(real_kde; normalization=true)
@assert ordering_table.typicality == vec(permutedims(real_kde.normalization_scores))
ordering_table[1:3, :]
```

```@example estimation
weights = ReferenceWeights(vcat(ones(19), 0.0))
weighted = fit_typicality(pairs; method=:kde, weights=weights,
    bandwidth=(0.2, 0.2))
@assert all(isfinite, weighted.pair_scores)
analysis_metadata(weighted).effective_sample_size
```

### Mahalanobis depth

Depth is global center-outward centrality, not local density. Weighted
reference moments include the focal entity when its weight is positive.
The two-dimensional mean, covariance, determinant, and depth are evaluated in
exact rational arithmetic, then depth is converted once to Float64. Exact
algebraic rank determines whether the support is a plane, line, or point.
Exactly off-support queries receive zero; a point reference has depth one on
its support. Very thin but full-rank data remains rank two.

**The spectral truncation policy has been removed.** Omit `depth_rtol`;
supplying a numeric value raises an actionable error rather than silently
ignoring it. Per-pair metadata records the solver, effective rank, exact support
policy, and absence of truncation.

For source-backed `:none` and `:zscore` depth, fitting uses original source
coordinates: invertible standardization does not change Mahalanobis depth.
Rank fits use exact weighted reference ranks. This preserves affine identities
when coordinate changes preserve represented inputs exactly; it cannot recover
information lost before the data reaches the package. Exact equal rational
depths round alike, but differences below Float64 output resolution may still
round to the same stored value. No tolerance merges nearby scores.

### Copula models

```@example estimation
gaussian = fit_copula(pairs; family=:gaussian)
beta = fit_copula(pairs; family=:beta)
@assert size(gaussian.pair_scores) == size(beta.pair_scores)
gaussian.method, beta.method
```

The Gaussian model estimates correlation from weighted normal scores; this is
not a claim to maximum-likelihood estimation. Covariance is evaluated with
exact rational centered moments of the represented normal scores, so a
factorized product reference retains exactly zero covariance without discarding
genuinely weak dependence. Complementary ranks use a common normal-quantile
evaluation with an explicit sign. Constant reference margins are rejected. Singular dependence is rejected by default; choose
`regularization` in `(0,1)` to explicitly shrink correlation toward zero.
Record and justify that modeling change.

Zero-weight queries outside a reference margin can map to rank endpoints.
Gaussian normal scores are infinite there: they are rejected unless
`query_clip` in `(0,0.5)` explicitly clips those query endpoints. Positive-weight
reference endpoints are always rejected. Clipping is a query convention, not a
measurement of extreme-tail probability.

The empirical beta copula is fitted in-sample and requires equal strictly
positive weights and untied margins. It rejects `regularization` and
`query_clip`. Both families return log densities, then use the same percentile
summaries as the other pairwise methods. Beta evaluation reflects complementary
rank/shape configurations canonically. These rules preserve specific proved
symmetries; they do not certify all equalities or orderings of arbitrary
transcendental density values.

## Weights, normalization, and summaries

KDE, depth, and Gaussian copulas apply weights to feature transformations,
reference fitting, and score normalization. Original represented weight ratios
are preserved through exact rational ranks, moments, and mid-CDF mass
accumulation rather than repeatedly rounded normalization. Zero-weight rows are
queries only. At least two rows must have positive reference mass.
Multiplying supplied Float64 weights by a nonbinary factor can change their
represented ratios; that input change is distinct from a broken exact identity.
LOF, kNN, and graph scoring reject nonuniform weights instead of silently ignoring
them.

Every pair's `normalization_scores` are converted to reference-weighted
mid-percentiles before aggregation.
The median is `conventionality`; the quantile selected by
`novelty_quantile` (default 0.1) is `novelty`; the IQR is `dispersion`.
Atypicality equals one minus conventionality. **Lower `novelty` means more
unusual lower-tail typicality.** It is not the minimum or historical newness.

`normalize=false` is no longer supported. Use `result.raw_scores` or
`scoretable(result; raw=true)` to inspect method-specific scores.
`scoretable(result)` returns normalized values.
`scoretable(result; normalization=true)` exports the values that define the
mid-CDF groups. The `raw=true` and `normalization=true` options are mutually
exclusive. Metadata's `normalization_basis` is `:scale_free_kernel_order` for
real-coordinate KDE and `:raw_score` otherwise. Save `analysis_metadata(result)`
with tables; it includes transforms, normalized reference weights, effective
sample size, effective `k`, pair parameters, and caller-supplied provenance.
Ordered `entity_ids` bind weights to rows, and `source_reference_weights`
preserves the original represented weights. Keep that mapping when sorting or
joining score tables. Neighbor fits record requested and effective backends,
geometry, and distance policy; depth parameters record exact rank/support.
Preserve the data, source revision, and Julia environment separately.

Normalization ranks the exposed `normalization_scores`, grouping only exact
equality. For real-coordinate KDE these are ordinal labels, so raw log-density
equality does not determine percentile equality. Preserving selected identities
does not turn every Float64 kernel, quantile, or exponential into an exact real
number. Genuine close-score controls remain; there is no blanket approximate-tie
rule.

## Full-feature graph scoring

```@example estimation
graph = fit_graph_typicality(es; k=4, max_edges=1000)
@assert graph isa GraphTypicalityResult
summarytable(graph)[1:3, :]
```

The graph uses all selected features together. It forms the union of exact
tie-inclusive nearest-neighbor edges and assigns self-tuning affinities
``\exp(-d_{ij}^2/(\sigma_i\sigma_j))``, where each local scale is the kth
distance. Coincident points have affinity one; other pairs involving a zero
scale have affinity zero. The squared geometric ratio is evaluated canonically
as an exact rational expression before its square root and exponential are
rounded. Thus equal geometric expressions share an evaluation, while affinity,
degree, and stationary outputs remain floating values. Both requested backends
use the same filtered exact scan described above.

Weighted degree divided by total degree supplies a stationary distribution of
the undirected random walk. For disconnected graphs this is a deliberate
selection among multiple stationary distributions, assigning component mass
proportional to component degree. Isolated vertices have zero mass if any edges
exist; an edgeless graph uses uniform mass. The percentile of stationary mass
is graph typicality, and its complement is atypicality.

These are connectivity scores, not pairwise summaries or anomaly probabilities.
Self-tuning can make a sparse cluster well connected internally. The graph
representation can retain higher-order structure, without guaranteeing its
detection. `GraphTypicalityResult` therefore has no conventionality/novelty
summary triple. Inspect `analysis_metadata(graph).n_components`.

## Descriptive outcome association

```@example estimation
outcome = collect(1.0:20.0)
association = fit_outcome(result, outcome; predictors=(:conventionality,))
@assert length(association.fitted) == 20
association.terms, association.coefficients
```

[`fit_outcome`](@ref) fits descriptive OLS with column-pivoted QR. Outcome
values must be finite and match result row order. Select measurement features
before estimation; the outcome must not be one of them.

Predictors can be `:conventionality`, `:novelty`, `:dispersion`, and
`:atypicality`. The defaults use the first three with an intercept.
`interaction=true` adds conventionality × (1 − novelty), requiring both main
effects. Constant or redundant designs are rejected; an intercept plus both
conventionality and atypicality is necessarily redundant. With only one pair,
several default predictors are redundant, so select a smaller design explicitly.

This routine does not construct historical reference windows, establish
causality, or return independent-observation standard errors. Generated scores
share a reference sample; study-specific uncertainty requires a resampling
design that refits measurement. See the
[Research workflow](research_workflow.md).

## Computational cost

Let ``N`` be the number of entities, ``p`` the number of features, and
``P=p(p-1)/2``. Pairwise results retain three Float64 score matrices:
`raw_scores`, `normalization_scores`, and `pair_scores`. They require
approximately `24*N*P` bytes before metadata, tables, and temporary storage.

| Stage | Working storage | Main cost and limitations |
| --- | --- | --- |
| Pair construction | Shared feature vectors plus raw snapshot, ``O(Np+P)`` | Rank transforms sort features |
| kNN/LOF | Linear scratch storage per active pair | Quadratic candidate scans; floating filters with rational fallback and neighbor reduction costs |
| KDE | Reference/query buffers plus rational moments | All query/reference contributions; sorting and exact-arithmetic allocation costs |
| Depth / Gaussian copula | Moment accumulators and score arrays | Exact rational centered moments; cost depends on integer bit lengths, not only row count |
| Empirical beta copula | Linear buffers per pair beyond scores | All query/reference contributions, with canonical reflection and sorted reductions |
| Graph | Features plus stored sparse edges | Quadratic candidate work; ``O(Nk)`` edges only when ties are limited |

No production neighbor method allocates a dense ``N\times N`` distance matrix.
Both requested backends currently scan candidates; do not assume a KD-tree
speedup. Rational arithmetic has variable-size integers and allocation costs,
so the score-matrix formula is not a total-memory estimate. `max_edges` bounds
retained undirected edges, not candidate work, temporary exact arithmetic,
dictionary overhead, or peak memory. High-dimensional graphs can remain expensive.

For pairwise `fit_typicality`, `parallel=true` distributes pairs over Julia
threads. Start Julia with `--threads=auto` or a chosen thread count; each active
pair needs its own working memory. Creating `scoretable` also materializes
`N*P` rows. Benchmark the whole workflow at the intended size, not only a single
neighbor query.
