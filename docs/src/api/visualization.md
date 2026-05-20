# Visualization API

All plotting functions are reserved by name in v0.1.0 and currently raise
`not implemented`. This page documents what each plot is *meant* to
render so that the public API is stable from the package's first release.
Implementations will land as [Plots.jl](https://docs.juliaplots.org)
recipes (or CairoMakie equivalents) in a follow-up revision.

The data inputs each plot consumes are already produced by the
exploratory and estimation functions:

| Plot | Consumes | Produced by |
|------|----------|-------------|
| `plot_feature_distributions` | an `EntitySet` | `load_entities` |
| `plot_pair_cloud` | a `FeaturePair` from a `PairSet` | `build_pairs` |
| `plot_correlation_heatmap` | `pair_correlations(::EntitySet)` | `pair_correlations` |
| `plot_pair_redundancy` | `pair_redundancy(::PairSet)` | `pair_redundancy` |
| `plot_typicality_profile` | one entity's row of `pair_scores` | `fit_typicality` |
| `plot_atypicality_distribution` | `atypicality(::TypicalityResult)` | `fit_typicality` |
| `plot_novelty_conventionality` | `summarytable(::TypicalityResult)` | `fit_typicality` |
| `plot_outcome_fit` | a fitted outcome model | `fit_outcome` (planned) |
| `plot_pair_heatmap` | the `pair_scores` matrix | `fit_typicality` |

Splitting "compute" from "draw" keeps the plotting surface independent of
any particular plotting backend.

## Descriptive plots

### `plot_feature_distributions`

Small-multiples of the marginal distribution of each feature. The first
diagnostic to run: it shows the shape that the copula-rank transform of
[`build_pairs`](@ref Platypus.build_pairs) will flatten to uniform.

```@docs
plot_feature_distributions
```

### `plot_pair_cloud`

Scatter plot of a single feature pair's two-dimensional cloud, with
entities shaded by local typicality. This is the elementary picture
behind the whole method — the cloud in which one entity's position is
scored — and is the natural figure for explaining a particular atypical
combination.

```@docs
plot_pair_cloud
```

### `plot_correlation_heatmap`

Heatmap of the feature-by-feature Spearman correlation matrix, the visual
companion to [`pair_correlations`](@ref Platypus.pair_correlations).
Clusters of strongly correlated features warn that the corresponding
feature pairs are redundant.

```@docs
plot_correlation_heatmap
```

### `plot_pair_redundancy`

Diagnostic plot of the per-pair redundancy weights from
[`pair_redundancy`](@ref Platypus.pair_redundancy): how much each feature
pair is discounted for near-collinearity. Reads directly on the caveat
that the ``n(n-1)/2`` pairwise scores are not an i.i.d. sample.

```@docs
plot_pair_redundancy
```

### `plot_typicality_profile`

The distribution of one entity's ``n(n-1)/2`` pairwise typicality scores,
with its conventionality (median), novelty (lower tail), and dispersion
(IQR) marked. This is the per-entity view that motivates the
[Uzzi et al. (2013)](../references.md) summary triple — the same picture
Uzzi et al. drew for a paper's distribution of co-citation ``z``-scores.

```@docs
plot_typicality_profile
```

## Estimation-output plots

### `plot_atypicality_distribution`

Histogram of per-entity atypicality across the whole collection — the
population view of where the focal entity sits.

```@docs
plot_atypicality_distribution
```

### `plot_novelty_conventionality`

Scatter of `novelty` against `conventionality`, the two-dimensional
summary plane of [Uzzi et al. (2013)](../references.md). Their central
empirical finding — that high-impact work occupies the high-conventionality,
high-novelty corner — is a statement about *position in this plane*, which
makes it the recommended plot for relating typicality to an outcome.

```@docs
plot_novelty_conventionality
```

### `plot_outcome_fit`

Fitted-versus-observed diagnostic for an outcome model produced by
[`fit_outcome`](@ref Platypus.fit_outcome). Because the hypothesised
relationship is non-monotone, the recommended overlay is the fitted
response surface over the novelty–conventionality plane.

```@docs
plot_outcome_fit
```

### `plot_pair_heatmap`

Heatmap of the full ``n_{\text{entities}} \times n_{\text{pairs}}``
pairwise typicality matrix. Useful for spotting whether an entity's
atypicality is concentrated in a few feature pairs or spread evenly —
the structure the `dispersion` summary quantifies.

```@docs
plot_pair_heatmap
```

## References

The visualization conventions follow [Uzzi et al. (2013)](../references.md)
for the conventionality–novelty plane and
[Breunig et al. (2000)](../references.md) for the local-density view of
the pair cloud. Full citations are on the [References](../references.md)
page.
