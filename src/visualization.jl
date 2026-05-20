"""
Visualization stubs.

Plot functions are exported under their final names but currently raise
`not implemented` errors. They are intended to be implemented as
[Plots.jl](https://docs.juliaplots.org) recipes (or CairoMakie
equivalents) in a follow-up revision; reserving the names here keeps the
public API stable from v0.1.0 onward.

Each plot consumes a data shape that the exploratory and estimation
functions already produce — [`feature_summary`](@ref),
[`pair_correlations`](@ref), [`pair_redundancy`](@ref),
[`scoretable`](@ref), and [`summarytable`](@ref) — so the plotting
surface stays independent of any particular backend.
"""

# Descriptive plots ----------------------------------------------------------

"""
    plot_feature_distributions(es::EntitySet)

Small-multiples of the marginal distribution of each feature. Not yet
implemented.
"""
plot_feature_distributions(::EntitySet) =
    throw(ErrorException("plot_feature_distributions is not implemented yet"))

"""
    plot_pair_cloud(ps::PairSet, pair)

Scatter plot of one feature pair's two-dimensional cloud, with the focal
entity (or all entities) shaded by local typicality. `pair` selects the
[`FeaturePair`](@ref) by index or by `(feature_a, feature_b)` name tuple.
Not yet implemented.
"""
plot_pair_cloud(::PairSet, ::Any) =
    throw(ErrorException("plot_pair_cloud is not implemented yet"))

"""
    plot_correlation_heatmap(es::EntitySet)

Heatmap of the feature-by-feature Spearman correlation matrix — the
visual companion to [`pair_correlations`](@ref). Not yet implemented.
"""
plot_correlation_heatmap(::EntitySet) =
    throw(ErrorException("plot_correlation_heatmap is not implemented yet"))

"""
    plot_pair_redundancy(ps::PairSet)

Diagnostic plot of per-pair redundancy weights from
[`pair_redundancy`](@ref): how much each feature pair is discounted for
near-collinearity. Not yet implemented.
"""
plot_pair_redundancy(::PairSet) =
    throw(ErrorException("plot_pair_redundancy is not implemented yet"))

"""
    plot_typicality_profile(r::TypicalityResult, entity)

The distribution of one entity's `n(n-1)/2` pairwise typicality scores,
with its conventionality, novelty, and dispersion marked. `entity`
selects the entity by row index or by entity ID. This is the
per-entity view that motivates the Uzzi-style summary triple. Not yet
implemented.
"""
plot_typicality_profile(::TypicalityResult, ::Any) =
    throw(ErrorException("plot_typicality_profile is not implemented yet"))

# Estimation-output plots ----------------------------------------------------

"""
    plot_atypicality_distribution(r::TypicalityResult)

Histogram of per-entity atypicality across the whole collection. Not yet
implemented.
"""
plot_atypicality_distribution(::TypicalityResult) =
    throw(ErrorException("plot_atypicality_distribution is not implemented yet"))

"""
    plot_novelty_conventionality(r::TypicalityResult)

Scatter of `novelty` against `conventionality`, the two-dimensional
summary plane of Uzzi et al. (2013). Not yet implemented.
"""
plot_novelty_conventionality(::TypicalityResult) =
    throw(ErrorException("plot_novelty_conventionality is not implemented yet"))

"""
    plot_outcome_fit(result)

Fitted-versus-observed diagnostic for an outcome model produced by
[`fit_outcome`](@ref). Not yet implemented.
"""
plot_outcome_fit(::Any) =
    throw(ErrorException("plot_outcome_fit is not implemented yet"))

"""
    plot_pair_heatmap(r::TypicalityResult)

Heatmap of the full `n_entities × n_pairs` pairwise typicality matrix.
Not yet implemented.
"""
plot_pair_heatmap(::TypicalityResult) =
    throw(ErrorException("plot_pair_heatmap is not implemented yet"))
