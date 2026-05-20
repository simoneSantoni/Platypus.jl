"""
    Platypus.jl - Entity atypicality in large, heterogeneous sets

A Julia package for measuring how *atypical* an entity is within a large
and heterogeneous collection — a setting in which there is no clean,
externally given reference set to compare the entity against.

Atypicality is built from the bottom up out of *pairs* of entity
features. Each pair of features defines a two-dimensional cloud of
points; the focal entity's position in that cloud is scored for *local*
typicality; and the distribution of those `n(n-1)/2` pairwise scores is
summarised, per entity, into a small number of statistics —
*conventionality* (central tendency), *novelty* (lower tail), and
*dispersion*. This is the design introduced by Uzzi, Mukherjee,
Stringer & Jones (2013), "Atypical Combinations and Scientific Impact",
generalised from discrete co-citation counts to continuous features.

The public surface follows the layout of `VacancyChains.jl`: domain
types (`Entity`, `EntitySet`, `FeaturePair`, `PairSet`,
`ReferenceWeights`), verb-named entry points (`load_entities`,
`build_pairs`, `fit_*`, `plot_*`), and a per-feature file layout
(`entities.jl`, `pairs.jl`, `eda.jl`, `estimation.jl`,
`visualization.jl`).
"""
module Platypus

using CSV
using DataFrames
using Printf
import Statistics: mean, std, median, quantile, cor
import StatsBase: tiedrank, iqr, corspearman

# Core types
export Entity, EntitySet, FeaturePair, PairSet, ReferenceWeights
export TYPICALITY_METHODS
export features, label, feature_names
export n_entities, n_features, n_pairs

# Data organisation
export load_entities
export build_pairs

# Exploratory data analysis
export feature_summary, rank_matrix, pair_correlations, pair_redundancy
export nearest_entities, compute_entity_statistics, pair_table

# Estimation
export TypicalityResult
export fit_typicality, fit_graph_typicality, fit_outcome, fit_copula
export atypicality, conventionality, novelty, dispersion
export scoretable, summarytable

# Visualization
export plot_feature_distributions, plot_pair_cloud, plot_correlation_heatmap
export plot_pair_redundancy, plot_typicality_profile
export plot_atypicality_distribution, plot_novelty_conventionality
export plot_outcome_fit, plot_pair_heatmap

include("types.jl")
include("entities.jl")
include("pairs.jl")
include("eda.jl")
include("estimation.jl")
include("visualization.jl")

end # module
