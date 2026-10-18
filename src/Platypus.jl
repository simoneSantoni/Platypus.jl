"""
    Platypus

Descriptive entity atypicality relative to a supplied collection and feature
schema. Pairwise neighborhood, bounded-density, depth and copula scores expose
both raw components and normalized summaries. A separate full-feature graph
estimator measures connectivity. Outcome fitting is descriptive OLS, and plot
functions provide backend-independent recipes.

The pairwise median/lower-tail design is inspired by Uzzi et al. (2013); IQR
is an additional summary. Percentile normalization is not null calibration.
See the manual and docs/research for method assumptions, validation and changes.
"""
module Platypus

using CSV
using DataFrames
using Printf
using LinearAlgebra
using SparseArrays
import NearestNeighbors
import Distributions
using RecipesBase
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
export analysis_metadata, OutcomeResult, TypicalityPlot, GraphTypicalityResult

# Visualization
export plot_feature_distributions, plot_pair_cloud, plot_correlation_heatmap
export plot_pair_redundancy, plot_typicality_profile
export plot_atypicality_distribution, plot_novelty_conventionality
export plot_outcome_fit, plot_pair_heatmap

include("types.jl")
include("entities.jl")
include("numerics.jl")
include("pairs.jl")
include("eda.jl")
include("neighbors.jl")
include("density.jl")
include("depth.jl")
include("estimation.jl")
include("graph.jl")
include("outcome.jl")
include("visualization.jl")

end # module
