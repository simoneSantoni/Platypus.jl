# Platypus.jl

[![Atypicality Analysis](https://img.shields.io/badge/Atypicality-Analysis-orange.svg)](https://github.com/simoneSantoni/Platypus.jl)
[![Build Status](https://github.com/simoneSantoni/Platypus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/simoneSantoni/Platypus.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://simoneSantoni.github.io/Platypus.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://simoneSantoni.github.io/Platypus.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.9+-purple.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
  <img src="docs/src/assets/logo.png" alt="Platypus.jl icon" width="160">
</p>

A Julia package for measuring **entity atypicality** in large, heterogeneous sets — and, in particular, in sets whose reference population has no clean boundary.

## Overview

How atypical is one entity relative to its peers — a song in a market with fuzzy genre boundaries, a firm in a loosely defined industry, a paper in a sprawling literature? Platypus.jl answers this by decomposing each entity into the `n(n-1)/2` *pairs* of its features and scoring the local typicality of each pairwise combination, following the design of Uzzi et al. (2013), *Atypical Combinations and Scientific Impact*. The package provides tools for:

- **Data organization**: ingesting entity records and projecting them onto feature pairs
- **Descriptive analysis**: feature summaries, correlation/redundancy diagnostics, nearest neighbours
- **Visualization**: feature distributions, pair clouds, typicality profiles
- **Estimation**: pairwise local-density atypicality with conventionality / novelty / dispersion summaries
- **Inference visualization**: atypicality distributions, the novelty–conventionality plane, outcome fits

When the reference set is fuzzy, Platypus.jl makes the reference *endogenous* (the collection is its own reference) and *local* (typicality is scored against local density, not a global baseline).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/simoneSantoni/Platypus.jl")
```

## Features

### 1. Data Organization

Functions for ingesting entity records and projecting them onto feature pairs.

```julia
load_entities(df)                    # Load entities from a DataFrame
load_entities("entities.csv")        # Load from CSV file
EntitySet(matrix, feature_names)     # Collection of entities sharing a schema
build_pairs(es; transform=:rank)     # Project onto all n(n-1)/2 feature pairs
ReferenceWeights(weights)            # Graded membership in the reference set
```

### 2. Exploratory Data Analysis

Summary statistics over features, pairs, and entities.

```julia
feature_summary(es)                  # Per-feature descriptive statistics
rank_matrix(es)                      # Copula-rank transform of the data matrix
pair_correlations(es)                # Spearman correlation of every feature pair
pair_redundancy(ps)                  # Per-pair redundancy weights
nearest_entities(es, idx; k=5)       # Nearest neighbours of an entity
compute_entity_statistics(es)        # Per-entity descriptive statistics
```

### 3. Descriptive Visualization

Plot recipes for features and feature-pair structure.

```julia
plot_feature_distributions(es)       # Marginal distribution of each feature
plot_pair_cloud(ps, pair)            # Scatter of one feature pair's cloud
plot_correlation_heatmap(es)         # Heatmap of feature correlations
plot_pair_redundancy(ps)             # Per-pair redundancy diagnostic
plot_typicality_profile(result, e)   # One entity's pairwise score distribution
```

### 4. Estimation

Statistical models for entity atypicality.

```julia
fit_typicality(ps; method=:lof)      # Pairwise local-density atypicality
fit_graph_typicality(es)             # Endogenous graph random-walk typicality
fit_outcome(result, y)               # Regression of an outcome on the summaries
fit_copula(ps)                       # Parametric copula-density typicality
```

### 5. Visualization of Estimation Outcomes

Plot recipes for fitted models.

```julia
plot_atypicality_distribution(result)  # Histogram of per-entity atypicality
plot_novelty_conventionality(result)   # The Uzzi novelty–conventionality plane
plot_outcome_fit(result)               # Fitted vs. observed outcome
plot_pair_heatmap(result)              # Entity × feature-pair typicality heatmap
```

## Usage

### Basic Example

```julia
using Platypus

# An n_entities × n_features matrix of continuous features
es = EntitySet(rand(500, 8),
               [:danceability, :energy, :acousticness, :valence,
                :tempo, :loudness, :speechiness, :liveness])

# Project onto all 8*7/2 = 28 feature pairs (copula-rank transform)
pairs = build_pairs(es)

# Score local typicality on each pair and aggregate per entity
result = fit_typicality(pairs; method=:lof)

# View results
println(result)
```

### Result Structure

`fit_typicality` returns a `TypicalityResult` with:

- `pair_scores::Matrix{Float64}`: `n_entities × n_pairs` local typicality scores
- `conventionality::Vector{Float64}`: Median pairwise typicality, per entity
- `novelty::Vector{Float64}`: Lower-tail (10th percentile) score, per entity
- `dispersion::Vector{Float64}`: Inter-quartile range of pairwise scores
- `atypicality::Vector{Float64}`: `1 .- conventionality`
- `method::Symbol`: Scoring rule used
- `n_entities::Int`, `n_pairs::Int`: Problem dimensions

```julia
# Access results
atypicality(result)        # Headline per-entity atypicality
conventionality(result)    # Median pairwise typicality
novelty(result)            # Lower-tail novelty
summarytable(result)       # Per-entity summary as DataFrame
scoretable(result)         # Long-form entity × feature-pair table
```

### Loading Data

```julia
using DataFrames

# From DataFrame — numeric columns become features
df = DataFrame(
    id           = [1, 2, 3],
    genre        = ["jazz", "rock", "pop"],
    danceability = [0.31, 0.74, 0.66],
    acousticness = [0.88, 0.12, 0.40],
)
es = load_entities(df; id_col=:id, label_col=:genre)

# From CSV file
es = load_entities("entities.csv"; id_col=:id)
```

### Graded Reference Membership

When the reference set has a fuzzy boundary, weight each entity's contribution.

```julia
# One non-negative weight per entity (e.g. by recency or genre proximity)
w = ReferenceWeights(weights)

# Honoured by the :kde and :depth methods
fit_typicality(pairs; method=:kde, weights=w)
```

### Scoring Methods

```julia
fit_typicality(pairs; method=:knn)     # k-th nearest-neighbour density
fit_typicality(pairs; method=:lof)     # Local Outlier Factor (default)
fit_typicality(pairs; method=:kde)     # Gaussian kernel density
fit_typicality(pairs; method=:depth)   # Mahalanobis data depth
```

## Utility Functions

```julia
# EntitySet accessors
n_entities(es)                       # Number of entities
n_features(es)                       # Number of features
feature_names(es)                    # Feature names

# PairSet accessors
n_pairs(ps)                          # Number of feature pairs

# Entity accessors
features(entity)                     # Feature vector
label(entity)                        # Category label (or nothing)
```

## Running Tests

```julia
include("test/runtests.jl")
```

## Documentation

For more detailed documentation, see:

- [Stable Documentation](https://simoneSantoni.github.io/Platypus.jl/stable/)
- [Development Documentation](https://simoneSantoni.github.io/Platypus.jl/dev/)

## References

1. Uzzi, B., Mukherjee, S., Stringer, M., & Jones, B. (2013). Atypical combinations and scientific impact. *Science*, 342(6157), 468–472.

2. Askin, N., & Mauskapf, M. (2017). What makes popular culture popular? Product features and optimal differentiation in music. *American Sociological Review*, 82(5), 910–944.

3. Breunig, M.M., Kriegel, H.-P., Ng, R.T., & Sander, J. (2000). LOF: Identifying density-based local outliers. *Proceedings of ACM SIGMOD 2000*, 93–104.

4. Rosch, E., & Mervis, C.B. (1975). Family resemblances: Studies in the internal structure of categories. *Cognitive Psychology*, 7(4), 573–605.

5. Hsu, G., Hannan, M.T., & Koçak, Ö. (2009). Multiple category memberships in markets. *American Sociological Review*, 74(1), 150–166.

6. Zuo, Y., & Serfling, R. (2000). General notions of statistical depth function. *The Annals of Statistics*, 28(2), 461–482.

7. Sklar, A. (1959). Fonctions de répartition à n dimensions et leurs marges. *Publications de l'Institut de Statistique de l'Université de Paris*, 8, 229–231.

8. Lovász, L. (1993). Random walks on graphs: A survey. *Combinatorics, Paul Erdős is Eighty*, 2, 1–46.

## License

MIT License — see [LICENSE](LICENSE) for details.
