# Estimation

The implemented estimator is [`fit_typicality`](@ref): it scores every
entity's local typicality on each feature pair and aggregates the
``n(n-1)/2`` scores per entity. More elaborate estimators
(`fit_graph_typicality`, `fit_outcome`, `fit_copula`) are reserved in the
public API but currently raise `not implemented`.

```julia
using Platypus

es = EntitySet(rand(400, 7), [:f1, :f2, :f3, :f4, :f5, :f6, :f7])
pairs = build_pairs(es)

result = fit_typicality(pairs; method=:lof)

result.pair_scores          # n_entities × n_pairs typicality matrix
result.conventionality      # median pairwise typicality, per entity
result.novelty              # lower-tail (10th percentile) per entity
result.dispersion           # inter-quartile range per entity
result.atypicality          # 1 .- conventionality
```

## Choosing a method

```julia
fit_typicality(pairs; method=:knn)     # k-th nearest-neighbour density
fit_typicality(pairs; method=:lof)     # Local Outlier Factor (default)
fit_typicality(pairs; method=:kde)     # Gaussian kernel density
fit_typicality(pairs; method=:depth)   # Mahalanobis data depth
```

`:lof` is the default: it scores an entity's density *relative to its
neighbours'* density, so it adapts to a collection whose density varies
from region to region — the heterogeneous case Platypus.jl is built for.

The neighbourhood size for `:knn` and `:lof` is set with `k`; the lower
quantile used for `novelty` is set with `novelty_quantile`:

```julia
fit_typicality(pairs; method=:lof, k=20, novelty_quantile=0.05)
```

## Graded reference membership

When membership in the reference set is graded rather than binary, pass
[`ReferenceWeights`](@ref). They are honoured by the `:kde` and `:depth`
methods:

```julia
w = ReferenceWeights(recency_weights)   # one non-negative weight per entity
fit_typicality(pairs; method=:kde, weights=w)
```

## Accessors and tables

```julia
atypicality(result)         # = result.atypicality
conventionality(result)
novelty(result)
dispersion(result)

summarytable(result)        # one row per entity
scoretable(result)          # long-form: one row per entity × feature-pair
```

API reference: [Estimation API](../api/estimation.md).
