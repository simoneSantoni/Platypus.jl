# Exploratory Analysis

Before fitting a typicality model it is worth looking at the marginal
shape of each feature and — more importantly — at how strongly the
features are correlated. The exploratory helpers operate on an
`EntitySet` (or a `PairSet`) and return `DataFrame`s.

```julia
using Platypus

es = EntitySet(rand(300, 6), [:f1, :f2, :f3, :f4, :f5, :f6])

feature_summary(es)            # per-feature mean, std, quartiles, range
rank_matrix(es)                # the copula-rank transform of the matrix
pair_correlations(es)          # Spearman correlation of every feature pair
compute_entity_statistics(es)  # per-entity mean rank, dispersion, extremeness
```

## The redundancy check

The ``n(n-1)/2`` pairwise typicality scores are **not** an i.i.d. sample:
features are correlated, so pairs built from correlated features carry
overlapping information. This matters when you summarise the pairwise
distribution — a quantile over redundant pairs is noisier than its count
suggests.

[`pair_correlations`](@ref) and [`pair_redundancy`](@ref) expose this
directly:

```julia
pairs = build_pairs(es)
pr = pair_redundancy(pairs)    # abs_correlation and a suggested weight per pair
```

A pair whose two features are nearly collinear gets `weight` close to
`0`; a pair of near-independent features gets `weight` close to `1`. The
weights are a suggested down-weighting when aggregating — pairs that say
the same thing should not be counted twice.

## Locating an entity's neighbours

```julia
nearest_entities(es, 1; k=5)   # the 5 entities closest to entity 1
```

Distance is Euclidean in the copula-rank space, so it is comparable
across features regardless of their original units.

## Listing the pairs

```julia
pair_table(pairs)              # one row per feature pair
```

API reference: [Types](../api/types.md) (for `EntitySet` and `PairSet`),
[Estimation API](../api/estimation.md) (for the model that consumes a
`PairSet`).
