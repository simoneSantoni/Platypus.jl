# Building Feature Pairs

[`build_pairs`](@ref) turns an `EntitySet` into the object the estimators
consume: the set of all ``n(n-1)/2`` two-dimensional projections of the
collection.

```julia
using Platypus

es = EntitySet(rand(300, 6), [:f1, :f2, :f3, :f4, :f5, :f6])
pairs = build_pairs(es)

n_pairs(pairs)     # 15  =  6*5/2
```

## The copula-rank transform

Before pairs are formed, each feature is transformed. The default
`transform=:rank` replaces a feature by its rank rescaled to ``(0, 1)``:

```math
u_i = \frac{\operatorname{rank}(x_i) - \tfrac{1}{2}}{n}.
```

This is the empirical copula transform. It makes every feature marginally
uniform, so a feature pair's cloud reflects only the **dependence
structure** of the two features — and an entity scores as atypical
because of an unusual *combination*, not because one feature is extreme
on its own. This is what makes the pairwise measure a measure of
*combination* atypicality.

Two other transforms are available:

```julia
build_pairs(es; transform=:rank)     # copula ranks (default, recommended)
build_pairs(es; transform=:zscore)   # standardised; keeps marginal shape
build_pairs(es; transform=:none)     # raw feature values
```

Use `:zscore` or `:none` when extreme *marginal* values should still
count as atypical; use `:rank` (the default) when only the combination
should.

## What a pair looks like

Each `FeaturePair` stores the two transformed coordinates, one point per
entity:

```julia
fp = pairs[1]
fp.name_i, fp.name_j      # the two feature names
fp.u, fp.v                # the coordinates of the cloud
```

A collection with fewer than two features has no pairs to form, and
`build_pairs` throws an `ArgumentError`. The redundancy among pairs —
pairs built from highly correlated features carry little independent
information — can be inspected with [`pair_redundancy`](@ref); see
[Exploratory Analysis](eda.md).

API reference: [Types](../api/types.md).
