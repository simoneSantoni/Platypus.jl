# Building Feature Pairs

[`build_pairs`](@ref) projects ``N`` entities with ``p`` features onto
``P=p(p-1)/2`` pairs. At least two features are required.

```@example pairs
using Platypus
es = EntitySet([0.1 0.2 0.9; 0.4 0.7 0.6; 0.8 0.5 0.2], [:a, :b, :c])
pairs = build_pairs(es)
@assert n_pairs(pairs) == 3
fp = pairs[1]
(fp.name_i, fp.name_j, fp.u, fp.v)
```

## Feature transformations

| `transform` | Operation | Interpretation |
| --- | --- | --- |
| `:rank` (default) | Empirical mid-rank divided by number of entities | Removes marginal units and spaces values by rank |
| `:zscore` | Center and divide by sample standard deviation | Preserves marginal shape; constant features are centered |
| `:none` | Keep raw values | Distances depend on supplied units |

For uniform reference weights, the rank transform is

```math
u_{ej} = \frac{\operatorname{tiedrank}(x_{\bullet j})_e - 1/2}{N}.
```

Untied continuous observations form a uniform grid. Ties retain probability
masses; constant features remain constant. Ranking is not a null model and
does not guarantee dependence-only anomaly detection. Local geometry, boundary
effects, and the selected score still matter.

```@example pairs
raw_pairs = build_pairs(es; transform=:none)
standardized_pairs = build_pairs(es; transform=:zscore)
@assert raw_pairs[1].u == es.matrix[:, 1]
(raw_pairs.transform, standardized_pairs.transform)
```

When a supported estimator receives nonuniform reference weights, it rebuilds
the transforms from the stored raw snapshot. Weighted ranks use the midpoint
of each tied value's reference mass. Zero-weight rows are queries and can map
to endpoints zero or one. Weighted standardization also uses reference mass.

## Storage and provenance

Pair coordinates share one vector per feature, avoiding a repeated copy for
every pair. The `PairSet` stores the transform and a detached source matrix for
weighted fitting. Editing the source `EntitySet` does not change existing pairs;
rebuild them after data changes. Avoid direct coordinate mutation: shared
vectors affect multiple pairs, and weighted fitting uses the source snapshot.

A manually constructed `PairSet` without a raw snapshot cannot reconstruct
weighted transforms and rejects nonuniform weights. Prefer `build_pairs`.

## Limits of projections

All pairs share features and observations; they are not independent replicates.
Duplicating a feature changes the distribution of pair scores and therefore the
summaries. Correlation diagnostics do not supply justified aggregation weights.

Two-dimensional projections can miss higher-order interactions. The separate
full-feature graph estimator changes the representation and estimand; it does
not guarantee recovery of every such interaction. See
[Estimation](estimation.md) and [Exploratory Analysis](eda.md).
