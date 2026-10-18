# Exploratory Analysis

Inspect feature distributions, feature selection, dependence, and nearby
observations before fitting. These helpers return matrices or tables; they do
not choose a reference population or validate the substantive construct.

```@example eda
using Platypus
es = EntitySet([0.1 0.2 0.9; 0.4 0.7 0.6; 0.8 0.5 0.2; 0.6 0.9 0.3],
    [:a, :b, :c])
feature_summary(es)
```

```@example eda
rank_matrix(es)
pair_correlations(es)
compute_entity_statistics(es)
```

`rank_matrix` returns the unweighted mid-rank matrix. `pair_correlations`
reports Spearman correlations. `compute_entity_statistics` reports mean rank,
rank dispersion, and marginal extremeness; these are distinct from fitted
typicality scores.

## Dependence diagnostics

```@example eda
pairs = build_pairs(es)
diagnostic = pair_redundancy(pairs)
@assert propertynames(diagnostic) == [:feature_a, :feature_b, :abs_correlation]
diagnostic
```

Despite the historical name, `pair_redundancy` measures only the absolute
correlation of the two coordinates within each pair. It does not measure
dependence between pair-score columns, and it does not return aggregation
weights. Undefined correlations, including constant features, remain `NaN`.

A highly correlated pair may be especially useful for detecting a violation of
its relationship. Automatically assigning it low weight is not justified.
The pair scores share features and observations, even if within-pair
correlations are small. Their IQR is descriptive variation, not inferential
uncertainty. Evaluate sensitivity to the selected feature schema.

## Neighbors and pair inventory

```@example eda
nearest_entities(es, 1; k=2)
pair_table(pairs)
```

`nearest_entities` takes a row index, not an entity ID. It computes Euclidean
distance across all features in unweighted rank space; it is an exploratory
helper, not the pair-specific LOF neighborhood. It returns up to `k` rows and
does not have the tie-inclusive scoring contract of LOF.

Plot recipes are described in [Visualization](../api/visualization.md).
