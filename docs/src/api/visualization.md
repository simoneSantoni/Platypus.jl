# Visualization API

For software attribution and BibTeX, see [Citing Platypus.jl](citation.md).

All nine plotting helpers are implemented and return [`TypicalityPlot`](@ref)
recipe objects containing detached data. Creating a recipe requires no graphics
backend. Rendering uses a separately installed compatible plotting package:

```julia
using Platypus, Plots
plot(plot_novelty_conventionality(result))
```

Install Plots in your project environment before running this rendering example.
Platypus depends only on RecipesBase for plotting and does not select a backend.
Use ordinary Plots keyword arguments on the outer `plot` call to style a figure.

## Accepted inputs

| Helper | Input | Values shown |
|---|---|---|
| `plot_feature_distributions` | `EntitySet` | Raw feature histograms |
| `plot_pair_cloud` | `PairSet`, pair index or ordered name tuple | Stored transformed pair coordinates |
| `plot_correlation_heatmap` | `EntitySet` | Spearman correlations |
| `plot_pair_redundancy` | `PairSet` | Absolute correlations of stored pair coordinates |
| `plot_typicality_profile` | `TypicalityResult`, row index or `:entity_id => id` | One entity's normalized pair scores |
| `plot_atypicality_distribution` | `TypicalityResult` or `GraphTypicalityResult` | Headline atypicality |
| `plot_novelty_conventionality` | `TypicalityResult` | Median and lower-quantile pairwise typicality |
| `plot_pair_heatmap` | `TypicalityResult` | Normalized entity × pair scores |
| `plot_outcome_fit` | `OutcomeResult` | Observed and fitted outcomes |

Only `plot_atypicality_distribution` accepts graph results. Pairwise score plots
use normalized scores, including the explicit KDE normalization ordering; they
do not plot physical log densities or the ordinal labels directly. To plot raw
or normalization values, export them with [`scoretable`](@ref) and pass the
table to your plotting package.

`TypicalityPlot.series` holds detached series descriptions and `layout` is a
`(rows, columns)` tuple. Editing source data or refitting does not update an
existing recipe; construct a new one for the new result.

```@docs
TypicalityPlot
```

## Features and pairs

`plot_pair_cloud` accepts a one-based pair index or an ordered tuple of feature
names matching the pair schema. Correlation plots retain undefined correlations
for constant features. The redundancy plot shows within-pair absolute
correlation; it does not apply aggregation weights.

```@docs
plot_feature_distributions
plot_pair_cloud
plot_correlation_heatmap
plot_pair_redundancy
```

## Fitted scores

`plot_typicality_profile(result, 3)` selects row 3.
`plot_typicality_profile(result, :entity_id => 30)` selects entity ID 30.
The novelty–conventionality plot labels the vertical axis as lower-tail
typicality, where smaller means more unusual.

The following example creates and inspects a recipe without installing a
rendering backend:

```@example api_visualization
using Platypus

es = EntitySet([0 0; 1 2; 2 1; 4 4], [:x, :y]; ids=[10, 20, 30, 40])
result = fit_typicality(build_pairs(es); method=:lof, k=2)
recipe = plot_typicality_profile(result, :entity_id => 30)
@assert recipe isa TypicalityPlot # hide
@assert recipe.series[1].args[1] == result.pair_scores[3, :] # hide
(layout=recipe.layout, title=recipe.series[1].title)
```

```@docs
plot_typicality_profile
plot_atypicality_distribution
plot_novelty_conventionality
plot_pair_heatmap
```

## Outcome association

`plot_outcome_fit` accepts the `OutcomeResult` returned by `fit_outcome`.
Its identity line compares observed and fitted values; it is not a confidence
band or an out-of-sample validation result.

```@docs
plot_outcome_fit
```
