# Platypus.jl

[![Build Status](https://github.com/simoneSantoni/Platypus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/simoneSantoni/Platypus.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-development-blue.svg)](https://simoneSantoni.github.io/Platypus.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-default_1.13.2-purple.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

<p align="center">
  <img src="docs/src/assets/logo.png" alt="Platypus.jl icon" width="160">
</p>

Platypus.jl computes **descriptive atypicality relative to a chosen collection**.
It scores two-feature projections with neighborhood, density, depth, or copula
methods and summarizes the scores for each entity. A separate graph estimator
measures connectivity using all selected features together.

The pairwise design is inspired by research on atypical combinations. Scores
depend on the selected features, reference collection, and method; they are not
anomaly probabilities, significance tests, or validated measures of human
judgments of novelty.

## Status and installation

The working tree targets **0.2.0**. Core scoring, copula models, graph scoring,
descriptive outcome regression, diagnostics, and plotting recipes are implemented.
This does not imply that a 0.2.0 release tag has been published.
See [CHANGELOG.md](CHANGELOG.md) for changes from 0.1.

The final 84-dataset numerical audit passed **14,637 checks with zero violations
and zero unexpected comparison errors**; the baseline also had zero unexpected
errors. The package suite passed 3,701/3,701 assertions, the independent
mathematical comparator checks 976/976, and the sociology checks 33/33.
Real-coordinate KDE now exposes a separate normalization ordering; the audit
checks its weighted mid-CDF contract explicitly. These results validate the
specified scenarios, not every possible input. See the
[numerical repair report](docs/research/numerical_repair.md) for evidence,
contract changes, preserved failures, and numerical limits.

The default development runtime is **Julia 1.13.2**, recorded in `.julia-version`
and pinned in CI and documentation builds. Julia 1.9 remains the minimum
supported version. From a checkout, configure Juliaup once:

```sh
juliaup add 1.13.2
juliaup override set "$(cat .julia-version)"
julia --version
```

The directory override makes `julia` select 1.13.2 in this checkout and its
subdirectories. On another checkout, repeat the setup, or use `julia +1.13.2`
explicitly. `.julia-version` records the choice; Juliaup needs the override command.

Install the development revision in a Julia project environment:

```julia
using Pkg
Pkg.add(url="https://github.com/simoneSantoni/Platypus.jl")
```

For research, retain the project's `Project.toml` and `Manifest.toml`, including
the resolved package commit. An unpinned repository installation follows the
revision available when installed. Record the dataset version and configuration.

## Reproducible example

Here 39 entities follow a diagonal relationship across three features. Entity 40
violates two of the three pairwise relationships.

```julia
using Platypus

x = collect(range(0.1, 0.9; length=39))
M = vcat(hcat(x, x, x), [0.95 0.05 0.95])
es = EntitySet(M, [:f1, :f2, :f3])
pairs = build_pairs(es; transform=:rank)
result = fit_typicality(pairs; method=:lof,
    provenance=(dataset="synthetic-diagonal-v1",))

@assert argmax(atypicality(result)) == 40
summarytable(result)
scoretable(result)             # normalized scores
scoretable(result; raw=true)   # method-specific scores
scoretable(result; normalization=true) # values used to order percentile groups
analysis_metadata(result)     # save alongside exported tables
```

This constructed example demonstrates the workflow, not general detection
performance. The [executed tutorial](docs/src/getting_started.md) explains the
output and data-loading choices.

## Score interpretation

| Output | Meaning | Direction |
| --- | --- | --- |
| `pair_scores` | Reference-weighted mid-percentiles of `normalization_scores` | Higher = more typical |
| `raw_scores` | Method-specific values, including physical log densities | Higher = more typical |
| `normalization_scores` | Raw values except real-coordinate KDE's ordinal ordering labels | Only order and ties matter |
| `conventionality` | Median over feature pairs | Higher = more typical |
| `novelty` | Lower quantile over feature pairs, default 0.1 | **Lower = more unusual** |
| `dispersion` | IQR over feature pairs | Greater variation across pairs |
| `atypicality` | `1 - conventionality` | Higher = more atypical |

`novelty` denotes lower-tail typicality, not historical newness or the minimum
pair score. Dispersion is not sampling uncertainty. Equal normalization values
share their reference-mass midpoint; if all are equal, every percentile is 0.5.

For Gaussian KDE on `:none` or `:zscore` coordinates, physical `raw_scores`
remain log densities. Percentiles instead use ordinal labels of dimensionless
kernel masses in `normalization_scores`. These labels have no density units or
comparable magnitudes across pairs/fits. Two rounded raw log densities can tie
while the ordering labels distinguish them, producing different percentiles.
Other methods, including rank-space KDE, use raw values for normalization.
`analysis_metadata(result).normalization_basis` records this choice.

Summaries always use normalized `pair_scores`; `normalize=false` is rejected.
No blanket tolerance merges nearby values. Underflow-affected KDE rows are
recomputed with a shift of the unweighted kernel exponents, retaining original
reference-weight ratios exactly, and compared using adaptive directed bounds.
These comparisons certify the represented expressions, not the exact original
Gaussian mixture; shifted terms can still round or underflow.

## Methods

| Entry point | Raw quantity | Nonuniform reference weights |
| --- | --- | --- |
| `fit_typicality(pairs; method=:lof)` | Reciprocal tie-inclusive LOF | Rejected |
| `method=:knn` | Negative kth-neighbor distance | Rejected |
| `method=:kde` | Log kernel density; reflected bounded kernels on ranks | Supported |
| `method=:depth` | Global Mahalanobis centrality | Supported |
| `fit_copula(pairs; family=:gaussian)` | Log Gaussian copula density | Supported |
| `fit_copula(pairs; family=:beta)` | Log empirical beta copula density | Uniform positive only; untied margins |
| `fit_graph_typicality(es)` | Stationary mass on a full-feature graph | Rejected |
| `fit_outcome(result, y)` | Descriptive least-squares association | Not a reference-weighting interface |

Supported weighted estimators apply `ReferenceWeights(w)` to transformations,
reference fitting, and score normalization. Zero-weight rows remain queries.
Labels are descriptive and never select peers automatically. Selecting a
collection still requires a substantive decision.

LOF/kNN and graph scoring use filtered exact squared-distance decisions: a
floating interval filter certifies clear comparisons, and rational arithmetic
resolves uncertain ones. Both requested backends, `:kdtree` and `:brute`,
currently use a conservative scan; no KD-tree pruning accelerator is active.
Candidate work is quadratic. There is no dense distance matrix, but ties can
make the graph dense. Metadata distinguishes requested and effective backends.

Depth uses exact rational two-dimensional moments, algebraic rank, and strict
affine support. Omit the retired `depth_rtol` argument; numeric values now raise
an error. Original represented reference weights and exact rank provenance are
retained through the relevant calculations. See the
[method guide](docs/src/guide/estimation.md) for geometry and precision limits.

`parallel=true` distributes typicality feature pairs across Julia threads.
Exact arithmetic adds allocation and integer-bit costs; KDE and empirical beta
copulas also require all-pairs work. The three retained Float64 score matrices
alone require approximately `24*N*P` bytes, where `P = p*(p-1)/2`. Benchmark the
complete workload.

## Data and visualization

Select measurement features explicitly to avoid including an outcome, timestamp,
or unrelated numeric covariate:

```julia
using DataFrames
df = DataFrame(id=1:3, genre=["jazz", "rock", "pop"],
               energy=[0.2, 0.9, 0.7], acousticness=[0.9, 0.1, 0.4],
               popularity=[30, 80, 50])
es = load_entities(df; id_col=:id, label_col=:genre,
                   feature_cols=[:energy, :acousticness])
```

IDs and feature names must be unique. Missing and nonfinite features are rejected.
`EntitySet.matrix` is authoritative; entity accessors are snapshots. Rebuild
pairs after editing features. Validation requires every unordered feature pair
exactly once. For pairwise fits, save metadata's ordered `entity_ids` with both
normalized `reference_weights` and original `source_reference_weights`; this preserves
their association when exported score tables are sorted or joined.

All nine `plot_*` functions return `TypicalityPlot` recipes with detached data.
Install Plots separately to render them:

```julia
using Plots
plot(plot_novelty_conventionality(result))
plot(plot_typicality_profile(result, :entity_id => 40))
```

Platypus depends on RecipesBase and does not choose a plotting backend.

## Documentation and validation

- [Getting started](docs/src/getting_started.md)
- [Methods and assumptions](docs/src/guide/estimation.md)
- [Research workflow](docs/src/guide/research_workflow.md)
- [Algorithm research and decision records](docs/research/README.md)
- [References](docs/src/references.md)
- [Development documentation](https://simoneSantoni.github.io/Platypus.jl/dev/)

From a checkout:

```sh
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Numerical tests and synthetic validation evaluate implementation properties and
specified scenarios. Human-perception validation is outside this development
phase. No causal or historical interpretation follows from a descriptive fit.

MIT License; see [LICENSE](LICENSE).
