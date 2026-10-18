# Getting Started

## Installation and reproducibility

Use **Julia 1.13.2**, the default development and documentation runtime.
Julia 1.9 remains the minimum supported version. With Juliaup, install the
default runtime using `juliaup add 1.13.2`. In a repository checkout, run
`juliaup override set 1.13.2` once so `julia` selects it in that directory and
its subdirectories. Alternatively, invoke `julia +1.13.2` explicitly.

Install the development revision in a Julia project environment:

```julia
using Pkg
Pkg.add(url="https://github.com/simoneSantoni/Platypus.jl")
```

This URL installation is not pinned to a release. Retain `Project.toml` and
`Manifest.toml` so the resolved revision and dependencies can be recovered.
The source tree targets 0.2.0; a release tag is a separate publication step.

When upgrading an existing checkout from Julia 1.10, preserve its generated
`Manifest.toml` as `Manifest-v1.10.toml` before resolving with Julia 1.13.2.
Do the same for `docs/Manifest.toml` if present. Run
`julia --project -e 'using Pkg; Pkg.instantiate()'` to resolve dependencies for
the new runtime. The new manifest can be named `Manifest-v1.13.toml` to keep
the environments separate. Old manifests can pin incompatible standard-library
versions even when the package itself supports the new Julia version.

## A worked example

This deterministic example executes when the documentation builds. Thirty-nine
entities follow a diagonal relationship in three features. The last entity
violates two of the three pairwise relationships, a pattern that the median over
pairs can retain.

```@example getting_started
using Platypus

x = collect(range(0.1, 0.9; length=39))
M = vcat(hcat(x, x, x), [0.95 0.05 0.95])
es = EntitySet(M, [:f1, :f2, :f3])
pairs = build_pairs(es; transform=:rank)
result = fit_typicality(pairs; method=:lof,
    provenance=(dataset="synthetic-diagonal-v1", selection="all 40 rows"))

@assert n_pairs(pairs) == 3
@assert argmax(atypicality(result)) == 40
summarytable(result)[[1, 20, 40], :]
```

This constructed example does not establish general detection performance.
Unusualness confined to one pair can disappear in a median over many pairs.
A compact unusual cluster may also look locally ordinary under LOF.

## Inspect the measurements

```@example getting_started
# Higher headline atypicality means less typical overall.
@assert atypicality(result)[40] > atypicality(result)[20]

# Lower novelty means more unusual lower-tail typicality.
tail_score = novelty(result)[40]
pair_index = argmin(result.pair_scores[40, :])
(result.entity_ids[40], tail_score, result.pair_names[pair_index])
```

`novelty` is a lower quantile, not a minimum or historical newness. The IQR
`dispersion` describes variation among pair scores, not uncertainty.

```@example getting_started
normalized_table = scoretable(result)
raw_table = scoretable(result; raw=true)
fit_metadata = analysis_metadata(result)
@assert fit_metadata.transform == :rank
@assert fit_metadata.normalization == :weighted_midcdf
(fit_metadata.method, fit_metadata.k, fit_metadata.provenance)
```

Save tables and metadata with the dataset identifier, analysis code, and Julia
environment. Normalized scores are reference-weighted mid-percentiles; raw LOF
scores here are reciprocals of LOF. Raw scales differ between methods.

## Load selected features

Numeric columns are detected automatically if `feature_cols` is omitted.
Explicit selection avoids treating an outcome or timestamp as a feature.

```@example getting_started
using DataFrames

df = DataFrame(
    id=1:3,
    genre=["jazz", "rock", "pop"],
    energy=[0.22, 0.91, 0.70],
    acousticness=[0.88, 0.12, 0.40],
    popularity=[30, 80, 50],
)
songs = load_entities(df; id_col=:id, label_col=:genre,
    feature_cols=[:energy, :acousticness])
@assert feature_names(songs) == [:energy, :acousticness]
feature_summary(songs)
```

CSV files use the same keywords:
`load_entities("songs.csv"; id_col=:id, feature_cols=[:energy, :acousticness])`.
Handle missing and nonfinite features before loading; they are rejected.
IDs must be unique integers, and feature names must be unique.

## Create a visualization

Plot helpers create detached recipe objects without loading a graphics backend:

```@example getting_started
recipe = plot_typicality_profile(result, :entity_id => 40)
@assert recipe isa TypicalityPlot
nothing
```

To render, install Plots in your environment and run:

```julia
using Plots
plot(recipe)
plot(plot_novelty_conventionality(result))
```

The second plot labels its vertical axis as lower-tail typicality: lower values
are more unusual. Read the [Estimation guide](guide/estimation.md) before changing
methods or adding reference weights.
