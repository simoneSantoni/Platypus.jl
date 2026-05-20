# Getting Started

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/simoneSantoni/Platypus.jl")
```

## A worked example

The example below builds a synthetic collection in which most entities
lie on a dependence ridge between two correlated features, plus a handful
of entities that violate that dependence. Those violators should come out
as the most atypical.

```julia
using Platypus

# 200 entities, 5 features. Features 1 and 2 are correlated; the last
# 10 entities break the correlation (high feature 1, low feature 2).
n = 200
M = rand(n, 5)
M[:, 2] .= 0.7 .* M[:, 1] .+ 0.3 .* M[:, 2]      # feature 2 tracks feature 1
M[end-9:end, 1] .= 0.95                          # the atypical entities
M[end-9:end, 2] .= 0.05

es = EntitySet(M, [:f1, :f2, :f3, :f4, :f5])

# Project onto all 5*4/2 = 10 feature pairs (copula-rank transform).
pairs = build_pairs(es)

# Score local typicality on each pair and aggregate per entity.
result = fit_typicality(pairs; method=:lof)
println(result)

# Per-entity summary table.
summary = summarytable(result)

# The headline quantity: higher = more atypical.
atyp = atypicality(result)
println("most atypical entity: ", argmax(atyp))   # one of the last 10
```

## What you get back

[`fit_typicality`](@ref) returns a [`TypicalityResult`](@ref) with:

- `pair_scores` — the full `n_entities × n_pairs` matrix of pairwise
  typicality scores.
- `conventionality`, `novelty`, `dispersion`, `atypicality` — the
  per-entity summary vectors.

Access them through the accessors and tables:

```julia
conventionality(result)     # median pairwise typicality, per entity
novelty(result)             # lower-tail (10th percentile) per entity
dispersion(result)          # inter-quartile range per entity

summarytable(result)        # one row per entity
scoretable(result)          # long-form: one row per entity × feature-pair
```

## Loading real data

```julia
using DataFrames

df = DataFrame(
    id           = 1:3,
    genre        = ["jazz", "rock", "pop"],
    danceability = [0.31, 0.74, 0.66],
    acousticness = [0.88, 0.12, 0.40],
    energy       = [0.22, 0.91, 0.70],
)

# `genre` is carried as a label; it is never used to define the
# reference set. Every numeric column becomes a feature.
es = load_entities(df; id_col=:id, label_col=:genre)
```

CSV files load the same way: `load_entities("songs.csv"; id_col=:id)`.
