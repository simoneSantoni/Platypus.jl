# Types and Data API

For software attribution and BibTeX, see [Citing Platypus.jl](citation.md).

```@docs
Platypus
```

The data model separates records, collections, and projections. An
`EntitySet` contains the analyst's selected collection and feature schema.
A `PairSet` contains all two-feature projections and a detached source snapshot.
Reference weights are graded contributions to supported estimators, not labels.

`EntitySet.matrix` is the authoritative mutable matrix. Entity accessors return
snapshots. Rebuild pairs after editing the matrix, and refit to produce a new
result. Pair coordinate vectors are shared within a `PairSet`; direct field
mutation is an advanced operation and is not an update API.

## Data types

### Construction and ownership

| Constructor | Required input | Optional input |
|---|---|---|
| `Entity(; id, features, label=nothing)` | Integer ID and a real-valued feature vector | String label |
| `EntitySet(matrix, feature_names; ids=nothing, labels=nothing)` | Nonempty real matrix and one name per column | IDs default to `1:N`; labels default to `nothing` |
| `EntitySet(entities, feature_names)` | `Vector{Entity}` with equal feature lengths and a `Vector{Symbol}` schema | None |
| `ReferenceWeights(weights)` | Finite nonnegative real vector with positive total mass | None |
| `ReferenceWeights(n)` | Positive entity count | Creates equal weights |

Feature values and reference weights are stored as Float64. IDs and feature
names must be unique; features must remain finite after conversion. A collection
may have one feature, but [`build_pairs`](@ref) requires at least two. Constructors
copy feature data; subsequent changes to the original matrix do not change the
collection. Public storage is mutable even though the enclosing structs are
immutable.

`es[i]` selects a **row**, not an entity ID, and returns a detached `Entity`.
Iteration, `es[:]`, integer-vector indexing, and Boolean-mask indexing return
entity snapshots. A multirow selection returns a `Vector{Entity}`, not another
`EntitySet`; use `EntitySet(selected, feature_names(es))` to create a collection.
`ps[i]` selects a `FeaturePair` in pair order and exposes its shared coordinates.

```@example api_types
using Platypus

input = [1.0 4.0; 2.0 2.0; 3.0 5.0]
es = EntitySet(input, [:x, :y]; ids=[10, 20, 30])
selected = EntitySet(es[[true, false, true]], feature_names(es))
@assert selected.entity_ids == [10, 30] # hide
@assert es[2].id == 20 # hide
input[1, 1] = 99.0
@assert es.matrix[1, 1] == 1.0 # hide
(original=es, selected=selected)
```

`ReferenceWeights` bind to input **row order**; they contain no entity IDs.
At fitting time their length must match the collection. Pairwise estimators
require at least two positive weights. Zero-weight queries are supported only
by KDE, depth, and Gaussian copulas; LOF, kNN, beta copulas, and graph scoring
require uniform positive weights. Supported fits recompute marginal transforms
using the original represented weights. A manually constructed `PairSet` without
source provenance cannot support nonuniform weighting.

```@docs
Entity
EntitySet
FeaturePair
PairSet
ReferenceWeights
TYPICALITY_METHODS
```

## Loading and projection

Explicitly choose `feature_cols` when the input also contains outcomes,
timestamps, or administrative variables. IDs and feature names must be unique;
missing and nonfinite feature values are rejected.

`load_entities` accepts a CSV path or a `DataFrame`. Both overloads use
`feature_cols=nothing`, `id_col=:id`, and `label_col=nothing`. With no explicit
feature list, numeric columns other than ID/label columns become features;
string-encoded features need explicit selection. Absent or disabled ID columns
produce sequential IDs. Missing labels become `nothing`; missing features raise.

`build_pairs(es; transform=:rank)` creates all `F*(F-1)/2` unordered pairs for
`F` features, in increasing column-index order. Other transforms are `:zscore`
and `:none`. The returned `source_matrix` is detached from `es.matrix`; rebuild
pairs after editing the collection. Estimators require each unordered pair
exactly once and reject duplicated, missing, or source-inconsistent projections.
Reordering a complete, consistent pair collection is allowed.

```@example api_types
using DataFrames

table = DataFrame(id=[10, 20, 30], x=[1, 2, 3], y=[4, 2, 5],
                  outcome=[7, 9, 8])
loaded = load_entities(table; feature_cols=[:x, :y])
pairs = build_pairs(loaded; transform=:rank)
@assert feature_names(loaded) == [:x, :y] # hide
@assert n_pairs(pairs) == 1 # hide
pair_table(pairs)
```

```@docs
load_entities
build_pairs
```

## Exploratory analysis

These return descriptive tables or matrices. Within-pair correlations are
diagnostics, not weights for aggregating pair scores. The exploratory nearest
neighbors use all features in unweighted rank space.

`nearest_entities(es, idx; k=5)` takes a row index and returns up to `N-1`
neighbors. This exploratory helper uses displayed Float64 ranks and returns
at most `k` rows; it is not the tie-inclusive neighborhood oracle used by LOF,
kNN, or graph estimators. Constant-feature correlations may be `NaN`, preserving
their undefined status.

```@docs
feature_summary
rank_matrix
pair_correlations
pair_redundancy
nearest_entities
compute_entity_statistics
pair_table
```

## Accessors

```@docs
features
label
feature_names
n_entities
n_features
n_pairs
```

## Parsing helper

```@docs
Platypus.parse_feature
```
