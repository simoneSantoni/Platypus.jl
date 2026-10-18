# Entities and Data

An [`Entity`](@ref) contains an integer ID, finite continuous feature values,
and an optional descriptive label. An [`EntitySet`](@ref) puts entities under
a common feature schema.

```@example entities
using Platypus, DataFrames
song = Entity(id=1, features=[0.74, 0.91, 0.12], label="rock")
features(song), label(song)
```

## Select measurement features

```@example entities
df = DataFrame(
    id=1:3, genre=["jazz", "rock", "pop"],
    energy=[0.22, 0.91, 0.70], acousticness=[0.88, 0.12, 0.40],
    popularity=[30, 80, 50], year=[2020, 2021, 2021],
)
es = load_entities(df; id_col=:id, label_col=:genre,
    feature_cols=[:energy, :acousticness])
@assert n_entities(es) == 3
@assert n_features(es) == 2
feature_names(es)
```

Without `feature_cols`, every numeric column except the selected ID and label
columns becomes a feature. Explicit selection prevents outcomes, timestamps, and
administrative covariates from entering the measurement unintentionally.

`load_entities("songs.csv"; feature_cols=[:energy, :acousticness])` applies the
same interface to CSV files. IDs default to row numbers if the configured ID
column is absent or `id_col=nothing`. IDs must be unique integers, and feature
names must be unique. Missing, infinite, and NaN feature values are rejected;
choose and report an imputation or exclusion procedure before loading.

Labels are carried into outputs for inspection. They do not define reference
groups, and calling a label “soft” does not provide graded memberships.
Use explicit [`ReferenceWeights`](@ref) with a supported estimator for that.

## Matrix construction and snapshots

```@example entities
matrix_set = EntitySet([0.2 0.8; 0.5 0.4; 0.9 0.1], [:energy, :acousticness];
    ids=[10, 20, 30], labels=["a", "b", "c"])
saved = build_pairs(matrix_set)
matrix_set.matrix[1, 1] = 0.3
@assert saved.source_matrix[1, 1] == 0.2
features(matrix_set[1])
```

`EntitySet.matrix` is the authoritative mutable matrix. Entity accessors and
feature-name accessors return snapshots; changing them does not edit the set.
Use the matrix for intentional data edits, then rebuild pairs and refit.
Construction and fitting validate the schema, IDs, and finite values.

`build_pairs` saves a detached data snapshot; results and plot recipes also
detach their inputs. Existing estimates do not update when the original data
changes. Direct mutation of public result fields remains possible in Julia;
treat fitted objects as records of an analysis and use a new fit for changes.

See [Types](../api/types.md) for exact constructors.
