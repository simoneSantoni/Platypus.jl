# Entities and Data

An `Entity` records one member of a collection: an `id`, a vector of
continuous `features`, and an optional category `label`.

```julia
using Platypus

song = Entity(id=1, features=[0.74, 0.91, 0.12], label="rock")
features(song)     # [0.74, 0.91, 0.12]
label(song)        # "rock"
```

The `label` is a soft tag — a genre, a sector, a field. It is carried
through so results can be inspected by group, but it is **never** used to
define the reference set: Platypus.jl treats the reference set as
endogenous (see [Home](../index.md)).

Collect entities into an `EntitySet`, which shares one feature schema
across all members and caches an `n_entities × n_features` value matrix.

```julia
es = EntitySet([
        Entity(id=1, features=[0.74, 0.91, 0.12]),
        Entity(id=2, features=[0.31, 0.22, 0.88]),
        Entity(id=3, features=[0.66, 0.70, 0.40]),
    ], [:danceability, :energy, :acousticness])

n_entities(es)     # 3
n_features(es)     # 3
feature_names(es)  # [:danceability, :energy, :acousticness]
```

An `EntitySet` can also be built directly from a numeric matrix, with
optional ids and labels:

```julia
es = EntitySet(rand(100, 6), [:f1, :f2, :f3, :f4, :f5, :f6])
```

Most often the data starts in a `DataFrame` or CSV file.
[`load_entities`](@ref) ingests it: by default every numeric column
becomes a feature, while the id and label columns are set aside.

```julia
using DataFrames

df = DataFrame(id=1:3, genre=["jazz", "rock", "pop"],
               danceability=[0.31, 0.74, 0.66],
               acousticness=[0.88, 0.12, 0.40])

es = load_entities(df; id_col=:id, label_col=:genre)
```

Missing feature values are rejected — impute or drop incomplete rows
before loading, because the density estimators in
[`fit_typicality`](@ref) have no missing-value semantics.

API reference: [Types](../api/types.md).
