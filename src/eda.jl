"""
Exploratory data analysis for entity collections.

These helpers operate on an [`EntitySet`](@ref) (or a [`PairSet`](@ref))
and return summary scalars, matrices, or `DataFrame`s. They are meant to
be run *before* [`fit_typicality`](@ref) — to inspect the marginal shape
of each feature, to see how strongly features are correlated (and hence
how redundant the feature pairs are), and to locate an entity's nearest
neighbours.
"""

"""
    feature_summary(es::EntitySet) -> DataFrame

Per-feature descriptive statistics. Columns: `feature`, `mean`, `std`,
`min`, `q25`, `median`, `q75`, `max`.
"""
function feature_summary(es::EntitySet)
    rows = NamedTuple[]
    for c in 1:es.n_features
        x = es.matrix[:, c]
        push!(rows, (
            feature = es.feature_names[c],
            mean    = mean(x),
            std     = std(x),
            min     = minimum(x),
            q25     = quantile(x, 0.25),
            median  = median(x),
            q75     = quantile(x, 0.75),
            max     = maximum(x),
        ))
    end
    return DataFrame(rows)
end

"""
    rank_matrix(es::EntitySet) -> Matrix{Float64}

The empirical copula transform of the whole data matrix: an
`n_entities × n_features` matrix in which every column is the focal
feature replaced by its rank rescaled to `(0, 1)`. This is the same
transform [`build_pairs`](@ref) applies with `transform=:rank`.
"""
function rank_matrix(es::EntitySet)
    R = Matrix{Float64}(undef, es.n_entities, es.n_features)
    for c in 1:es.n_features
        R[:, c] = _rank_transform(@view es.matrix[:, c])
    end
    return R
end

"""
    pair_correlations(es::EntitySet) -> DataFrame

Spearman (rank) correlation between every pair of features. Returns a
long-form `DataFrame` with columns `feature_a`, `feature_b`,
`correlation`, one row per unordered pair.

Strong correlations are the warning sign discussed in the manual: highly
correlated features make their pair carry little independent information,
so the `n(n-1)/2` pairwise typicality scores are not an i.i.d. sample.
See [`pair_redundancy`](@ref).
"""
function pair_correlations(es::EntitySet)
    C = corspearman(es.matrix)
    rows = NamedTuple[]
    for i in 1:(es.n_features - 1), j in (i + 1):es.n_features
        push!(rows, (
            feature_a   = es.feature_names[i],
            feature_b   = es.feature_names[j],
            correlation = C[i, j],
        ))
    end
    return DataFrame(rows)
end

"""
    pair_redundancy(ps::PairSet) -> DataFrame

Redundancy diagnostic for the feature pairs of a [`PairSet`](@ref).

For each pair, `abs_correlation` is the absolute correlation between the
pair's two (transformed) coordinates, and `weight = 1 - abs_correlation`
is a suggested down-weight: a pair whose two features are nearly collinear
adds little to the aggregate and can be discounted when summarising the
pairwise scores. Columns: `feature_a`, `feature_b`, `abs_correlation`,
`weight`.
"""
function pair_redundancy(ps::PairSet)
    rows = NamedTuple[]
    for fp in ps.pairs
        r = abs(cor(fp.u, fp.v))
        isnan(r) && (r = 0.0)
        push!(rows, (
            feature_a       = fp.name_i,
            feature_b       = fp.name_j,
            abs_correlation = r,
            weight          = 1.0 - r,
        ))
    end
    return DataFrame(rows)
end

"""
    nearest_entities(es::EntitySet, idx::Integer; k::Integer=5) -> DataFrame

The `k` entities closest to entity `idx` (a row index), measured by
Euclidean distance in the copula-rank space of [`rank_matrix`](@ref).
Returns a `DataFrame` with columns `entity_id`, `label`, `distance`,
sorted by increasing distance.
"""
function nearest_entities(es::EntitySet, idx::Integer; k::Integer=5)
    1 <= idx <= es.n_entities ||
        throw(BoundsError("entity index $idx out of range 1:$(es.n_entities)"))
    R = rank_matrix(es)
    d = Vector{Float64}(undef, es.n_entities)
    for a in 1:es.n_entities
        d[a] = sqrt(sum(abs2, @view(R[a, :]) .- @view(R[idx, :])))
    end
    order = sortperm(d)
    order = filter(!=(idx), order)
    take = order[1:min(k, length(order))]
    return DataFrame(
        entity_id = es.entity_ids[take],
        label     = es.labels[take],
        distance  = d[take],
    )
end

"""
    compute_entity_statistics(es::EntitySet) -> DataFrame

Per-entity descriptive statistics computed without fitting a model.
Works on the copula-rank representation, so the statistics are
comparable across features. Columns: `entity_id`, `label`, `mean_rank`,
`rank_dispersion`, `extremeness`.

- `mean_rank` — average rank of the entity across features.
- `rank_dispersion` — standard deviation of those ranks.
- `extremeness` — the largest deviation `|rank - 0.5|`, i.e. how
  marginally extreme the entity is on its most unusual single feature.
"""
function compute_entity_statistics(es::EntitySet)
    R = rank_matrix(es)
    rows = NamedTuple[]
    for i in 1:es.n_entities
        r = @view R[i, :]
        push!(rows, (
            entity_id       = es.entity_ids[i],
            label           = es.labels[i],
            mean_rank       = mean(r),
            rank_dispersion = es.n_features > 1 ? std(r) : 0.0,
            extremeness     = maximum(abs.(r .- 0.5)),
        ))
    end
    return DataFrame(rows)
end

"""
    pair_table(ps::PairSet) -> DataFrame

The feature pairs of a [`PairSet`](@ref) as a `DataFrame`, one row per
pair. Columns: `pair_index`, `feature_a`, `feature_b`.
"""
function pair_table(ps::PairSet)
    rows = NamedTuple[]
    for (p, fp) in enumerate(ps.pairs)
        push!(rows, (
            pair_index = p,
            feature_a  = fp.name_i,
            feature_b  = fp.name_j,
        ))
    end
    return DataFrame(rows)
end
