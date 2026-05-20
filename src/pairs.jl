"""
Construction of feature-pair projections.

The pairwise decomposition is the heart of Platypus.jl. Rather than score
an entity's typicality in the full `n`-dimensional feature space — where
density estimates are fragile and uninterpretable — the collection is
projected onto each of its `n(n-1)/2` two-dimensional feature pairs.
Every pair is a small, well-populated cloud in which the focal entity's
position has a concrete, reportable meaning ("unusual for a song this
danceable to be this acoustic").

[`build_pairs`](@ref) first applies a per-feature `transform` and then
enumerates the pairs. The default `:rank` transform replaces each feature
with its *copula rank*, which strips out marginal rarity so the
subsequent typicality score reflects the *combination* of two features,
not one feature being extreme on its own.
"""

"""
    build_pairs(es::EntitySet; transform::Symbol=:rank) -> PairSet

Project an [`EntitySet`](@ref) onto all of its `n(n-1)/2` feature pairs.

# Keyword Arguments
- `transform::Symbol=:rank`: per-feature transform applied before pairs
  are formed.
  - `:rank` — replace each feature by its rank rescaled to `(0, 1)` via
    `(tiedrank(x) - 0.5) / n`. This is the empirical *copula* transform:
    every feature becomes marginally uniform, so a pair's cloud reflects
    only the dependence structure of the two features. Recommended, and
    the default.
  - `:zscore` — standardise each feature to zero mean and unit variance.
    Keeps the marginal shape; useful when extreme marginal values should
    still count as atypical.
  - `:none` — use the raw feature values unchanged.

Throws `ArgumentError` for any other `transform`, and for an
`EntitySet` with fewer than two features (no pair can be formed).
"""
function build_pairs(es::EntitySet; transform::Symbol=:rank)
    transform in (:rank, :zscore, :none) || throw(ArgumentError(
        "unknown transform :$transform; expected :rank, :zscore, or :none"))
    es.n_features < 2 && throw(ArgumentError(
        "build_pairs needs at least two features, got $(es.n_features)"))

    n, p = es.n_entities, es.n_features
    R = Matrix{Float64}(undef, n, p)
    for c in 1:p
        col = @view es.matrix[:, c]
        R[:, c] = transform === :rank   ? _rank_transform(col) :
                  transform === :zscore ? _zscore(col) :
                                          collect(col)
    end

    pairs = FeaturePair[]
    sizehint!(pairs, p * (p - 1) ÷ 2)
    for i in 1:(p - 1), j in (i + 1):p
        push!(pairs, FeaturePair(i, j, es.feature_names[i], es.feature_names[j],
                                 R[:, i], R[:, j]))
    end
    return PairSet(pairs, es.feature_names, es.entity_ids, es.labels,
                   n, length(pairs))
end

"""
    _rank_transform(x) -> Vector{Float64}

Empirical copula transform of a single feature: `(tiedrank(x) - 0.5) / n`,
mapping the values into `(0, 1)` with ties averaged.
"""
function _rank_transform(x::AbstractVector)
    n = length(x)
    return (tiedrank(x) .- 0.5) ./ n
end

"""
    _zscore(x) -> Vector{Float64}

Standardise a feature to zero mean and unit variance. A feature with zero
variance is returned centred (all zeros) rather than divided by zero.
"""
function _zscore(x::AbstractVector)
    m = mean(x)
    s = std(x)
    s == 0 && return zeros(Float64, length(x))
    return (x .- m) ./ s
end
