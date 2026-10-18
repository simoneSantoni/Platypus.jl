"""
Construction of two-feature projections. Pairwise scoring is inspectable, but
bivariate densities do not identify the full joint distribution. Rank preprocessing removes
measurement units, not finite-sample boundary effects or the need to specify
which notion of typicality is being estimated. Ties remain probability masses.
"""

"""
    build_pairs(es::EntitySet; transform::Symbol=:rank) -> PairSet

Project an [`EntitySet`](@ref) onto all of its `n(n-1)/2` feature pairs.

# Keyword Arguments
- `transform::Symbol=:rank`: per-feature transform applied before pairs
  are formed.
  - `:rank` — replace each feature by its rank rescaled to `(0, 1)` via
    `(tiedrank(x) - 0.5) / n`. This is the empirical *copula* transform:
    continuous untied features become marginally uniform. Tied features
    retain probability masses, so constant or discrete features require care. Recommended, and
    the default.
  - `:zscore` — standardise each feature to zero mean and unit variance.
    Keeps the marginal shape; useful when extreme marginal values should
    still count as atypical.
  - `:none` — use the raw feature values unchanged.

Pair coordinates share one vector per feature. The returned object keeps a
detached raw-data snapshot and transform provenance for weighted fitting.
Direct field mutation is an advanced operation; rebuild pairs after editing
features. `EntitySet.matrix` is authoritative; entity accessors are snapshots.

Throws `ArgumentError` for any other `transform`, and for an
`EntitySet` with fewer than two features (no pair can be formed).
"""
function build_pairs(es::EntitySet; transform::Symbol=:rank)
    _validate_entityset(es)
    transform in (:rank, :zscore, :none) || throw(ArgumentError(
        "unknown transform :$transform; expected :rank, :zscore, or :none"))
    es.n_features < 2 && throw(ArgumentError(
        "build_pairs needs at least two features, got $(es.n_features)"))

    n, p = es.n_entities, es.n_features
    columns = _transform_columns(es.matrix, transform)

    pairs = FeaturePair[]
    sizehint!(pairs, p * (p - 1) ÷ 2)
    for i in 1:(p - 1), j in (i + 1):p
        push!(pairs, FeaturePair(i, j, es.feature_names[i], es.feature_names[j],
                                 columns[i], columns[j]))
    end
    return PairSet(pairs, copy(es.feature_names), copy(es.entity_ids), copy(es.labels),
                   n, length(pairs), transform, copy(es.matrix))
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
    length(x) <= 1 && return zeros(Float64, length(x))
    return _weighted_zscore(x, ones(Float64,length(x)))
end

# One contiguous coordinate vector per feature, shared by its FeaturePairs.
# With weights, ties receive the midpoint of their reference probability mass.
# Zero-weight rows are queries only and may map to the CDF boundaries 0 or 1.
function _transform_columns(matrix::AbstractMatrix, transform::Symbol, weights=nothing)
    transform in (:rank, :zscore, :none) || throw(ArgumentError("unknown feature transform"))
    all(isfinite, matrix) || throw(ArgumentError("feature values must be finite"))
    w = weights === nothing ? nothing : Float64.(weights)
    w === nothing || _validated_weights(ReferenceWeights(w),size(matrix,1))
    columns = Vector{Vector{Float64}}(undef, size(matrix, 2))
    for c in axes(matrix, 2)
        x = @view matrix[:, c]
        columns[c] = if transform === :none
            collect(x)
        elseif transform === :rank
            w === nothing ? _rank_transform(x) : _weighted_rank_transform(x, w)
        else
            w === nothing ? _zscore(x) : _weighted_zscore(x, w)
        end
        all(isfinite, columns[c]) || throw(ArgumentError("feature transform produced nonfinite coordinates"))
    end
    return columns
end

function _weighted_rank_transform(x::AbstractVector, w::Vector{Float64})
    return _rational_float.(_exact_rank_transform(x,w))
end

# Retain the intended empirical mass geometry until the numerical estimator
# needs floating values. In particular reflection maps r exactly to 1-r.
function _exact_rank_transform(x::AbstractVector, weights::AbstractVector)
    w=Rational{BigInt}.(weights)
    order = sortperm(x)
    result = zeros(Rational{BigInt}, length(x))
    total = sum(w)
    total > 0 || throw(ArgumentError("positive reference mass is required"))
    cumulative = zero(total)
    first = 1
    while first <= length(order)
        last = first
        while last < length(order) && x[order[last + 1]] == x[order[first]]
            last += 1
        end
        mass = sum(w[order[j]] for j in first:last)
        rank = (cumulative + mass / 2) / total
        for j in first:last
            result[order[j]] = rank
        end
        cumulative += mass
        first = last + 1
    end
    return result
end

function _weighted_zscore(x::AbstractVector, w::Vector{Float64})
    support = findall(>(0), w)
    isempty(support) && throw(ArgumentError("positive reference mass is required"))
    values=Rational{BigInt}.(x)
    p=Rational{BigInt}.(w); p ./= sum(p)
    mu=sum(p[i]*values[i] for i in support)
    denom=1-sum(abs2,p)
    variance=denom>0 ? sum(p[i]*(values[i]-mu)^2 for i in support)/denom : zero(mu)
    # Exact centered moments avoid losing small deviations beside a large
    # offset. Round the dimensionless squared ratio only after cancellation.
    if iszero(variance)
        return _rational_float.(values .- mu)
    end
    [v < mu ? -_rational_sqrt((v-mu)^2/variance) :
        _rational_sqrt((v-mu)^2/variance) for v in values]
end

function _exact_transform_columns(matrix::AbstractMatrix, transform::Symbol, w)
    transform in (:rank,:none) || throw(ArgumentError("exact coordinate provenance supports rank or raw coordinates"))
    [transform===:rank ? _exact_rank_transform(view(matrix,:,c),w) :
        Rational{BigInt}.(view(matrix,:,c)) for c in axes(matrix,2)]
end

# Common translation and positive scale cancel from neighbor membership and
# LOF. kNN distances are converted back by the caller using the factor 2n.
function _neighbor_rank_coordinates(matrix::AbstractMatrix)
    points=Matrix{Float64}(undef,size(matrix,2),size(matrix,1))
    for c in axes(matrix,2)
        points[c,:]=2 .* tiedrank(view(matrix,:,c))
    end
    points
end
