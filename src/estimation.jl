"""
Estimation of entity atypicality.

This file contains the implemented estimator, [`fit_typicality`](@ref),
which scores every entity's local typicality on each feature pair and
aggregates the `n(n-1)/2` scores per entity into *conventionality*,
*novelty*, and *dispersion* — the summary triple of Uzzi et al. (2013).

More elaborate estimators — an endogenous graph random-walk model for
the genuinely fuzzy-reference case, a regression of an outcome on the
typicality summaries, and a parametric copula fit — are exposed by name
in [`Platypus`](@ref) but raise `not implemented` until added in a later
revision.
"""

# ----------------------------------------------------------------------
# Result type
# ----------------------------------------------------------------------

"""
    TypicalityResult

Result of [`fit_typicality`](@ref).

The core output is `pair_scores`, an `n_entities × n_pairs` matrix of
local typicality scores (high = typical). Each entity's row is then
summarised into the Uzzi-style triple:

- `conventionality` — the median pairwise score: how typical the entity
  is *overall*.
- `novelty` — a low quantile of the pairwise scores: how striking the
  entity's *most unusual* feature combination is.
- `dispersion` — the inter-quartile range of the pairwise scores: whether
  the entity mixes conventional and unconventional pairs.
- `atypicality` — `1 - conventionality`, provided for convenience.

# Fields
- `pair_scores::Matrix{Float64}`: `n_entities × n_pairs` typicality scores
- `conventionality::Vector{Float64}`: Median pairwise score, per entity
- `novelty::Vector{Float64}`: `novelty_quantile` of the pairwise scores
- `dispersion::Vector{Float64}`: IQR of the pairwise scores, per entity
- `atypicality::Vector{Float64}`: `1 .- conventionality`
- `pair_names::Vector{Tuple{Symbol, Symbol}}`: Feature pair of each column
- `entity_ids::Vector{Int}`: Entity IDs, in row order
- `labels::Vector{Union{String, Nothing}}`: Entity labels, in row order
- `method::Symbol`: Scoring rule used (one of [`TYPICALITY_METHODS`](@ref))
- `k::Int`: Neighbourhood size (for `:knn` / `:lof`)
- `novelty_quantile::Float64`: Quantile used for the `novelty` field
- `normalized::Bool`: Whether per-pair scores were rank-normalised
- `n_entities::Int`: Number of entities
- `n_pairs::Int`: Number of feature pairs
"""
struct TypicalityResult
    pair_scores::Matrix{Float64}
    conventionality::Vector{Float64}
    novelty::Vector{Float64}
    dispersion::Vector{Float64}
    atypicality::Vector{Float64}
    pair_names::Vector{Tuple{Symbol, Symbol}}
    entity_ids::Vector{Int}
    labels::Vector{Union{String, Nothing}}
    method::Symbol
    k::Int
    novelty_quantile::Float64
    normalized::Bool
    n_entities::Int
    n_pairs::Int
end

function Base.show(io::IO, r::TypicalityResult)
    println(io, "TypicalityResult")
    println(io, "  method      : ", r.method,
                r.method in (:knn, :lof) ? " (k = $(r.k))" : "")
    println(io, "  n_entities  : ", r.n_entities)
    println(io, "  n_pairs     : ", r.n_pairs)
    println(io, "  normalized  : ", r.normalized)
    println(io, "  atypicality : ",
            @sprintf("min %.3f / mean %.3f / max %.3f",
                     minimum(r.atypicality), mean(r.atypicality),
                     maximum(r.atypicality)))
end

"""
    atypicality(r::TypicalityResult) -> Vector{Float64}

Per-entity atypicality, `1 - conventionality`. This is the package's
headline quantity: higher means the entity sits, on average across its
feature pairs, in lower-density (more unusual) regions of the collection.
"""
atypicality(r::TypicalityResult) = r.atypicality

"""
    conventionality(r::TypicalityResult) -> Vector{Float64}

Per-entity conventionality: the median typicality across the entity's
feature pairs.
"""
conventionality(r::TypicalityResult) = r.conventionality

"""
    novelty(r::TypicalityResult) -> Vector{Float64}

Per-entity novelty: the `novelty_quantile` (default `0.1`) of the
entity's pairwise typicality scores — the signature of its single most
unusual feature combination.
"""
novelty(r::TypicalityResult) = r.novelty

"""
    dispersion(r::TypicalityResult) -> Vector{Float64}

Per-entity dispersion: the inter-quartile range of the entity's pairwise
typicality scores. A high value flags an entity that combines very
conventional pairs with very unconventional ones.
"""
dispersion(r::TypicalityResult) = r.dispersion

"""
    summarytable(r::TypicalityResult) -> DataFrame

Per-entity summary, one row per entity. Columns: `entity_id`, `label`,
`conventionality`, `novelty`, `dispersion`, `atypicality`.
"""
function summarytable(r::TypicalityResult)
    return DataFrame(
        entity_id       = r.entity_ids,
        label           = r.labels,
        conventionality = r.conventionality,
        novelty         = r.novelty,
        dispersion      = r.dispersion,
        atypicality     = r.atypicality,
    )
end

"""
    scoretable(r::TypicalityResult) -> DataFrame

Long-form table of the raw pairwise scores, one row per
entity × feature-pair. Columns: `entity_id`, `feature_a`, `feature_b`,
`typicality`. Suitable for joining against entity-level covariates or
for plotting an entity's pairwise typicality profile.
"""
function scoretable(r::TypicalityResult)
    rows = NamedTuple[]
    for i in 1:r.n_entities, p in 1:r.n_pairs
        push!(rows, (
            entity_id  = r.entity_ids[i],
            feature_a  = r.pair_names[p][1],
            feature_b  = r.pair_names[p][2],
            typicality = r.pair_scores[i, p],
        ))
    end
    return DataFrame(rows)
end

# ----------------------------------------------------------------------
# Per-pair local-density scoring rules
# ----------------------------------------------------------------------

# Pairwise Euclidean distance matrix for a 2-D cloud.
function _dist_matrix(u::Vector{Float64}, v::Vector{Float64})
    n = length(u)
    D = zeros(Float64, n, n)
    for a in 1:n, b in (a + 1):n
        d = hypot(u[a] - u[b], v[a] - v[b])
        D[a, b] = d
        D[b, a] = d
    end
    return D
end

# :knn — inverse distance to the k-th nearest neighbour.
function _knn_density(D::Matrix{Float64}, k::Int)
    n = size(D, 1)
    dens = Vector{Float64}(undef, n)
    for a in 1:n
        d = sort(D[a, :])          # d[1] == 0 (self)
        kd = d[min(k + 1, n)]
        dens[a] = 1.0 / (kd + 1e-12)
    end
    return dens
end

# :lof — Local Outlier Factor (Breunig, Kriegel, Ng & Sander 2000).
# Returns the LOF score itself (≈ 1 inlier, ≫ 1 outlier).
function _lof(D::Matrix{Float64}, k::Int)
    n = size(D, 1)
    kk = clamp(k, 1, n - 1)
    neigh = Vector{Vector{Int}}(undef, n)
    kdist = Vector{Float64}(undef, n)
    for a in 1:n
        order = sortperm(D[a, :])
        nb = Int[o for o in order if o != a]
        neigh[a] = nb[1:kk]
        kdist[a] = D[a, neigh[a][end]]
    end
    lrd = Vector{Float64}(undef, n)
    for a in 1:n
        s = 0.0
        for o in neigh[a]
            s += max(kdist[o], D[a, o])
        end
        lrd[a] = length(neigh[a]) / (s + 1e-12)
    end
    lof = Vector{Float64}(undef, n)
    for a in 1:n
        s = 0.0
        for o in neigh[a]
            s += lrd[o] / (lrd[a] + 1e-12)
        end
        lof[a] = s / length(neigh[a])
    end
    return lof
end

# Silverman rule-of-thumb bandwidth for one dimension.
function _silverman(x::AbstractVector, n::Int)
    s = std(x)
    s = s == 0 ? 1e-6 : s
    h = 1.06 * s * n^(-1 / 5)
    return h <= 0 ? 1e-6 : h
end

# :kde — weighted 2-D Gaussian kernel density (leave-one-out).
function _kde_density(u::Vector{Float64}, v::Vector{Float64}, w::Vector{Float64})
    n = length(u)
    hu = _silverman(u, n)
    hv = _silverman(v, n)
    dens = Vector{Float64}(undef, n)
    for a in 1:n
        s = 0.0
        wsum = 0.0
        for b in 1:n
            b == a && continue
            du = (u[a] - u[b]) / hu
            dv = (v[a] - v[b]) / hv
            s += w[b] * exp(-0.5 * (du * du + dv * dv))
            wsum += w[b]
        end
        dens[a] = s / (wsum * 2π * hu * hv + 1e-12)
    end
    return dens
end

# :depth — weighted Mahalanobis data depth, depth = 1 / (1 + dᴹ²).
function _maha_depth(u::Vector{Float64}, v::Vector{Float64}, w::Vector{Float64})
    n = length(u)
    wsum = sum(w)
    mu = sum(w .* u) / wsum
    mv = sum(w .* v) / wsum
    c11 = c12 = c22 = 0.0
    for a in 1:n
        du = u[a] - mu
        dv = v[a] - mv
        c11 += w[a] * du * du
        c12 += w[a] * du * dv
        c22 += w[a] * dv * dv
    end
    c11 = c11 / wsum + 1e-9
    c22 = c22 / wsum + 1e-9
    c12 = c12 / wsum
    det = c11 * c22 - c12 * c12
    depth = Vector{Float64}(undef, n)
    for a in 1:n
        du = u[a] - mu
        dv = v[a] - mv
        q = (c22 * du * du - 2 * c12 * du * dv + c11 * dv * dv) / (det + 1e-12)
        depth[a] = 1.0 / (1.0 + max(q, 0.0))
    end
    return depth
end

# Score one feature pair: raw local density, optionally rank-normalised
# so that the scores of different pairs are comparable.
function _pair_typicality(fp::FeaturePair, method::Symbol, k::Int,
                          w::Vector{Float64}, normalize::Bool)
    raw = if method === :knn
        _knn_density(_dist_matrix(fp.u, fp.v), k)
    elseif method === :lof
        1.0 ./ (_lof(_dist_matrix(fp.u, fp.v), k) .+ 1e-12)
    elseif method === :kde
        _kde_density(fp.u, fp.v, w)
    elseif method === :depth
        _maha_depth(fp.u, fp.v, w)
    else
        throw(ArgumentError(
            "unknown method :$method; expected one of $TYPICALITY_METHODS"))
    end
    normalize || return raw
    return tiedrank(raw) ./ length(raw)
end

# ----------------------------------------------------------------------
# fit_typicality
# ----------------------------------------------------------------------

"""
    fit_typicality(ps::PairSet; kwargs...) -> TypicalityResult

Score every entity's atypicality from the local density of its position
in each feature-pair cloud, then aggregate the `n(n-1)/2` pairwise scores
per entity.

The reference set is the [`PairSet`](@ref) itself: each entity is scored
relative to the other members of the same collection, so no externally
defined population is required.

# Keyword Arguments
- `method::Symbol=:lof`: local-density scoring rule, one of
  [`TYPICALITY_METHODS`](@ref) (`:knn`, `:lof`, `:kde`, `:depth`).
- `k::Integer=0`: neighbourhood size for `:knn` and `:lof`. The default
  `0` selects `round(Int, sqrt(n_entities))`, clamped to `1:(n-1)`.
- `novelty_quantile::Real=0.1`: quantile of the pairwise scores stored in
  the `novelty` field of the result.
- `normalize::Bool=true`: if `true`, each pair's raw density scores are
  replaced by their ranks rescaled to `(0, 1]`, making scores comparable
  across pairs before they are aggregated. Recommended.
- `weights::Union{ReferenceWeights, Nothing}=nothing`: graded membership
  in the reference set. Honoured by `:kde` and `:depth`; the
  neighbour-based rules `:knn` and `:lof` ignore it (a warning-free
  no-op) because they have no natural weighted form.

# Aggregation
For each entity, the pairwise scores are summarised into
`conventionality` (median), `novelty` (the `novelty_quantile`),
`dispersion` (inter-quartile range), and `atypicality`
(`1 - conventionality`). This is the summary design of
Uzzi et al. (2013), generalised from discrete co-citation `z`-scores to
continuous-feature local densities.

Throws `ArgumentError` if the `PairSet` has no pairs (fewer than two
features) or if `method` is not recognised.
"""
function fit_typicality(ps::PairSet;
                        method::Symbol=:lof,
                        k::Integer=0,
                        novelty_quantile::Real=0.1,
                        normalize::Bool=true,
                        weights::Union{ReferenceWeights, Nothing}=nothing)
    method in TYPICALITY_METHODS || throw(ArgumentError(
        "unknown method :$method; expected one of $TYPICALITY_METHODS"))
    ps.n_pairs == 0 && throw(ArgumentError(
        "fit_typicality needs at least one feature pair (two or more features)"))
    0 <= novelty_quantile <= 1 ||
        throw(ArgumentError("novelty_quantile must lie in [0, 1]"))

    n = ps.n_entities
    kk = k > 0 ? Int(k) : clamp(round(Int, sqrt(n)), 1, max(n - 1, 1))

    w = weights === nothing ? ones(Float64, n) : weights.weights
    length(w) == n || throw(DimensionMismatch(
        "weights has length $(length(w)) but the PairSet has $n entities"))

    scores = Matrix{Float64}(undef, n, ps.n_pairs)
    for (p, fp) in enumerate(ps.pairs)
        scores[:, p] = _pair_typicality(fp, method, kk, w, normalize)
    end

    conv = Vector{Float64}(undef, n)
    nov  = Vector{Float64}(undef, n)
    disp = Vector{Float64}(undef, n)
    for i in 1:n
        row = scores[i, :]
        conv[i] = median(row)
        nov[i]  = quantile(row, novelty_quantile)
        disp[i] = ps.n_pairs > 1 ? iqr(row) : 0.0
    end
    atyp = 1.0 .- conv

    pair_names = Tuple{Symbol, Symbol}[(fp.name_i, fp.name_j) for fp in ps.pairs]
    return TypicalityResult(scores, conv, nov, disp, atyp, pair_names,
                            ps.entity_ids, ps.labels, method, kk,
                            Float64(novelty_quantile), normalize,
                            n, ps.n_pairs)
end

# ----------------------------------------------------------------------
# Reserved estimators (not yet implemented)
# ----------------------------------------------------------------------

"""
    fit_graph_typicality(es::EntitySet; kwargs...)

Endogenous, graph-based typicality for the genuinely fuzzy-reference
case. Builds a similarity graph over entities and takes typicality to be
the stationary distribution of a random walk on it — the principled form
of the "emergent, self-organised, no imposed reference set" idea. Not yet
implemented.
"""
function fit_graph_typicality(::EntitySet; kwargs...)
    throw(ErrorException("fit_graph_typicality is not implemented yet"))
end

"""
    fit_outcome(r::TypicalityResult, y::AbstractVector; kwargs...)

Regress an entity-level outcome `y` (e.g. popularity, impact) on the
typicality summaries — `conventionality`, `novelty`, `dispersion` — to
test the Uzzi-style hypothesis that the sweet spot is high conventionality
paired with tail novelty. Not yet implemented.
"""
function fit_outcome(::TypicalityResult, ::AbstractVector; kwargs...)
    throw(ErrorException("fit_outcome is not implemented yet"))
end

"""
    fit_copula(ps::PairSet; kwargs...)

Fit a parametric copula to each feature pair and score typicality from
the fitted copula density, as a smooth, marginal-free alternative to the
non-parametric rules in [`fit_typicality`](@ref). Not yet implemented.
"""
function fit_copula(::PairSet; kwargs...)
    throw(ErrorException("fit_copula is not implemented yet"))
end
