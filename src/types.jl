"""
Core types for entity atypicality analysis.

An *entity* is a row of continuous feature values (a song with its audio
features, a firm with its financials, a paper with its bibliometric
scores). An [`EntitySet`](@ref) is the full collection of entities — and,
crucially, *its own reference set*: when no external population is
available, typicality is measured relative to the observed cloud itself.

A [`FeaturePair`](@ref) is the projection of the collection onto two
features, a two-dimensional point cloud in which one entity's position
can be scored for local typicality. A [`PairSet`](@ref) is the set of all
`n(n-1)/2` such projections, the object the estimators in
[`Platypus`](@ref) consume.
"""

"""
    Entity

A single entity: an `id`, a vector of continuous `features`, and an
optional `label`.

The `label` is a soft, possibly fuzzy category tag (a genre, a sector, a
field). It is never used to *define* the reference set — Platypus.jl
treats the reference set as endogenous — but it is carried through so
that results can be inspected by group.

# Fields
- `id::Int`: Identifier of the entity
- `features::Vector{Float64}`: Continuous feature values, aligned to the
  feature schema of the enclosing [`EntitySet`](@ref)
- `label::Union{String, Nothing}`: Optional category label (`nothing` if
  unlabelled)
"""
struct Entity
    id::Int
    features::Vector{Float64}
    label::Union{String, Nothing}

    function Entity(id::Int, features::Vector{Float64}, label::Union{String, Nothing})
        isempty(features) &&
            throw(ArgumentError("an Entity must have at least one feature"))
        new(id, features, label)
    end
end

function Entity(; id::Int,
                features::AbstractVector{<:Real},
                label::Union{AbstractString, Nothing}=nothing)
    Entity(id, Float64.(features), label === nothing ? nothing : String(label))
end

"""
    features(e::Entity) -> Vector{Float64}

Feature vector of the entity.
"""
features(e::Entity) = e.features

"""
    label(e::Entity) -> Union{String, Nothing}

Category label of the entity, or `nothing` if unlabelled.
"""
label(e::Entity) = e.label

function Base.show(io::IO, e::Entity)
    lab = e.label === nothing ? "" : ", label=\"$(e.label)\""
    print(io, "Entity(id=$(e.id), $(length(e.features)) features$lab)")
end

"""
    EntitySet

A collection of [`Entity`](@ref) records sharing one feature schema.

An `EntitySet` is the unit of analysis *and* the reference set: when the
boundary of the population is fuzzy or unknown, every estimator in
Platypus.jl scores an entity relative to the other members of the same
`EntitySet` (optionally re-weighted by [`ReferenceWeights`](@ref)).

# Fields
- `entities::Vector{Entity}`: The entities
- `feature_names::Vector{Symbol}`: Names of the features, in column order
- `matrix::Matrix{Float64}`: Cached `n_entities × n_features` value matrix
- `entity_ids::Vector{Int}`: Entity IDs, in row order
- `labels::Vector{Union{String, Nothing}}`: Entity labels, in row order
- `n_entities::Int`: Number of entities
- `n_features::Int`: Number of features
"""
struct EntitySet
    entities::Vector{Entity}
    feature_names::Vector{Symbol}
    matrix::Matrix{Float64}
    entity_ids::Vector{Int}
    labels::Vector{Union{String, Nothing}}
    n_entities::Int
    n_features::Int

    function EntitySet(entities::Vector{Entity}, feature_names::Vector{Symbol})
        isempty(entities) &&
            throw(ArgumentError("an EntitySet must contain at least one entity"))
        p = length(feature_names)
        p == 0 && throw(ArgumentError("an EntitySet must have at least one feature"))
        for e in entities
            length(e.features) == p || throw(DimensionMismatch(
                "entity $(e.id) has $(length(e.features)) features, " *
                "but $p feature names were given"))
        end
        n = length(entities)
        M = Matrix{Float64}(undef, n, p)
        for (i, e) in enumerate(entities)
            M[i, :] = e.features
        end
        ids = Int[e.id for e in entities]
        labs = Union{String, Nothing}[e.label for e in entities]
        new(entities, feature_names, M, ids, labs, n, p)
    end
end

"""
    EntitySet(matrix::AbstractMatrix, feature_names; ids=nothing, labels=nothing)

Build an `EntitySet` directly from an `n_entities × n_features` numeric
matrix. `ids` defaults to `1:n`; `labels` defaults to all `nothing`.
"""
function EntitySet(matrix::AbstractMatrix{<:Real}, feature_names::AbstractVector;
                   ids::Union{AbstractVector{<:Integer}, Nothing}=nothing,
                   labels::Union{AbstractVector, Nothing}=nothing)
    n, p = size(matrix)
    length(feature_names) == p || throw(DimensionMismatch(
        "matrix has $p columns but $(length(feature_names)) feature names were given"))
    idv = ids === nothing ? collect(1:n) : collect(Int.(ids))
    labv = labels === nothing ? fill(nothing, n) : labels
    ents = Entity[]
    for i in 1:n
        lab = labv[i]
        labstr = (lab === nothing || (lab isa Missing)) ? nothing : String(string(lab))
        push!(ents, Entity(idv[i], Float64.(matrix[i, :]), labstr))
    end
    EntitySet(ents, Symbol.(feature_names))
end

Base.length(es::EntitySet) = es.n_entities
Base.iterate(es::EntitySet, state=1) =
    state > es.n_entities ? nothing : (es.entities[state], state + 1)
Base.getindex(es::EntitySet, i) = es.entities[i]
Base.lastindex(es::EntitySet) = es.n_entities

function Base.show(io::IO, es::EntitySet)
    print(io, "EntitySet($(es.n_entities) entities × $(es.n_features) features)")
end

"""
    FeaturePair

The projection of an [`EntitySet`](@ref) onto two of its features.

Each entity contributes one point `(u, v)` to a two-dimensional cloud.
By default the coordinates are *copula ranks* (each feature replaced by
its rank, rescaled to `(0, 1)`) so that marginal rarity is removed and
only the *joint* configuration of the two features carries information —
see [`build_pairs`](@ref).

# Fields
- `i::Int`, `j::Int`: Column indices of the two features (`i < j`)
- `name_i::Symbol`, `name_j::Symbol`: Names of the two features
- `u::Vector{Float64}`: Coordinate of feature `i`, one value per entity
- `v::Vector{Float64}`: Coordinate of feature `j`, one value per entity
"""
struct FeaturePair
    i::Int
    j::Int
    name_i::Symbol
    name_j::Symbol
    u::Vector{Float64}
    v::Vector{Float64}
end

function Base.show(io::IO, fp::FeaturePair)
    print(io, "FeaturePair($(fp.name_i) ✕ $(fp.name_j), $(length(fp.u)) entities)")
end

"""
    PairSet

The set of all `n(n-1)/2` [`FeaturePair`](@ref) projections of an
[`EntitySet`](@ref). This is the object [`fit_typicality`](@ref)
consumes.

# Fields
- `pairs::Vector{FeaturePair}`: All feature pairs
- `feature_names::Vector{Symbol}`: Names of the underlying features
- `entity_ids::Vector{Int}`: Entity IDs, in row order
- `labels::Vector{Union{String, Nothing}}`: Entity labels, in row order
- `n_entities::Int`: Number of entities
- `n_pairs::Int`: Number of feature pairs
"""
struct PairSet
    pairs::Vector{FeaturePair}
    feature_names::Vector{Symbol}
    entity_ids::Vector{Int}
    labels::Vector{Union{String, Nothing}}
    n_entities::Int
    n_pairs::Int
end

Base.length(ps::PairSet) = ps.n_pairs
Base.iterate(ps::PairSet, state=1) =
    state > ps.n_pairs ? nothing : (ps.pairs[state], state + 1)
Base.getindex(ps::PairSet, i) = ps.pairs[i]
Base.lastindex(ps::PairSet) = ps.n_pairs

function Base.show(io::IO, ps::PairSet)
    print(io, "PairSet($(ps.n_pairs) pairs over $(length(ps.feature_names)) " *
              "features, $(ps.n_entities) entities)")
end

"""
    ReferenceWeights

Per-entity, non-negative weights describing graded membership in the
reference set.

When the boundary of the population is fuzzy, an entity is not simply
"in" or "out" of the comparison set. `ReferenceWeights` lets each entity
contribute to the reference cloud in proportion to its weight — by
recency, by similarity to a focal genre, by sampling probability. A
weight of `0` excludes an entity; uniform weights recover the plain
"the collection is its own reference" case.

Weights are honoured by the `:kde` and `:depth` methods of
[`fit_typicality`](@ref); see that function for details.

# Fields
- `weights::Vector{Float64}`: One non-negative weight per entity
"""
struct ReferenceWeights
    weights::Vector{Float64}

    function ReferenceWeights(weights::AbstractVector{<:Real})
        any(<(0), weights) &&
            throw(ArgumentError("reference weights must be non-negative"))
        all(==(0), weights) &&
            throw(ArgumentError("at least one reference weight must be positive"))
        new(Float64.(weights))
    end
end

"""
    ReferenceWeights(n::Integer)

Uniform reference weights for `n` entities (every entity contributes
equally) — the default reference model.
"""
ReferenceWeights(n::Integer) = ReferenceWeights(ones(Float64, n))

Base.length(rw::ReferenceWeights) = length(rw.weights)
Base.getindex(rw::ReferenceWeights, i) = rw.weights[i]

function Base.show(io::IO, rw::ReferenceWeights)
    print(io, "ReferenceWeights($(length(rw.weights)) entities, " *
              "$(count(>(0), rw.weights)) with positive weight)")
end

"""
    TYPICALITY_METHODS

The local-density scoring rules accepted by [`fit_typicality`](@ref).

- `:knn` — inverse distance to the `k`-th nearest neighbour. The simplest
  density estimate; fast and assumption-free.
- `:lof` — Local Outlier Factor (Breunig et al. 2000). Scores an entity's
  density *relative to its neighbours' density*, so it adapts to a cloud
  whose density varies from region to region — the natural choice for a
  heterogeneous set with no uniform scale.
- `:kde` — two-dimensional Gaussian kernel density. Smooth; honours
  [`ReferenceWeights`](@ref).
- `:depth` — Mahalanobis data depth. A centre-outward ordering of the
  cloud that needs no distributional assumption; honours
  [`ReferenceWeights`](@ref).
"""
const TYPICALITY_METHODS = (:knn, :lof, :kde, :depth)

"""
    feature_names(es::EntitySet) -> Vector{Symbol}
    feature_names(ps::PairSet)   -> Vector{Symbol}

Names of the features, in column order.
"""
feature_names(es::EntitySet) = es.feature_names
feature_names(ps::PairSet) = ps.feature_names

"""
    n_entities(es::EntitySet) -> Int
    n_entities(ps::PairSet)   -> Int

Number of entities.
"""
n_entities(es::EntitySet) = es.n_entities
n_entities(ps::PairSet) = ps.n_entities

"""
    n_features(es::EntitySet) -> Int

Number of features.
"""
n_features(es::EntitySet) = es.n_features

"""
    n_pairs(ps::PairSet) -> Int

Number of feature pairs, equal to `n_features * (n_features - 1) ÷ 2`.
"""
n_pairs(ps::PairSet) = ps.n_pairs
