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
        all(isfinite, features) || throw(ArgumentError("feature values must be finite"))
        new(id, copy(features), label)
    end
end

function Entity(; id::Int,
                features::AbstractVector{<:Real},
                label::Union{AbstractString, Nothing}=nothing)
    Entity(id, Float64.(features), label === nothing ? nothing : String(label))
end

"""
    features(e::Entity) -> Vector{Float64}

Detached copy of the entity feature vector.
"""
features(e::Entity) = copy(e.features)

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
Platypus.jl scores an entity relative to the supplied members of the same
`EntitySet` (optionally re-weighted by [`ReferenceWeights`](@ref)).

# Fields
- `entities::Vector{Entity}`: Detached entity snapshots (a computed property)
- `feature_names::Vector{Symbol}`: Names of the features, in column order
- `matrix::Matrix{Float64}`: Authoritative `n_entities × n_features` value matrix
- `entity_ids::Vector{Int}`: Entity IDs, in row order
- `labels::Vector{Union{String, Nothing}}`: Entity labels, in row order
- `n_entities::Int`: Number of entities
- `n_features::Int`: Number of features

`matrix` is the authoritative value store. Constructors copy input data;
indexing, iteration, and the `entities` property return detached `Entity`
snapshots. Edit `matrix` explicitly to update a collection, then rebuild its
pairs. Analysis functions validate direct field edits before using them.
"""
struct EntitySet
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
        allunique(feature_names) || throw(ArgumentError("feature names must be unique"))
        n = length(entities)
        M = Matrix{Float64}(undef, n, p)
        for (i, e) in enumerate(entities)
            M[i, :] = e.features
        end
        ids = Int[e.id for e in entities]
        labs = Union{String, Nothing}[e.label for e in entities]
        allunique(ids) || throw(ArgumentError("entity IDs must be unique"))
        all(isfinite, M) || throw(ArgumentError("feature values must be finite"))
        new(copy(feature_names), M, ids, labs, n, p)
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
    length(idv) == n || throw(DimensionMismatch("ids must have one entry per row"))
    length(labv) == n || throw(DimensionMismatch("labels must have one entry per row"))
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
    state > es.n_entities ? nothing : (es[state], state + 1)
Base.getindex(es::EntitySet, i::Integer) = Entity(es.entity_ids[i], es.matrix[i, :], es.labels[i])
Base.getindex(es::EntitySet, ::Colon) = [es[i] for i in 1:es.n_entities]
Base.getindex(es::EntitySet, indices) = [es[i] for i in (1:es.n_entities)[indices]]
Base.eltype(::Type{EntitySet}) = Entity
function Base.getproperty(es::EntitySet, name::Symbol)
    name === :entities && return [es[i] for i in 1:es.n_entities]
    return getfield(es, name)
end
Base.propertynames(::EntitySet, private::Bool=false) =
    (:entities, :feature_names, :matrix, :entity_ids, :labels, :n_entities, :n_features)
Base.lastindex(es::EntitySet) = es.n_entities

# Direct field arrays may be edited by advanced callers; every analysis entry
# point validates them. Entity indexing/accessors return detached snapshots.
function _validate_entityset(es::EntitySet)
    size(es.matrix) == (es.n_entities, es.n_features) || throw(DimensionMismatch("entity matrix dimensions changed"))
    length(es.feature_names) == es.n_features || throw(DimensionMismatch("feature schema length changed"))
    length(es.entity_ids) == es.n_entities || throw(DimensionMismatch("entity ID length changed"))
    length(es.labels) == es.n_entities || throw(DimensionMismatch("entity label length changed"))
    allunique(es.entity_ids) || throw(ArgumentError("entity IDs must be unique"))
    allunique(es.feature_names) || throw(ArgumentError("feature names must be unique"))
    all(isfinite, es.matrix) || throw(ArgumentError("feature values must be finite"))
    return es
end

function Base.show(io::IO, es::EntitySet)
    print(io, "EntitySet($(es.n_entities) entities × $(es.n_features) features)")
end

"""
    FeaturePair

The projection of an [`EntitySet`](@ref) onto two of its features.

Each entity contributes one point `(u, v)` to a two-dimensional cloud.
By default the coordinates are empirical midranks in `(0, 1)`, removing
measurement units while preserving marginal order and ties. How boundaries
and dependence affect typicality depends on the scoring method —
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
- `transform::Symbol`: Preprocessing used by `build_pairs`, or `:unknown`
- `source_matrix::Union{Nothing, Matrix{Float64}}`: Detached raw features,
  allowing reference weights to be applied before preprocessing at fit time

Fitting validates that every unordered feature pair occurs exactly once;
the order of the pairs is unrestricted. Known transforms require source
provenance. Direct field edits are validated before analysis, but concurrent
mutation while an analysis is running is unsupported.
"""
struct PairSet
    pairs::Vector{FeaturePair}
    feature_names::Vector{Symbol}
    entity_ids::Vector{Int}
    labels::Vector{Union{String, Nothing}}
    n_entities::Int
    n_pairs::Int
    transform::Symbol
    source_matrix::Union{Nothing, Matrix{Float64}}
end

# Preserve the historical constructor for manually supplied coordinates. Their
# original marginal transform cannot be recovered, so provenance is unknown.
PairSet(pairs::Vector{FeaturePair}, names::Vector{Symbol}, ids::Vector{Int},
        labels::Vector{Union{String, Nothing}}, n::Int, np::Int) =
    PairSet(pairs, copy(names), copy(ids), copy(labels), n, np, :unknown, nothing)

function _validate_pairset(ps::PairSet)
    ps.n_entities > 0 || throw(ArgumentError("a PairSet needs at least one entity"))
    length(ps.pairs) == ps.n_pairs || throw(DimensionMismatch("pair count changed"))
    length(ps.entity_ids) == ps.n_entities || throw(DimensionMismatch("entity ID length changed"))
    length(ps.labels) == ps.n_entities || throw(DimensionMismatch("entity label length changed"))
    allunique(ps.entity_ids) || throw(ArgumentError("entity IDs must be unique"))
    allunique(ps.feature_names) || throw(ArgumentError("feature names must be unique"))
    p = length(ps.feature_names)
    ps.n_pairs == p * (p - 1) ÷ 2 || throw(ArgumentError("PairSet must contain every unordered feature pair exactly once"))
    allunique((fp.i, fp.j) for fp in ps.pairs) || throw(ArgumentError("PairSet contains duplicate feature pairs"))
    if ps.source_matrix === nothing
        ps.transform === :unknown || throw(ArgumentError("known transforms require source-matrix provenance; rebuild with build_pairs"))
    else
        ps.transform in (:rank, :none, :zscore) || throw(ArgumentError("source-backed PairSet requires a known transform"))
    end
    for fp in ps.pairs
        1 <= fp.i < fp.j <= length(ps.feature_names) || throw(ArgumentError("invalid feature pair indices"))
        (fp.name_i, fp.name_j) == (ps.feature_names[fp.i], ps.feature_names[fp.j]) ||
            throw(ArgumentError("pair names do not match feature schema"))
        length(fp.u) == length(fp.v) == ps.n_entities || throw(DimensionMismatch("pair coordinate lengths differ"))
        all(isfinite, fp.u) && all(isfinite, fp.v) || throw(ArgumentError("pair coordinates must be finite"))
    end
    if ps.source_matrix !== nothing
        size(ps.source_matrix) == (ps.n_entities, length(ps.feature_names)) ||
            throw(DimensionMismatch("source matrix dimensions differ"))
        all(isfinite, ps.source_matrix) || throw(ArgumentError("source features must be finite"))
        expected = _transform_columns(ps.source_matrix, ps.transform)
        for fp in ps.pairs
            fp.u == expected[fp.i] && fp.v == expected[fp.j] ||
                throw(ArgumentError("pair coordinates disagree with source data; rebuild with build_pairs"))
        end
    end
    return ps
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
recency, by a documented membership rule. They are not survey-design weights. A
weight of `0` excludes an entity from reference mass while retaining it as a
query where supported; uniform weights recover the plain
"the collection is its own reference" case.

Nonuniform weights are supported by the `:kde` and `:depth` methods of
[`fit_typicality`](@ref) and the Gaussian family of [`fit_copula`](@ref).
Other estimators reject nonuniform weights; see their documented contracts.

# Fields
- `weights::Vector{Float64}`: One non-negative weight per entity
"""
struct ReferenceWeights
    weights::Vector{Float64}

    function ReferenceWeights(weights::AbstractVector{<:Real})
        values = Float64.(weights)
        all(isfinite, values) || throw(ArgumentError("reference weights must be finite"))
        any(<(0), weights) &&
            throw(ArgumentError("reference weights must be non-negative"))
        any(i -> weights[i] > 0 && values[i] == 0, eachindex(weights)) &&
            throw(ArgumentError("positive reference weights must be representable as Float64"))
        all(==(0), values) &&
            throw(ArgumentError("at least one reference weight must be positive"))
        new(values)
    end
end

# Revalidate mutable storage at fit time and normalize without overflowing sums.
function _validated_weights(rw::ReferenceWeights, n::Integer; min_positive::Integer=1)
    min_positive >= 1 || throw(ArgumentError("min_positive must be at least one"))
    length(rw) == n || throw(DimensionMismatch("weights must have one entry per entity"))
    values = rw.weights
    all(isfinite, values) || throw(ArgumentError("reference weights must be finite"))
    all(>=(0), values) || throw(ArgumentError("reference weights must be non-negative"))
    count(>(0), values) >= min_positive || throw(ArgumentError("insufficient positive reference weights"))
    # Preserve the supplied binary weights until the single final conversion.
    # Floating reductions would make normalization depend on row order or on
    # zero-mass rows. Estimators retain the original weights for later geometry.
    exact = Rational{BigInt}.(values)
    normalized = _rational_float.(exact ./ sum(exact))
    any(i -> values[i] > 0 && normalized[i] == 0, eachindex(values)) &&
        throw(ArgumentError("reference weight dynamic range is not representable"))
    return normalized
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

The descriptive pairwise rules accepted by [`fit_typicality`](@ref).

- `:knn` — negative kth-neighbor distance; exact search, uniform weights.
- `:lof` — reciprocal tie-inclusive Local Outlier Factor; uniform weights.
- `:kde` — weighted leave-one-out log density: reflected compact kernels in
  rank space, Gaussian kernels on real coordinates.
- `:depth` — weighted global Mahalanobis centrality with an explicit singular
  affine-support convention. This is not a local density.

Gaussian and empirical beta copula families use [`fit_copula`](@ref); graph
connectivity uses [`fit_graph_typicality`](@ref).
"""
const TYPICALITY_METHODS = (:knn, :lof, :kde, :depth)

"""
    feature_names(es::EntitySet) -> Vector{Symbol}
    feature_names(ps::PairSet)   -> Vector{Symbol}

Detached copy of the feature names, in column order.
"""
feature_names(es::EntitySet) = copy(es.feature_names)
feature_names(ps::PairSet) = copy(ps.feature_names)

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
