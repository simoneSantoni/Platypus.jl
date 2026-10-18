"""
    GraphTypicalityResult

Full-feature graph connectivity scores, distinct from pairwise typicality.
`adjacency` is the symmetric, sparse self-tuning affinity graph; `degree` is
its weighted degree; `stationary` is the selected degree-proportional stationary
distribution. `typicality` is its empirical mid-percentile and `atypicality` is
one minus that value. Disconnected graphs do not have a unique stationary
distribution; see `configuration.stationary_selection`.
"""
struct GraphTypicalityResult
    typicality::Vector{Float64}
    atypicality::Vector{Float64}
    stationary::Vector{Float64}
    degree::Vector{Float64}
    adjacency::SparseArrays.SparseMatrixCSC{Float64,Int}
    entity_ids::Vector{Int}
    labels::Vector{Union{String,Nothing}}
    n_entities::Int
    configuration::NamedTuple
end

atypicality(r::GraphTypicalityResult) = copy(r.atypicality)
analysis_metadata(r::GraphTypicalityResult) = deepcopy(r.configuration)

function summarytable(r::GraphTypicalityResult)
    DataFrame(entity_id=r.entity_ids, label=r.labels, typicality=r.typicality,
              atypicality=r.atypicality, stationary=r.stationary, degree=r.degree)
end

function scoretable(r::GraphTypicalityResult; raw::Bool=false)
    DataFrame(entity_id=r.entity_ids, typicality=raw ? r.stationary : r.typicality)
end

function Base.show(io::IO, r::GraphTypicalityResult)
    print(io, "GraphTypicalityResult($(r.n_entities) entities, ",
          "$(SparseArrays.nnz(r.adjacency) ÷ 2) edges, k=$(r.configuration.k), ",
          "$(r.configuration.n_components) components)")
end

# Express the scale-free exponent with one canonical exact rational key.
# Equal geometric relations therefore produce identical affinities even when
# their distances/radii have different floating square-root roundoff.
function _graph_affinity(search::_NeighborSearch, i::Int, j::Int, squared_radii)
    squared = _neighbor_squared(search, i, j)
    iszero(squared) && return 1.0
    left, right = squared_radii[i], squared_radii[j]
    (iszero(left) || iszero(right)) && return 0.0
    ratio_squared = squared^2 / (left * right)
    return exp(-_neighbor_sqrt(ratio_squared))
end

function _graph_component_count(adjacency)
    n = size(adjacency, 1)
    seen, stack = falses(n), Int[]
    components = 0
    rows = SparseArrays.rowvals(adjacency)
    for start in 1:n
        seen[start] && continue
        components += 1
        push!(stack, start)
        seen[start] = true
        while !isempty(stack)
            i = pop!(stack)
            for p in SparseArrays.nzrange(adjacency, i)
                j = rows[p]
                if !seen[j]
                    seen[j] = true
                    push!(stack, j)
                end
            end
        end
    end
    return components
end

"""
    fit_graph_typicality(es::EntitySet; k=0, transform=:rank, weights=nothing,
                         backend=:kdtree, max_edges=1_000_000, provenance=(;))

Fit a graph on **all selected features together**. Exact Euclidean k-nearest
neighborhoods include every boundary tie, then their union is symmetrized.
An edge has affinity `exp(-distance²/(σᵢ*σⱼ))`, with each local scale σ equal to
that vertex's kth distance (Zelnik-Manor and Perona, 2004). At zero local scale,
coincident points have affinity one and noncoincident points affinity zero.

Weighted degree divided by total degree is a stationary distribution of the
undirected random walk (Lovász, 1993). For disconnected graphs this explicitly
selects component mass proportional to component degree; it is not unique.
Isolated vertices have zero mass when any edges exist; their conceptual
transition is a self-loop, but no loop is inserted in `adjacency`. An edgeless
graph uses the uniform stationary distribution of these self-loops.

These scores describe graph connectivity, not a probability of anomaly or a
pairwise score. Self-tuning affinities adapt to density, so low density does not
necessarily imply low connectivity. The full-feature representation can retain
interactions absent from pairwise projections, but does not guarantee their
detection. Nonuniform reference weights are explicitly unsupported.

`k=0` selects `min(max(2, round(Int, sqrt(n))), n-1)`; explicit k must lie in
1:n-1. `max_edges` bounds stored undirected edges and throws if exact ties or
large k would exceed it. Both requested search backends currently use a conservative filtered exact
linear scan, because floating KD-tree pruning is not certified for exact binary
geometry. Configuration records the requested and effective backends. The graph
has O(nk) edges only when boundary ties are limited.
"""
function fit_graph_typicality(es::EntitySet; k::Int=0, transform::Symbol=:rank,
        weights::Union{ReferenceWeights,Nothing}=nothing, backend::Symbol=:kdtree,
        max_edges::Int=1_000_000, provenance::NamedTuple=(;))
    _validate_entityset(es)
    n = es.n_entities
    n >= 2 || throw(ArgumentError("graph typicality requires at least two observations"))
    k >= 0 || throw(ArgumentError("k must be nonnegative; zero selects the automatic value"))
    kk = k == 0 ? min(max(2, round(Int, sqrt(n))), n - 1) : k
    1 <= kk < n || throw(ArgumentError("k must satisfy 1 ≤ k < number of observations"))
    max_edges > 0 || throw(ArgumentError("max_edges must be positive"))
    w = weights === nothing ? fill(1.0 / n, n) : _validated_weights(weights, n)
    all(==(w[1]), w) || throw(ArgumentError("graph typicality supports uniform positive reference weights only"))
    points = if transform == :rank
        _neighbor_rank_coordinates(es.matrix)
    else
        columns = _transform_columns(es.matrix, transform)
        reduce(vcat, permutedims.(columns))
    end
    search = _neighbor_setup_points(points, backend)
    radii = _neighbor_kdistances(search, kk)
    squared_radii = [_neighbor_squared(search, i, search.boundary[i]) for i in 1:n]
    edges = Dict{Tuple{Int,Int},Float64}()
    indices = Int[]
    for i in 1:n
        _neighbor_indices!(indices, search, i, radii[i])
        for j in indices
            key = minmax(i, j)
            haskey(edges, key) && continue
            affinity = _graph_affinity(search, i, j, squared_radii)
            iszero(affinity) && continue
            length(edges) < max_edges || throw(ArgumentError(
                "exact graph exceeds max_edges=$max_edges (boundary ties can make it dense); reduce data/k or explicitly raise the limit"))
            edges[key] = affinity
        end
    end
    row, col, value = Int[], Int[], Float64[]
    sizehint!(row, 2length(edges)); sizehint!(col, 2length(edges)); sizehint!(value, 2length(edges))
    for ((i, j), affinity) in edges
        push!(row, i, j); push!(col, j, i); push!(value, affinity, affinity)
    end
    adjacency = SparseArrays.sparse(row, col, value, n, n)
    degree = [_stable_sum(view(SparseArrays.nonzeros(adjacency),
                              SparseArrays.nzrange(adjacency, i))) for i in 1:n]
    total = _stable_sum(degree)
    stationary = total == 0 ? fill(1.0 / n, n) : degree ./ total
    typicality = _score_percentiles(stationary, fill(1.0 / n, n))
    components = _graph_component_count(adjacency)
    configuration = (method=:graph, graph_family=:self_tuning_knn,
        transform=transform, feature_names=copy(es.feature_names), k=kk,
        requested_k=k, requested_backend=backend, backend=:brute,
        effective_backend=:brute, distance_policy=:filtered_exact,
        distance_predicate=:exact_binary_squared,
        geometry=transform == :rank ? :integer_midrank : :transformed_binary,
        entity_ids=copy(es.entity_ids),
        symmetrization=:union, boundary_ties=:inclusive, max_edges=max_edges,
        n_edges=length(edges), n_components=components,
        stationary_unique=components == 1,
        stationary_selection=total == 0 ? :uniform_edgeless : :degree_weighted_components,
        isolated_transition=:self_loop, zero_scale_affinity=:coincident_only,
        score_scale=:empirical_mid_percentile, raw_scale=:stationary_mass,
        package_version="0.2.0", normalization=:weighted_midcdf,
        reference_weights=copy(w), provenance=deepcopy(provenance))
    return GraphTypicalityResult(typicality, 1 .- typicality, stationary, degree,
        adjacency, copy(es.entity_ids), copy(es.labels), n, configuration)
end
