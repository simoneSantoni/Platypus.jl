# Filtered exact Euclidean predicates. No n×n distance matrix is formed.
# A rounded norm is an output, never a boundary-membership predicate.
struct _NeighborSearch
    points::Matrix{Float64}
    exact_points::Matrix{Rational{BigInt}}
    tree::Nothing
    scale_exponent::Int
    boundary::Vector{Int}
    requested_backend::Symbol
end

function _neighbor_setup(u::AbstractVector, v::AbstractVector, k::Int, backend::Symbol)
    n = length(u)
    length(v) == n || throw(DimensionMismatch("neighbor coordinate lengths differ"))
    n >= 2 || throw(ArgumentError("neighbor scoring requires at least two observations"))
    1 <= k < n || throw(ArgumentError("k must satisfy 1 ≤ k < number of observations"))
    points = Matrix{Float64}(undef, 2, n)
    points[1, :] = u
    points[2, :] = v
    return _neighbor_setup_points(points, backend)
end

function _neighbor_setup_points(points::Matrix{Float64}, backend::Symbol)
    backend in (:kdtree, :brute) || throw(ArgumentError("backend must be :kdtree or :brute"))
    maxexp, minexp = typemin(Int), typemax(Int)
    for y in points
        isfinite(y) || throw(ArgumentError("neighbor coordinates must be finite Float64 values"))
        if !iszero(y)
            e = exponent(abs(y))
            maxexp, minexp = max(maxexp, e), min(minexp, e)
        end
    end
    # This common power-of-two scaling is exact: the bounded exponent spread
    # prevents underflow. Wider ranges retain the original binary coordinates.
    extreme = maxexp != typemin(Int) && maxexp - minexp > 400
    shift = maxexp == typemin(Int) || extreme ? 0 : maxexp
    scaled = shift == 0 ? copy(points) : ldexp.(points, -shift)
    # KD-tree floating-point pruning has no certified bound here. Until it does,
    # both requested backends conservatively scan every candidate, then apply
    # the filtered exact predicate. Metadata exposes this fallback explicitly.
    return _NeighborSearch(scaled, Rational{BigInt}.(scaled), nothing, shift,
                           zeros(Int, size(points, 2)), backend)
end

function _neighbor_squared(s::_NeighborSearch, i::Int, j::Int)
    result = BigInt(0) // BigInt(1)
    for c in axes(s.points, 1)
        delta = s.exact_points[c, i] - s.exact_points[c, j]
        result += delta * delta
    end
    return result
end

# Enclose the exact squared distance using outward rounding after EACH operation.
# This is an error filter, not a tolerance defining equality. Subnormal rounding
# and overflow simply widen the enclosure and trigger the exact fallback.
function _neighbor_squared_interval(s::_NeighborSearch, i::Int, j::Int)
    lower, upper = 0.0, 0.0
    for c in axes(s.points, 1)
        x, y = s.points[c, i], s.points[c, j]
        x == y && continue
        delta = abs(x - y)
        lo, hi = max(0.0, prevfloat(delta)), nextfloat(delta)
        lower = max(0.0, prevfloat(lower + max(0.0, prevfloat(lo * lo))))
        upper = nextfloat(upper + nextfloat(hi * hi))
    end
    return lower, upper
end

@inline function _neighbor_compare(s::_NeighborSearch, i::Int, j::Int, h::Int,
                                    jbounds, hbounds)
    jbounds[1] == jbounds[2] && hbounds[1] == hbounds[2] &&
        return cmp(jbounds[1], hbounds[1])
    jbounds[2] < hbounds[1] && return -1
    hbounds[2] < jbounds[1] && return 1
    return cmp(_neighbor_squared(s, i, j), _neighbor_squared(s, i, h))
end

_neighbor_sqrt(q::Rational{BigInt}) = _rational_sqrt(q)

_neighbor_distance(s::_NeighborSearch, i::Int, j::Int) =
    _neighbor_sqrt(_neighbor_squared(s, i, j))

function _neighbor_indices!(indices::Vector{Int}, s::_NeighborSearch,
                            i::Int, radius::Float64)
    empty!(indices)
    boundary = s.boundary[i]
    boundary != 0 || throw(ArgumentError("compute kth distances before selecting neighbors"))
    bound = _neighbor_squared_interval(s, i, boundary)
    for j in axes(s.points, 2)
        i == j && continue
        candidate = _neighbor_squared_interval(s, i, j)
        _neighbor_compare(s, i, j, boundary, candidate, bound) <= 0 && push!(indices, j)
    end
    return indices
end

function _neighbor_kdistances(s::_NeighborSearch, k::Int)
    n = size(s.points, 2)
    1 <= k < n || throw(ArgumentError("k must satisfy 1 ≤ k < number of observations"))
    radii = Vector{Float64}(undef, n)
    candidates = Vector{Int}(undef, n - 1)
    bounds = Vector{Tuple{Float64,Float64}}(undef, n)
    for i in 1:n
        position = 0
        for j in 1:n
            i == j && continue
            position += 1
            candidates[position] = j
            bounds[j] = _neighbor_squared_interval(s, i, j)
        end
        boundary = partialsort!(candidates, k;
            lt=(j, h) -> _neighbor_compare(s, i, j, h, bounds[j], bounds[h]) < 0)
        s.boundary[i] = boundary
        radii[i] = _neighbor_distance(s, i, boundary)
        isfinite(radii[i]) || throw(ArgumentError("neighbor distances exceed Float64 range; rescale the features"))
    end
    return radii
end

"""Internal exact kth nonself Euclidean distance; both requested backends use a certified filtered scan."""
function _neighbor_kdistance(u::AbstractVector, v::AbstractVector, k::Int;
                             backend::Symbol=:kdtree)
    s = _neighbor_setup(u, v, k, backend)
    radii = _neighbor_kdistances(s, k)
    if s.scale_exponent != 0
        power = BigInt(1) << abs(2s.scale_exponent)
        for i in eachindex(radii)
            squared = _neighbor_squared(s, i, s.boundary[i])
            squared = s.scale_exponent > 0 ? squared * power : squared / power
            # Rescale the exact square before rounding, avoiding double rounding
            # when the output distance is subnormal.
            radii[i] = _neighbor_sqrt(squared)
        end
    end
    all(isfinite, radii) || throw(ArgumentError("neighbor distances exceed Float64 range; rescale the features"))
    return radii
end

# Stable against permutation and overflow of an otherwise representable mean.
function _neighbor_mean!(values::Vector{Float64})
    sort!(values)
    scale = last(values)
    iszero(scale) && return 0.0
    isinf(scale) && return Inf
    total = 0.0
    for x in values
        total += x / scale
    end
    return (total / length(values)) * scale
end

"""
Internal LOF using every neighbor at or inside the kth nonself distance.

For zero-reachability duplicate clusters we explicitly extend LOF to 1.
A positive-reachability observation with a zero-reachability neighbor has LOF
Inf. The original paper excludes this singular case; no epsilon or jitter is
introduced here. Inverse LOF maps the latter case to typicality zero.
"""
function _neighbor_lof(u::AbstractVector, v::AbstractVector, k::Int;
                       backend::Symbol=:kdtree)
    s = _neighbor_setup(u, v, k, backend)
    radii = _neighbor_kdistances(s, k)
    n = length(radii)
    reach_mean, lof = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
    indices, values = Int[], Float64[]
    for i in 1:n
        if iszero(radii[i])
            reach_mean[i] = 0.0
            continue
        end
        _neighbor_indices!(indices, s, i, radii[i])
        empty!(values)
        for j in indices
            push!(values, max(radii[j], _neighbor_distance(s, i, j)))
        end
        reach_mean[i] = _neighbor_mean!(values)
        isfinite(reach_mean[i]) || throw(ArgumentError("reachability distances exceed Float64 range; rescale the features"))
    end
    for i in 1:n
        if iszero(reach_mean[i])
            # All members of this zero-radius neighborhood are identical and
            # have zero reachability, so their density ratios are assigned 1.
            lof[i] = 1.0
            continue
        end
        _neighbor_indices!(indices, s, i, radii[i])
        empty!(values)
        for j in indices
            push!(values, reach_mean[i] / reach_mean[j])
        end
        lof[i] = _neighbor_mean!(values)
    end
    return lof
end
