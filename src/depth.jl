# Two-dimensional Mahalanobis depth admits an exact rational solution. Binary
# floating inputs (and exact rank coordinates) are rational numbers, so exact
# moments avoid cancellation, numerical rank decisions, and broken score ties.
# Conversion to Float64 occurs only after the depth has been evaluated.
function _maha_depth(u::AbstractVector{<:Real}, v::AbstractVector{<:Real},
                     w::AbstractVector{<:Real}; rtol=nothing,
                     return_metadata::Bool=false)
    rtol === nothing || throw(ArgumentError(
        "depth_rtol spectral truncation has been removed; depth now uses exact covariance rank and strict affine support. Omit depth_rtol."))
    length(u) == length(v) == length(w) ||
        throw(DimensionMismatch("depth coordinate/weight lengths differ"))
    all(isfinite, u) && all(isfinite, v) ||
        throw(ArgumentError("depth coordinates must be finite"))
    all(x -> isfinite(x) && x >= 0, w) ||
        throw(ArgumentError("depth weights must be finite and nonnegative"))
    support = findall(>(0), w)
    isempty(support) && throw(ArgumentError("depth requires positive reference mass"))

    R = Rational{BigInt}
    rw = R.(w)
    total = sum(rw)
    # Exact subtraction also handles opposite-signed floatmax values without
    # overflowing. The anchor does not introduce order dependence: all following
    # arithmetic is exact and the final quantities are translation invariant.
    anchor = first(support)
    x = R.(u) .- R(u[anchor])
    y = R.(v) .- R(v[anchor])
    mx = sum(rw[i] * x[i] for i in support) / total
    my = sum(rw[i] * y[i] for i in support) / total
    x .-= mx
    y .-= my
    a = sum(rw[i] * x[i]^2 for i in support) / total
    b = sum(rw[i] * x[i] * y[i] for i in support) / total
    c = sum(rw[i] * y[i]^2 for i in support) / total
    determinant = a*c - b*b
    determinant >= 0 || error("exact reference covariance is not positive semidefinite")
    rank = determinant > 0 ? 2 : (a > 0 || c > 0 ? 1 : 0)
    depth = zeros(Float64, length(u))
    if rank == 2
        for i in eachindex(depth)
            numerator = c*x[i]^2 - 2*b*x[i]*y[i] + a*y[i]^2
            depth[i] = _rational_float(determinant / (determinant + numerator))
        end
    elseif rank == 1
        for i in eachindex(depth)
            # A zero covariance eigenvalue is a strict affine-support rule,
            # not a tolerance tube. Use either nonconstant coordinate on the
            # line; these give the same intrinsic Mahalanobis distance.
            if a > 0
                a*y[i] == b*x[i] || continue
                depth[i] = _rational_float(a / (a + x[i]^2))
            else
                x[i] == 0 || continue
                depth[i] = _rational_float(c / (c + y[i]^2))
            end
        end
    else
        for i in eachindex(depth)
            depth[i] = x[i] == 0 && y[i] == 0 ? 1.0 : 0.0
        end
    end
    parameters = (raw_scale=:mahalanobis_depth, solver=:exact_rational_2x2,
                  effective_rank=rank, rank_policy=:exact_covariance_rank,
                  support_policy=:exact_affine_span, rank_truncation=false)
    return_metadata ? (depth, parameters) : depth
end
