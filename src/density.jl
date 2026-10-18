"""
Weighted density helpers. Bounded rank-space KDE uses an exactly normalized
reflected Epanechnikov kernel; real-space KDE uses Gaussian kernels. Copula
helpers return log densities. See `docs/research/statistics.md` for derivations,
adoption decisions, and limitations.
"""

function _density_reference_weights(w::AbstractVector{<:Real}, n::Int)
    length(w) == n || throw(DimensionMismatch("one weight is required per entity"))
    all(x -> isfinite(x) && x >= 0, w) ||
        throw(ArgumentError("density reference weights must be finite and nonnegative"))
    maximum(w; init=0.0) > 0 || throw(ArgumentError("positive reference mass is required"))
    exact = Rational{BigInt}.(w)
    p = _density_float.(exact ./ sum(exact))
    count(>(0), p) >= 2 ||
        throw(ArgumentError("at least two positive-weight reference entities are required"))
    return p
end

# Summing the same multiset in a fixed order makes row permutations reproducible.
_density_sum(x) = sum(sort!(collect(x)))

function _density_inputs(u, v, w)
    n = length(u)
    length(v) == n || throw(DimensionMismatch("pair coordinates must have equal lengths"))
    all(isfinite, u) && all(isfinite, v) ||
        throw(ArgumentError("density coordinates must be finite"))
    p = _density_reference_weights(w, n)
    return p, findall(>(0), p)
end

# Exact rational moments operate on the supplied, represented inputs. They avoid
# cancellation without inventing an epsilon tie policy, and preserve zero
# covariance on Cartesian product references. No global BigFloat precision is set.
_density_exact(x::Real) = Rational{BigInt}(x)
_density_float(x::Real) = Float64(x)
_density_float(x::Rational{BigInt}) = _rational_float(x)
# Shared rational square root uses exact rounding-boundary comparisons.
_density_sqrt_rational(q::Rational{BigInt}) = _rational_sqrt(q)

function _density_exact_weights(w)
    q = _density_exact.(w)
    return q ./ sum(q)
end
function _density_exact_variance(x, q, active)
    origin = _density_exact(x[first(active)])
    centered = [_density_exact(x[i]) - origin for i in active]
    mu = sum(q[i] * centered[j] for (j,i) in enumerate(active))
    return sum(q[i] * (centered[j]-mu)^2 for (j,i) in enumerate(active))
end

# Three images suffice exactly because the kernel support is [-1,1], h<=1.
_epanechnikov(z::Real) = abs(z) < 1 ? 0.75 * (1 - z * z) : 0.0
function _reflected_kernel(x::Real, center::Real, h::Real)
    0 <= x <= 1 && 0 <= center <= 1 && isfinite(h) && 0 < h <= 1 ||
        throw(ArgumentError("reflected kernel needs coordinates in [0,1] and bandwidth in (0,1]"))
    return exp(_reflected_logkernel(x, center, h))
end

function _reflected_logkernel(x, center, h)
    # Reflection exchanges the two image distances. Form them exactly before
    # conversion, and reduce them in a canonical order.
    a, b = _density_exact(x), _density_exact(center)
    distances = sort!([abs(a-b), a+b, 2-a-b])
    total = sum(_epanechnikov(_density_float(d)/h) for d in distances)
    return log(total) - log(h)
end

function _density_bandwidth(x, p, active, support, supplied; exact_weights=_density_exact_weights(p))
    constant = all(i -> x[i] == x[first(active)], active)
    if supplied !== nothing
        supplied isa Real && isfinite(supplied) && supplied > 0 ||
            throw(ArgumentError("bandwidth must be finite and positive"))
        supplied = _density_float(supplied)
        isfinite(supplied) && supplied > 0 ||
            throw(ArgumentError("bandwidth must be representable as a finite positive Float64 value"))
        support === :unit && supplied > 1 &&
            throw(ArgumentError("bounded Epanechnikov bandwidth must be at most one"))
    end
    constant && return (true, 1.0)
    supplied !== nothing && return (false, _density_float(supplied))
    # Center represented inputs exactly before taking the square root. This
    # preserves small differences at large offsets and coordinate reflection.
    variance = _density_exact_variance(x, exact_weights, active)
    sigma = _density_sqrt_rational(variance)
    neff = _density_float(inv(sum(abs2, exact_weights)))
    # sqrt(5) matches the Epanechnikov variance to the Gaussian's variance;
    # this is an explicit smoothing heuristic, not an optimal bandwidth claim.
    h = sigma * neff^(-1/6) * (support === :unit ? sqrt(5.0) : 1.0)
    support === :unit && (h = min(h, 1.0))
    isfinite(h) && h > 0 || throw(ArgumentError(
        "automatic bandwidth is numerically degenerate; supply a positive bandwidth"))
    return (false, h)
end

# Ordinary subtraction preserves accuracy for nearby values. Only its overflow
# needs an alternative order: opposite-sign inputs can have a representable
# standardized difference even when their unscaled difference is infinite.
function _density_standardized_difference(x, center, bandwidth)
    difference = x - center
    difference isa Rational && return _density_float(difference / _density_exact(bandwidth))
    return isfinite(difference) ? difference / bandwidth :
        x / bandwidth - center / bandwidth
end

"""
    _density_scores(u, v, w; support=:unit, bandwidth=nothing, logscale=false)

Return `(scores, parameters)` for kernel-contribution leave-one-out scoring.
Ranks and bandwidths are held fixed
at their full-reference values. `bandwidth` is an absolute scalar or a tuple
of two absolute coordinate bandwidths. Unit support uses reflected compact
kernels; real support uses Gaussian kernels. A constant reference coordinate
is dropped; a query differing on that coordinate receives density zero.
Two constant coordinates give density one on their joint reference support.
"""
function _density_scores(u::AbstractVector{<:Real}, v::AbstractVector{<:Real},
                         w::AbstractVector{<:Real}; support::Symbol=:unit,
                         bandwidth=nothing, logscale::Bool=false,
                         return_ordering::Bool=false)
    support in (:unit, :real) || throw(ArgumentError("support must be :unit or :real"))
    p, active = _density_inputs(u, v, w)
    if support === :unit
        all(x -> 0 <= x <= 1, u) && all(x -> 0 <= x <= 1, v) ||
            throw(ArgumentError("rank-space KDE coordinates must lie in [0,1]"))
    end
    bw = if bandwidth === nothing || bandwidth isa Real
        (bandwidth, bandwidth)
    elseif bandwidth isa Tuple && length(bandwidth) == 2
        bandwidth
    else
        throw(ArgumentError("bandwidth must be a positive scalar or a two-element tuple"))
    end
    exact_weights = _density_exact_weights(w)
    cu, hu = _density_bandwidth(u, p, active, support, bw[1]; exact_weights)
    cv, hv = _density_bandwidth(v, p, active, support, bw[2]; exact_weights)
    result = Vector{Float64}(undef, length(u))
    terms = Vector{Float64}(undef, length(active))
    exponents = similar(terms)
    # Kernel masses retain exact sums of represented Gaussian exponential
    # terms. Keeping these separate from dimensional log normalizers avoids
    # destroying a real ordering when a unit change adds a large constant.
    masses = zeros(Rational{BigInt}, length(u))
    logmasses = fill(-Inf, length(u))
    mass_offsets = zeros(Float64,length(u))
    normconstant = 0.5 * log(2pi)
    lognormalizer = (cu ? 0.0 : log(hu)+normconstant) +
                    (cv ? 0.0 : log(hv)+normconstant)
    for i in eachindex(u)
        remaining = _density_sum(p[j] for j in active if j != i)
        remaining > 0 || throw(ArgumentError("entity $i has no leave-one-out reference mass"))
        if (cu && u[i] != u[first(active)]) || (cv && v[i] != v[first(active)])
            result[i] = -Inf
            continue
        end
        if cu && cv
            masses[i] = 1
            logmasses[i] = 0.0
            result[i] = 0.0
            continue
        end
        kernel_underflow = false
        for (a, j) in enumerate(active)
            if j == i
                terms[a] = -Inf
                exponents[a] = -Inf
                continue
            end
            if support === :unit
                lu = cu ? 0.0 : _reflected_logkernel(u[i], u[j], hu)
                lv = cv ? 0.0 : _reflected_logkernel(v[i], v[j], hv)
                terms[a] = log(p[j]) + (lu + lv)
            else
                du = cu ? 0.0 : _density_standardized_difference(u[i],u[j],hu)
                dv = cv ? 0.0 : _density_standardized_difference(v[i],v[j],hv)
                exponent = -0.5 * (du*du + dv*dv)
                exponents[a] = exponent
                terms[a] = log(p[j]) + exponent
                if return_ordering
                    kernel = exp(exponent)
                    kernel_underflow |= iszero(kernel)
                    masses[i] += exact_weights[j] * _density_exact(kernel)
                end
            end
        end
        logmasses[i] = _density_logsumexp(terms) - log(remaining)
        if support === :real
            if return_ordering
                remaining_exact = 1-exact_weights[i]
                if kernel_underflow
                    # Shift kernel exponents, keeping exact reference weights
                    # outside exp/log. This preserves tiny weight contrasts
                    # near unit conditional mass as well as single-reference ties.
                    offset = maximum(exponents)
                    if isfinite(offset)
                        mass_offsets[i] = offset
                        masses[i] = sum(exact_weights[j] *
                            _density_shifted_exp(exponents[a]-offset)
                            for (a,j) in enumerate(active) if j!=i) / remaining_exact
                    else
                        masses[i] = 0
                    end
                else
                    masses[i] /= remaining_exact
                end
            end
            result[i] = logmasses[i] - lognormalizer
        else
            result[i] = logmasses[i]
        end
    end
    scores = logscale ? result : exp.(result)
    !logscale && any(!isfinite, scores) && throw(ArgumentError(
        "raw KDE density overflows Float64; use logscale=true"))
    parameters = (support=support, bandwidth=(hu,hv), constant=(cu,cv),
                  effective_n=_density_float(inv(sum(abs2, exact_weights))), logscale=logscale,
                  raw_scale=logscale ? :log_density : :density,
                  normalization_basis=support===:real ? :scale_free_kernel_order : :raw_score,
                  loo=:kernel_contributions)
    if return_ordering
        ordering = support===:real ? _density_kernel_order(masses,mass_offsets) : copy(result)
        return scores, parameters, ordering
    end
    return scores, parameters
end

# Extend the represented exponential range without a floating product. After
# shifting, z<=0. Squaring a rational half-exponential retains contributions
# below Float64's minimum; this matters when a nearer kernel has tiny weight.
function _density_shifted_exp(z::Float64)
    value = exp(z)
    if iszero(value) && isfinite(z)
        return _density_exact(exp(z/2))^2
    end
    return _density_exact(value)
end

# Outward bounds constructed with explicit destination precision and rounding.
# Unlike generic BigFloat arithmetic these operations never read or alter the
# caller's global precision/rounding context (including Julia 1.9).
function _density_rational_bounds(q::Rational{BigInt}, bits::Int)
    down, up = Base.MPFR.MPFRRoundDown, Base.MPFR.MPFRRoundUp
    nlo = BigFloat(numerator(q),down;precision=bits)
    nhi = BigFloat(numerator(q),up;precision=bits)
    dlo = BigFloat(denominator(q),down;precision=bits)
    dhi = BigFloat(denominator(q),up;precision=bits)
    lower = BigFloat(0;precision=bits)
    upper = BigFloat(0;precision=bits)
    ld, ud = q >= 0 ? (dhi,dlo) : (dlo,dhi)
    ccall((:mpfr_div,Base.MPFR.libmpfr),Cint,
        (Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),
        lower,nlo,ld,down)
    ccall((:mpfr_div,Base.MPFR.libmpfr),Cint,
        (Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),
        upper,nhi,ud,up)
    return lower,upper
end

# Sign of log(q)-target. Apart from q=1,target=0 (handled exactly), equality
# cannot hold for a positive rational q and nonzero rational target. Refine
# enclosing intervals until their separation proves the ordering. Only logs
# are evaluated: huge negative exponents never become huge rational integers.
function _density_logratio_compare(q::Rational{BigInt}, target::Rational{BigInt})
    q > 0 || throw(DomainError(q,"positive kernel mass ratio required"))
    q == 1 && return target > 0 ? -1 : target < 0 ? 1 : 0
    iszero(target) && return q < 1 ? -1 : 1
    bits = 128
    down, up = Base.MPFR.MPFRRoundDown, Base.MPFR.MPFRRoundUp
    while true
        qlo,qhi = _density_rational_bounds(q,bits)
        tlo,thi = _density_rational_bounds(target,bits)
        lower = BigFloat(0;precision=bits)
        upper = BigFloat(0;precision=bits)
        ccall((:mpfr_log,Base.MPFR.libmpfr),Cint,
            (Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),lower,qlo,down)
        ccall((:mpfr_log,Base.MPFR.libmpfr),Cint,
            (Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),upper,qhi,up)
        upper < tlo && return -1
        lower > thi && return 1
        bits = Base.checked_mul(bits,2)
    end
end

function _density_kernel_compare(a,oa,b,ob)
    iszero(a) && return iszero(b) ? 0 : -1
    iszero(b) && return 1
    oa == ob && return a < b ? -1 : a > b ? 1 : 0
    return _density_logratio_compare(a/b,_density_exact(ob)-_density_exact(oa))
end

# Dense ordinal labels express only order/equality, not density magnitudes.
# Every row represents the positive real exp(offset)*mass. Normal rows retain
# exact direct sums (offset zero); underflow-affected rows retain shifted sums.
# Comparing these same real values supplies one strict weak order across both
# representations, without imposing an artificial underflow bucket.
function _density_kernel_order(masses, offsets)
    compare(i,j) = _density_kernel_compare(masses[i],offsets[i],masses[j],offsets[j])
    permutation = sortperm(eachindex(masses);lt=(i,j)->compare(i,j)<0)
    ordering = zeros(Float64,length(masses))
    level = 0.0
    previous = 0
    for i in permutation
        if previous == 0 || compare(i,previous) != 0
            level += 1
        end
        ordering[i] = level
        previous = i
    end
    return ordering
end

function _density_logsumexp(x)
    m = maximum(x)
    m == -Inf && return -Inf
    return m + log(_density_sum(exp(t-m) for t in x))
end

function _gaussian_copula_logdensity(z1::Real, z2::Real, rho::Real)
    isfinite(rho) && abs(rho) < 1 ||
        throw(ArgumentError("Gaussian copula correlation must lie strictly between -1 and 1"))
    rho == 0 && return 0.0
    return -0.5 * log1p(-rho*rho) -
        (rho*rho*(z1*z1+z2*z2)-2rho*(z1*z2)) / (2*(1-rho*rho))
end

"""
    _copula_scores(u,v,w; family=:gaussian, regularization=0.0, query_clip=nothing)

Return `(log_scores, parameters)`. Gaussian correlation is estimated from
weighted normal scores, not by maximum likelihood. `regularization` explicitly
shrinks correlation by `1-regularization`. Endpoint queries are rejected by
default; explicit `query_clip` in `(0,0.5)` clips zero-weight endpoint queries
only. Positive-reference endpoints always raise an error. Constant margins have no fitted
dependence and are rejected. The beta family requires equal positive weights
and untied marginal ranks; its density is the full in-sample empirical beta
copula, not a leave-one-out approximation.
"""
function _copula_scores(u::AbstractVector{<:Real}, v::AbstractVector{<:Real},
                        w::AbstractVector{<:Real}; family::Symbol=:gaussian,
                        regularization::Real=0.0, query_clip=nothing)
    family in (:gaussian, :beta) || throw(ArgumentError("copula family must be :gaussian or :beta"))
    isfinite(regularization) && 0 <= regularization < 1 ||
        throw(ArgumentError("regularization must lie in [0,1)"))
    query_clip === nothing || (query_clip isa Real && isfinite(query_clip) &&
        0 < query_clip < 0.5) || throw(ArgumentError("query_clip must lie in (0,0.5)"))
    p, active = _density_inputs(u, v, w)
    all(x -> 0 <= x <= 1, u) && all(x -> 0 <= x <= 1, v) ||
        throw(ArgumentError("copula coordinates must be rank coordinates in [0,1]"))
    if family === :beta
        regularization == 0 || throw(ArgumentError("regularization applies only to the Gaussian copula"))
        query_clip === nothing || throw(ArgumentError("query_clip applies only to the Gaussian copula"))
        length(active) == length(w) && all(==(first(w)), w) ||
            throw(ArgumentError("empirical beta copula requires equal strictly positive weights"))
        length(unique(u)) == length(u) && length(unique(v)) == length(v) ||
            throw(ArgumentError("empirical beta copula requires untied marginal ranks"))
        n = length(u)
        ru, rv = invperm(sortperm(u)), invperm(sortperm(v))
        function beta_logkernel(rank, x)
            # b(n,r,x) == b(n,n+1-r,1-x), evaluated on the same
            # canonical half interval even when the original rank is rational.
            x = _density_exact(x)
            if x > 1//2
                rank, x = n+1-rank, 1-x
            end
            return Distributions.logpdf(Distributions.Beta(rank,n+1-rank),_density_float(x))
        end
        scores = [_density_logsumexp([beta_logkernel(ru[j],u[i]) +
                   beta_logkernel(rv[j],v[i]) for j in 1:n]) - log(n) for i in 1:n]
        return scores, (family=:beta, n=n, in_sample=true)
    end
    normal = Distributions.Normal()
    function normal_score(x, i)
        if x == 0 || x == 1
            p[i] == 0 || throw(ArgumentError("positive-reference copula ranks must be strictly inside (0,1)"))
            query_clip === nothing && throw(ArgumentError(
                "endpoint copula query needs explicit query_clip or an interior reference rank"))
            x = clamp(x, query_clip, 1-query_clip)
        end
        # One quantile evaluation for complementary exact ranks preserves the
        # odd normal-score identity, rather than asking two approximations to agree.
        x = _density_exact(x)
        return x <= 1//2 ? Distributions.quantile(normal, _density_float(x)) :
            -Distributions.quantile(normal, _density_float(1-x))
    end
    zu = [normal_score(x,i) for (i,x) in enumerate(u)]
    zv = [normal_score(x,i) for (i,x) in enumerate(v)]
    all(isfinite, zu) && all(isfinite, zv) || throw(ArgumentError(
        "reference weights yield unrepresentable normal-score endpoints"))
    q = _density_exact_weights(w)
    au, av = _density_exact.(zu), _density_exact.(zv)
    mu = sum(q[i]*au[i] for i in active)
    mv = sum(q[i]*av[i] for i in active)
    vu = sum(q[i]*(au[i]-mu)^2 for i in active)
    vv = sum(q[i]*(av[i]-mv)^2 for i in active)
    vu > 0 && vv > 0 || throw(ArgumentError("Gaussian copula needs two nonconstant reference margins"))
    covariance = sum(q[i]*(au[i]-mu)*(av[i]-mv) for i in active)
    rho = clamp(_density_sqrt_rational(covariance^2/(vu*vv)), 0.0, 1.0)
    covariance < 0 && (rho = -rho)
    if regularization == 0 && 1-abs(rho) <= 8eps(Float64)
        throw(ArgumentError("singular Gaussian copula fit; choose positive regularization explicitly"))
    end
    rho *= 1-_density_float(regularization)
    abs(rho) < 1 || throw(ArgumentError("regularization is too small to resolve singular dependence"))
    scores = [_gaussian_copula_logdensity(zu[i],zv[i],rho) for i in eachindex(u)]
    return scores, (family=:gaussian, rho=rho, regularization=_density_float(regularization),
                    query_clip=query_clip, in_sample=true)
end
