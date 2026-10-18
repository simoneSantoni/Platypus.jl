"""
    TypicalityResult

Snapshot of pairwise descriptive scores. `pair_scores` contains reference-weighted
mid-percentiles, with larger values meaning more typical. `raw_scores` preserves
method-specific scores, also oriented toward typicality. `normalization_scores`
contains the values whose weighted mid-CDF defines `pair_scores`. These equal
`raw_scores` except for real-coordinate KDE, which uses ordinal labels for the
scale-free kernel ordering. These labels preserve ties and order; their
magnitudes are not density estimates and cannot be compared across fits or
pairs. Physical-unit adjustments can round distinct densities to the same
Float64 raw log score; such raw ties may have different percentiles when the
normalization ordering resolves the distinction. Rank-space KDE uses raw scores.
Summaries always use
`pair_scores`: median conventionality, lower-quantile `novelty` (smaller means more
unusual), IQR dispersion, and `1 - conventionality` atypicality. Neither scores nor
IQR are probabilities of anomaly or sampling uncertainty.

`configuration` records transformation, normalized reference weights, effective
sample size, algorithm settings, per-pair parameters, version, and caller-supplied
provenance. Metadata and numerical arrays are detached from the inputs.
"""
struct TypicalityResult
    pair_scores::Matrix{Float64}
    raw_scores::Matrix{Float64}
    normalization_scores::Matrix{Float64}
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
    configuration::NamedTuple
end

function Base.show(io::IO, r::TypicalityResult)
    print(io, "TypicalityResult($(r.n_entities) entities × $(r.n_pairs) pairs, ",
          "method=$(r.method), transform=$(r.configuration.transform)")
    r.k > 0 && print(io, ", k=$(r.k)")
    print(io, ", normalized mid-percentiles)")
end

"""`atypicality(r)`: detached headline scores, `1 - median(pair_scores)`. Higher is more atypical."""
atypicality(r::TypicalityResult) = copy(r.atypicality)
"""`conventionality(r)`: detached median pairwise typicality, higher means more typical."""
conventionality(r::TypicalityResult) = copy(r.conventionality)
"""`novelty(r)`: detached lower-tail typicality; **smaller** means more unusual. This is not a minimum or historical newness."""
novelty(r::TypicalityResult) = copy(r.novelty)
"""`dispersion(r)`: detached IQR across pair scores, not a standard error or confidence interval."""
dispersion(r::TypicalityResult) = copy(r.dispersion)
"""`analysis_metadata(r)`: detached analysis configuration to save with exported scores."""
analysis_metadata(r::TypicalityResult) = deepcopy(r.configuration)

"""`summarytable(r)`: one row per entity. Save `analysis_metadata(r)` alongside this table."""
function summarytable(r::TypicalityResult)
    DataFrame(entity_id=r.entity_ids, label=r.labels, conventionality=r.conventionality,
              novelty=r.novelty, dispersion=r.dispersion, atypicality=r.atypicality)
end

"""
    scoretable(r::TypicalityResult; raw=false, normalization=false)

Long-form fitted pairwise scores, normalized by default. With `raw=true`, the
`typicality` column contains method-specific `raw_scores`; their units and scale
are recorded in `analysis_metadata(r)` and must not be compared between methods.
With `normalization=true`, it contains `normalization_scores`, the values used
to order and tie observations for the weighted mid-CDF. For real-coordinate KDE
these are ordinal labels of the scale-free kernel ordering, not density values;
they can distinguish observations whose Float64 raw log densities coincide.
Their magnitudes are not comparable across pairs or fits. All other methods,
including rank-space KDE, use raw scores. The two options are mutually exclusive.
"""
function scoretable(r::TypicalityResult; raw::Bool=false, normalization::Bool=false)
    raw && normalization && throw(ArgumentError("raw and normalization exports are mutually exclusive"))
    values = raw ? r.raw_scores : normalization ? r.normalization_scores : r.pair_scores
    n, p = size(values)
    DataFrame(entity_id=repeat(r.entity_ids; inner=p),
              feature_a=repeat(first.(r.pair_names); outer=n),
              feature_b=repeat(last.(r.pair_names); outer=n),
              typicality=vec(permutedims(values)))
end

# Sorting makes a reduction reproducible for a fixed multiset; it does not
# establish mathematical equality of independently computed summands.
_stable_sum(values) = sum(sort!(collect(Float64, values)))

# Reference-weighted mid-CDF. Zero-weight entities are queries, never reference
# mass. Infinite log scores are valid boundaries; NaNs are always invalid.
function _score_percentiles(raw::AbstractVector, w::Vector{Float64})
    length(raw) == length(w) || throw(DimensionMismatch("score/weight lengths differ"))
    any(isnan, raw) && throw(ArgumentError("estimator produced NaN scores"))
    order = sortperm(raw)
    result = similar(w)
    exact_weights = Rational{BigInt}.(w)
    total = sum(exact_weights)
    cumulative = zero(total)
    i = 1
    while i <= length(order)
        j = i
        while j < length(order) && raw[order[j+1]] == raw[order[i]]
            j += 1
        end
        mass = sum(exact_weights[order[l]] for l in i:j)
        value = _rational_float((cumulative + mass/2)/total)
        for l in i:j
            result[order[l]] = value
        end
        cumulative += mass
        i = j + 1
    end
    result
end

function _fit_setup(ps, weights, novelty_quantile, normalize)
    _validate_pairset(ps)
    ps.n_pairs > 0 || throw(ArgumentError("at least one feature pair is required"))
    ps.n_entities >= 2 || throw(ArgumentError("at least two entities are required"))
    isfinite(novelty_quantile) && 0 <= novelty_quantile <= 1 ||
        throw(ArgumentError("novelty_quantile must lie in [0,1]"))
    normalize || throw(ArgumentError(
        "normalize=false no longer creates ambiguous raw summaries; use result.raw_scores or scoretable(result; raw=true)"))
    rw = weights === nothing ? ReferenceWeights(ps.n_entities) : weights
    w = _validated_weights(rw, ps.n_entities; min_positive=2)
    uniform = all(==(first(rw.weights)), rw.weights)
    ps.source_matrix === nothing && !uniform && throw(ArgumentError(
        "weighted fitting requires build_pairs provenance; manually constructed PairSet has unknown source transform"))
    columns = ps.source_matrix === nothing ? nothing :
        _transform_columns(ps.source_matrix, ps.transform, uniform ? nothing : rw.weights)
    return w, uniform, columns
end

function _summarize(ps, raw, method, k, q, w, configuration;
        score_weights=w, normalization_scores=raw, normalization_basis::Symbol=:raw_score)
    n, p = size(raw)
    size(normalization_scores) == size(raw) ||
        throw(DimensionMismatch("normalization and raw score dimensions differ"))
    scores = similar(raw)
    for j in 1:p
        scores[:,j] = _score_percentiles(@view(normalization_scores[:,j]), score_weights)
    end
    conv = [median(@view scores[i,:]) for i in 1:n]
    nov = [quantile(@view(scores[i,:]), q) for i in 1:n]
    disp = [p > 1 ? iqr(@view(scores[i,:])) : 0.0 for i in 1:n]
    names = [(fp.name_i, fp.name_j) for fp in ps.pairs]
    config = merge((package_version="0.2.0", transform=ps.transform,
                    feature_names=copy(ps.feature_names),
                    entity_ids=copy(ps.entity_ids),
                    reference_weights=copy(w), effective_sample_size=1/_stable_sum(abs2.(w)),
                    source_reference_weights=copy(score_weights),
                    normalization=:weighted_midcdf, normalization_basis=normalization_basis,
                    method=method, k=k,
                    novelty_quantile=Float64(q)), configuration)
    TypicalityResult(scores, raw, copy(normalization_scores), conv, nov, disp, 1 .- conv, names,
        copy(ps.entity_ids), copy(ps.labels), method, k, Float64(q), true, n, p, deepcopy(config))
end

"""
    fit_typicality(ps::PairSet; method=:lof, k=0, novelty_quantile=0.1,
                   weights=nothing, normalize=true, bandwidth=nothing,
                backend=:kdtree, parallel=false, depth_rtol=nothing,
                   provenance=(;))

Fit descriptive pairwise typicality relative to the supplied collection.

- `:knn`: negative kth-neighbor distance (larger is more typical), using exact
  search. `:lof`: reciprocal original, tie-inclusive Local Outlier Factor.
- `:kde`: leave-one-out product log density; reflected Epanechnikov kernels on rank
  coordinates, Gaussian kernels on unbounded coordinates. `bandwidth` is a
  positive scalar or a pair of bandwidths; automatic bandwidth uses weighted
  effective sample size and the bivariate rule. No marginal-uniformity guarantee
  is claimed for this density estimate.
- `:depth`: weighted global Mahalanobis centrality, including the focal entity;
  not a local-density measure. Exact rational covariance determines algebraic
  rank and strict affine support. Off-support depth is zero. Spectral truncation
  has been removed; omit the retired `depth_rtol` argument. Rank geometry is
  retained exactly, and raw/z-score depth uses original coordinate provenance.

`k=0` selects rounded sqrt(N); explicit k must lie in `1:N-1`. Nonuniform
weights are supported by KDE/depth and applied to marginal transformations,
reference fitting and score normalization. kNN/LOF reject nonuniform weights.
Zero-weight rows are scored as queries but do not affect weighted references.
Each pair is normalized to reference-weighted mid-percentiles before aggregation.
All-equal normalization values map to 0.5. Raw method scores remain in
`result.raw_scores`. Real-coordinate KDE normalization uses ordinal labels of
the scale-free kernel ordering stored in `result.normalization_scores`. These
labels are not density estimates. Raw Float64 log densities can therefore tie
while their percentiles differ;
export the ordering values with `scoretable(result; normalization=true)`.
Other methods, including rank-space KDE, use raw scores for normalization.

`parallel=true` independently scores pairs with Julia threads. It does not use
approximate neighbors. Both requested neighbor backends currently execute a
filtered exact linear scan, with quadratic work over all entities; configuration
records requested and effective backends. `provenance` is an optional NamedTuple for dataset version,
cohort selection, time cutoff and preprocessing notes; it is copied into results.
"""
function fit_typicality(ps::PairSet; method::Symbol=:lof, k::Integer=0,
        novelty_quantile::Real=0.1, normalize::Bool=true,
        weights::Union{ReferenceWeights,Nothing}=nothing, bandwidth=nothing,
        backend::Symbol=:kdtree, parallel::Bool=false, depth_rtol=nothing,
        provenance::NamedTuple=(;))
    method in TYPICALITY_METHODS || throw(ArgumentError("unknown method :$method; expected $TYPICALITY_METHODS"))
    depth_rtol === nothing || throw(ArgumentError("depth_rtol spectral truncation has been removed; omit it to use exact covariance rank and strict affine support"))
    backend in (:kdtree,:brute) || throw(ArgumentError("backend must be :kdtree or :brute"))
    w, uniform, columns = _fit_setup(ps,weights,novelty_quantile,normalize)
    source_weights=weights===nothing ? ones(ps.n_entities) : weights.weights
    neighbor = method in (:knn,:lof)
    rank_lattice = neighbor && ps.transform===:rank && ps.source_matrix!==nothing ?
        _neighbor_rank_coordinates(ps.source_matrix) : nothing
    if ps.source_matrix !== nothing && (ps.transform===:rank || method===:depth)
        columns = _exact_transform_columns(ps.source_matrix,ps.transform===:rank ? :rank : :none,source_weights)
    end
    neighbor && !uniform && throw(ArgumentError("$method does not support nonuniform reference weights; use :kde or :depth"))
    k >= 0 || throw(ArgumentError("k must be non-negative; zero selects the default"))
    !neighbor && k != 0 && throw(ArgumentError("k is only meaningful for :knn/:lof"))
    method !== :kde && bandwidth !== nothing && throw(ArgumentError("bandwidth is only meaningful for :kde"))
    kk = neighbor ? (k == 0 ? clamp(round(Int,sqrt(ps.n_entities)),1,ps.n_entities-1) : Int(k)) : 0
    neighbor && !(1 <= kk < ps.n_entities) && throw(ArgumentError("k must lie in 1:N-1"))
    raw = Matrix{Float64}(undef,ps.n_entities,ps.n_pairs)
    ordering = similar(raw)
    parameters = Vector{NamedTuple}(undef,ps.n_pairs)
    function score_pair(p)
        fp = ps.pairs[p]
        u,v = columns === nothing ? (fp.u,fp.v) : (columns[fp.i],columns[fp.j])
        if rank_lattice !== nothing
            u,v=collect(view(rank_lattice,fp.i,:)),collect(view(rank_lattice,fp.j,:))
        end
        values, pars, keys = if method === :kde
            _density_scores(u,v,source_weights;
                support=ps.transform === :rank ? :unit : :real,
                bandwidth=bandwidth, logscale=true, return_ordering=true)
        else
            values, pars = if method === :knn
                distance_scale=rank_lattice===nothing ? 1 : 2ps.n_entities
                (-_neighbor_kdistance(u,v,kk;backend=backend)./distance_scale, (raw_scale=:negative_distance,))
            elseif method === :lof
                (1 ./ _neighbor_lof(u,v,kk;backend=backend), (raw_scale=:inverse_lof,))
            else
                _maha_depth(u,v,source_weights;return_metadata=true)
            end
            (values, pars, values)
        end
        all(x -> isfinite(x) || (method === :kde && x == -Inf),values) ||
            throw(ArgumentError("$method produced invalid raw scores in pair $(fp.name_i), $(fp.name_j)"))
        raw[:,p] = values
        all(x -> isfinite(x) || (method === :kde && x == -Inf), keys) ||
            throw(ArgumentError("$method produced invalid normalization scores in pair $(fp.name_i), $(fp.name_j)"))
        ordering[:,p] = keys
        parameters[p] = pars
    end
    if parallel
        Threads.@threads for p in 1:ps.n_pairs
            score_pair(p)
        end
    else
        for p in 1:ps.n_pairs
            score_pair(p)
        end
    end
    _summarize(ps,raw,method,kk,novelty_quantile,w,
        (backend=neighbor ? :brute : :direct, requested_backend=neighbor ? backend : :direct,
         effective_backend=neighbor ? :brute : :direct,
         distance_policy=neighbor ? :filtered_exact : :not_applicable,
         geometry=ps.transform===:rank ? :exact_reference_midrank :
             method===:depth ? :source_binary : :transformed_binary, parallel=parallel,
         pair_parameters=parameters, provenance=provenance);
        score_weights=source_weights, normalization_scores=ordering,
        normalization_basis=method === :kde && ps.transform !== :rank ?
            :scale_free_kernel_order : :raw_score)
end

"""
    fit_copula(ps::PairSet; family=:gaussian, weights=nothing,
               regularization=0.0, novelty_quantile=0.1, normalize=true,
               query_clip=nothing, provenance=(;))

Score genuine fitted bivariate copula densities on `transform=:rank` pairs.
`:gaussian` fits a normal-score correlation; singular fits require explicit
correlation shrinkage `regularization` in (0,1). `:beta` uses the empirical beta
copula of Segers, Sibuya and Tsukahara, requiring untied continuous coordinates
and uniform positive weights. Beta fitting is in-sample, not leave-one-out.
These model assumptions differ from LOF and nonparametric KDE. Raw scores are
log copula densities; normalized summaries have the same contract as
[`fit_typicality`](@ref). This is not a null-hypothesis test.

Gaussian endpoint queries are rejected unless `query_clip` is supplied in
`(0,0.5)`; it clips only zero-weight endpoint queries, with the choice recorded
in metadata. Positive-reference endpoints always raise an error. The beta
family does not support clipping or regularization.
"""
function fit_copula(ps::PairSet; family::Symbol=:gaussian,
        weights::Union{ReferenceWeights,Nothing}=nothing, regularization::Real=0.0,
        novelty_quantile::Real=0.1, normalize::Bool=true, query_clip=nothing,
        provenance::NamedTuple=(;))
    ps.transform === :rank || throw(ArgumentError("fit_copula requires build_pairs(...; transform=:rank)"))
    family in (:gaussian,:beta) || throw(ArgumentError("family must be :gaussian or :beta"))
    w, _, columns = _fit_setup(ps,weights,novelty_quantile,normalize)
    source_weights=weights===nothing ? ones(ps.n_entities) : weights.weights
    columns=_exact_transform_columns(ps.source_matrix,:rank,source_weights)
    raw = Matrix{Float64}(undef,ps.n_entities,ps.n_pairs)
    parameters = Vector{NamedTuple}(undef,ps.n_pairs)
    for (p,fp) in enumerate(ps.pairs)
        u,v = columns[fp.i],columns[fp.j]
        raw[:,p],parameters[p] = _copula_scores(u,v,source_weights;family=family,regularization=regularization,query_clip=query_clip)
        any(isnan,@view(raw[:,p])) && throw(ArgumentError("copula produced NaN scores"))
    end
    _summarize(ps,raw,Symbol(family,:_copula),0,novelty_quantile,w,
        (family=family, regularization=Float64(regularization), query_clip=query_clip,
         pair_parameters=parameters, provenance=provenance);score_weights=source_weights)
end
