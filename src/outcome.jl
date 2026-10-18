"""
    OutcomeResult

Descriptive least-squares fit. Stores coefficients, term names, observed/fitted
outcomes, residuals, R², entity IDs, design matrix and measurement provenance.
No independent-observation standard errors, p-values, causal claims, or temporal
validation are implied. Dependent generated predictors need study-specific
resampling that refits the measurement pipeline.
"""
struct OutcomeResult
    coefficients::Vector{Float64}
    terms::Vector{Symbol}
    observed::Vector{Float64}
    fitted::Vector{Float64}
    residuals::Vector{Float64}
    r2::Float64
    entity_ids::Vector{Int}
    design::Matrix{Float64}
    configuration::NamedTuple
end

function Base.show(io::IO, r::OutcomeResult)
    print(io,"OutcomeResult($(length(r.observed)) observations, $(length(r.terms)) terms, descriptive OLS, R²=$(round(r.r2;digits=4)))")
end
analysis_metadata(r::OutcomeResult) = deepcopy(r.configuration)

"""
    fit_outcome(r::TypicalityResult, y; predictors=(:conventionality,:novelty,:dispersion),
                intercept=true, interaction=false)

Descriptive OLS fitted with column-pivoted QR, avoiding normal equations.
`y` must be finite, in result row order, and have one entry per entity. Supported
predictors are `:conventionality`, `:novelty`, `:dispersion`, `:atypicality`.
With `interaction=true`, adds conventionality × (1 - novelty), requiring both
base terms. Rank-deficient designs are rejected rather than silently selecting
one of many coefficient vectors. Including conventionality and atypicality with
an intercept is therefore invalid.

This function does not construct a historical reference, establish causality,
or provide uncertainty estimates. Explicitly select measurement features before
fitting typicality; never include the outcome in that feature schema.
"""
function fit_outcome(r::TypicalityResult, y::AbstractVector;
        predictors=(:conventionality,:novelty,:dispersion),
        intercept::Bool=true, interaction::Bool=false)
    length(y) == r.n_entities || throw(DimensionMismatch("outcome length differs from entities"))
    all(x -> x isa Real && isfinite(x),y) || throw(ArgumentError("outcomes must be finite real values"))
    observed = Float64.(y)
    all(isfinite,observed) || throw(ArgumentError("outcomes cannot be represented as finite Float64 values"))
    selected = Symbol.(collect(predictors))
    allunique(selected) || throw(ArgumentError("predictor terms must be unique"))
    allowed = (:conventionality,:novelty,:dispersion,:atypicality)
    all(in(allowed),selected) || throw(ArgumentError("unsupported predictor term"))
    terms = Symbol[]
    columns = Vector{Float64}[]
    if intercept
        push!(terms,:intercept)
        push!(columns,ones(r.n_entities))
    end
    for term in selected
        push!(terms,term)
        push!(columns,copy(getfield(r,term)))
    end
    if interaction
        :conventionality in selected && :novelty in selected ||
            throw(ArgumentError("interaction requires conventionality and novelty main effects"))
        push!(terms,:conventionality_x_tail_atypicality)
        push!(columns,r.conventionality .* (1 .- r.novelty))
    end
    isempty(columns) && throw(ArgumentError("at least one model term is required"))
    X = reduce(hcat,columns)
    # reduce(hcat, ...) may leave a vector for one input.
    X = reshape(X,r.n_entities,length(columns))
    all(isfinite,X) || throw(ArgumentError("result contains nonfinite predictor values"))
    size(X,1) > size(X,2) || throw(ArgumentError("more observations than coefficients are required"))
    rank(X) == size(X,2) || throw(ArgumentError("outcome design is rank deficient; remove redundant or constant terms"))
    scale = maximum(abs,observed)
    # Orthogonal projection can overflow on a finite outcome vector (for
    # example, an intercept-only response of 1e308). Solve in bounded outcome
    # units, and form residuals before returning to the original units.
    outcome_scale = scale == 0 ? 1.0 : scale
    scaled_observed = observed ./ outcome_scale
    scaled_beta = qr(X,ColumnNorm()) \ scaled_observed
    scaled_fitted = X * scaled_beta
    scaled_residuals = scaled_observed .- scaled_fitted
    beta = scaled_beta .* outcome_scale
    fitted = scaled_fitted .* outcome_scale
    residuals = scaled_residuals .* outcome_scale
    all(isfinite,beta) && all(isfinite,fitted) && all(isfinite,residuals) ||
        throw(ArgumentError("outcome coefficients, fitted values or residuals exceed finite Float64 range; rescale the outcome or predictors"))
    centered = scaled_observed .- mean(scaled_observed)
    total = sum(abs2,centered)
    r2 = total == 0 ? NaN : 1 - sum(abs2,scaled_residuals)/total
    config = (estimator=:descriptive_ols, solver=:pivoted_qr,
              measurement=analysis_metadata(r), predictors=copy(selected),
              intercept=intercept, interaction=interaction)
    OutcomeResult(beta,terms,observed,fitted,residuals,r2,copy(r.entity_ids),copy(X),config)
end
