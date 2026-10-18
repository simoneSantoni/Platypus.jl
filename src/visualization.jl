"""
    TypicalityPlot

Backend-independent plotting recipe with detached data. Render with
`using Plots; plot(plot_feature_distributions(es))`, or inspect `.series`.
Platypus depends only on RecipesBase; install a Plots backend separately.
"""
struct TypicalityPlot
    series::Vector{NamedTuple}
    layout::Tuple{Int,Int}
end

_plotseries(kind,args; label="",xlabel="",ylabel="",title="",subplot=1) =
    (kind=kind,args=deepcopy(args),label=label,xlabel=xlabel,ylabel=ylabel,title=title,subplot=subplot)
TypicalityPlot(series::Vector{<:NamedTuple}) = TypicalityPlot(NamedTuple[series...],(1,1))

@recipe function f(p::TypicalityPlot)
    layout --> p.layout
    legend --> false
    for s in p.series
        @series begin
            seriestype := s.kind
            subplot := s.subplot
            label --> s.label
            # Omit unspecified labels. Multi-series helpers also repeat shared
            # subplot labels to protect them from downstream Plots defaults.
            if !isempty(s.xlabel)
                xlabel --> s.xlabel
            end
            if !isempty(s.ylabel)
                ylabel --> s.ylabel
            end
            if !isempty(s.title)
                title --> s.title
            end
            if s.kind === :heatmap
                # With legend=false, Plots' automatic colorbar is hidden.
                # A matrix of scores still needs an explicit numeric color key.
                colorbar --> true
            end
            s.args
        end
    end
    nothing
end

"""`plot_feature_distributions(es)`: histogram recipe, one panel per raw feature."""
function plot_feature_distributions(es::EntitySet)
    _validate_entityset(es)
    series = [_plotseries(:histogram,(es.matrix[:,j],);xlabel=string(es.feature_names[j]),
              ylabel="Count",title=string(es.feature_names[j]),subplot=j) for j in 1:es.n_features]
    TypicalityPlot(NamedTuple[series...],(cld(es.n_features,2),min(es.n_features,2)))
end

function _pair_index(ps::PairSet,pair)
    if pair isa Integer
        1 <= pair <= ps.n_pairs || throw(BoundsError(ps.pairs,pair))
        return Int(pair)
    elseif pair isa Tuple && length(pair) == 2
        names = Symbol.(pair)
        index = findfirst(fp -> (fp.name_i,fp.name_j) == names,ps.pairs)
        index === nothing && throw(ArgumentError("unknown ordered feature pair $pair"))
        return index
    end
    throw(ArgumentError("pair must be an index or a tuple of feature names"))
end

"""`plot_pair_cloud(ps, pair)`: scatter recipe for a pair index or `(feature_a,feature_b)` tuple."""
function plot_pair_cloud(ps::PairSet,pair)
    _validate_pairset(ps)
    fp = ps.pairs[_pair_index(ps,pair)]
    TypicalityPlot([_plotseries(:scatter,(fp.u,fp.v);xlabel=string(fp.name_i),
                               ylabel=string(fp.name_j),title="Pair cloud ($(ps.transform))")])
end

"""`plot_correlation_heatmap(es)`: Spearman correlation heatmap; constant-feature correlations are undefined."""
function plot_correlation_heatmap(es::EntitySet)
    _validate_entityset(es)
    C = corspearman(es.matrix)
    names = string.(es.feature_names)
    TypicalityPlot([_plotseries(:heatmap,(names,names,C);title="Spearman correlations")])
end

"""`plot_pair_redundancy(ps)`: absolute within-pair correlations, a diagnostic only; no aggregation weights are applied."""
function plot_pair_redundancy(ps::PairSet)
    df = pair_redundancy(ps)
    names = string.(df.feature_a) .* " × " .* string.(df.feature_b)
    TypicalityPlot([_plotseries(:bar,(names,df.abs_correlation);ylabel="|correlation|",
                               title="Within-pair dependence (not aggregation weights)")])
end

"""
    plot_typicality_profile(r, entity)

Histogram of normalized pair scores. Integer `entity` means row index; use
`entity_id => id` written as `:entity_id => id` to select by identifier explicitly.
"""
function plot_typicality_profile(r::TypicalityResult,entity)
    i = if entity isa Integer
        1 <= entity <= r.n_entities || throw(BoundsError(r.entity_ids,entity))
        Int(entity)
    elseif entity isa Pair && first(entity) === :entity_id
        row = findfirst(==(last(entity)),r.entity_ids)
        row === nothing && throw(ArgumentError("unknown entity ID"))
        row
    else
        throw(ArgumentError("use a row index or :entity_id => id"))
    end
    TypicalityPlot([_plotseries(:histogram,(r.pair_scores[i,:],);
        xlabel="Pairwise typicality (higher = typical)",ylabel="Count",title="Entity $(r.entity_ids[i])")])
end

"""`plot_atypicality_distribution(r)`: histogram of headline descriptive scores."""
function plot_atypicality_distribution(r::Union{TypicalityResult,GraphTypicalityResult})
    TypicalityPlot([_plotseries(:histogram,(r.atypicality,);xlabel="Atypicality (higher = unusual)",ylabel="Count")])
end

"""`plot_novelty_conventionality(r)`: conventionality versus lower-tail typicality; lower vertical values are more unusual."""
function plot_novelty_conventionality(r::TypicalityResult)
    TypicalityPlot([_plotseries(:scatter,(r.conventionality,r.novelty);
        xlabel="Conventionality (higher = typical)",ylabel="Lower-tail typicality (lower = unusual)")])
end

"""`plot_outcome_fit(r::OutcomeResult)`: observed versus fitted outcomes with identity line."""
function plot_outcome_fit(r::OutcomeResult)
    lo = min(minimum(r.observed),minimum(r.fitted))
    hi = max(maximum(r.observed),maximum(r.fitted))
    TypicalityPlot([_plotseries(:scatter,(r.observed,r.fitted);xlabel="Observed",ylabel="Fitted",title="Descriptive outcome fit"),
                    # Plots applies subplot defaults while processing each
                    # series, so repeat these labels on the overlay as well.
                    _plotseries(:path,([lo,hi],[lo,hi]);label="Identity",
                        xlabel="Observed",ylabel="Fitted",title="Descriptive outcome fit")])
end

"""`plot_pair_heatmap(r)`: entity × feature-pair normalized typicality heatmap."""
function plot_pair_heatmap(r::TypicalityResult)
    names = [string(a," × ",b) for (a,b) in r.pair_names]
    TypicalityPlot([_plotseries(:heatmap,(names,string.(r.entity_ids),r.pair_scores);
        xlabel="Feature pair",ylabel="Entity ID",title="Pairwise typicality")])
end
