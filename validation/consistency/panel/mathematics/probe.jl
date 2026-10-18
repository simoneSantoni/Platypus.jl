# Review probe: records current behavior; does not change estimators or audit results.
# julia --project=. --startup-file=no validation/consistency/panel/mathematics/probe.jl
using Platypus, LinearAlgebra, CSV, DataFrames

function fit_depth(X; weights=nothing, transform=:none, rtol=1e-12)
    fit_typicality(build_pairs(EntitySet(X, [:x, :y]); transform);
        method=:depth, depth_rtol=rtol,
        weights=weights === nothing ? nothing : ReferenceWeights(weights))
end

function exact_depth(X, weights=ones(Int, size(X,1)))
    R = Rational{BigInt}.(X)
    w = Rational{BigInt}.(weights); w ./= sum(w)
    Z = R .- reshape(vec(w' * R), 1, :)
    C = Z' * Diagonal(w) * Z
    [1 / (1 + dot(Z[i,:], C \ Z[i,:])) for i in axes(Z,1)]
end

function exact_percentiles(raw, weights=ones(Int,length(raw)))
    w = Rational{BigInt}.(weights); w ./= sum(w)
    [sum(w[j] * (raw[j] < raw[i] ? 1 : raw[j] == raw[i] ? 1//2 : 0)
        for j in eachindex(raw)) for i in eachindex(raw)]
end

function main(outdir=joinpath(@__DIR__, "results"))
    mkpath(outdir)
    observations = NamedTuple[]
    inputs = NamedTuple[]
    function record(example, X, actual, expected; weights=ones(size(X,1)))
        scores = exact_percentiles(expected, weights)
        for i in axes(X,1)
            push!(inputs, (; example, entity_id=i, x=X[i,1], y=X[i,2], weight=weights[i]))
            push!(observations, (; example, entity_id=i,
                raw_expected=Float64(expected[i]), raw_observed=actual.raw_scores[i,1],
                percentile_expected=Float64(scores[i]), percentile_observed=actual.pair_scores[i,1],
                raw_error=abs(Float64(expected[i])-actual.raw_scores[i,1]),
                percentile_error=abs(Float64(scores[i])-actual.pair_scores[i,1])))
        end
    end

    triangle = [0. 0.; 1. 0.; 0. 1.]
    oracle = exact_depth(triangle)
    @assert all(==(1//3), oracle)
    for transform in (:none, :rank, :zscore)
        record("triangle_$(transform)", triangle, fit_depth(triangle; transform), oracle)
    end
    weights = [1,2,3]
    weighted_oracle = exact_depth(triangle, weights)
    @assert weighted_oracle == Rational{BigInt}.(weights) ./ sum(weights)
    record("weighted_triangle_control", triangle, fit_depth(triangle; weights), weighted_oracle; weights)

    offset = [1.421875 -0.1875; 0.5625 2.078125; 0.421875 -1.65625;
              1.625 0.3125; 0.171875 -0.71875]
    shifted = offset .+ 1e14
    @assert shifted .- 1e14 == offset
    oracle = exact_depth(offset)
    @assert exact_depth(shifted) == oracle
    record("offset_centered", offset, fit_depth(offset), oracle)
    record("offset_shifted", shifted, fit_depth(shifted), oracle)
    anchored = shifted .- shifted[1:1,:]
    @assert exact_depth(anchored) == oracle
    record("offset_anchor_probe", anchored, fit_depth(anchored), oracle)

    # A representable perturbation creates two genuine score groups separated
    # by about 1e-13. A blanket 1e-12 score tie tolerance would erase them.
    near = [-1. -1.; 1. -1.; 1. 1.; -1. 1. + 2.0^-40]
    oracle = exact_depth(near)
    @assert length(unique(oracle)) == 2
    @assert 0 < maximum(oracle)-minimum(oracle) < 1//10^12
    @assert exact_percentiles(oracle) == [1//4, 3//4, 3//4, 1//4]
    record("genuine_near_ties_control", near, fit_depth(near), oracle)

    # Positive reference lies exactly on a line. Query is exactly off it.
    # The nonsingular integer shear maps the support to a constant coordinate.
    weights = [1,1,1,0]
    shear = [1. -1.; 0. 1.]
    @assert det(shear) == 1
    for delta in (2.0^-30, 2.0^-24, 2.0^-20, 2.0^-16)
        line = [-1. -1.; 0. 0.; 1. 1.; 0. delta]
        axis = line * shear
        @assert Rational{BigInt}.(line) * Rational{BigInt}.(shear) == Rational{BigInt}.(axis)
        @assert all(line[1:3,1] .== line[1:3,2]) && line[4,1] != line[4,2]
        @assert all(axis[1:3,2] .== 0) && axis[4,2] != 0
        # On-support squared distances are (3/2,0,3/2); off-support depth is 0.
        oracle = Rational{BigInt}[2//5, 1, 2//5, 0]
        record("line_query_$(delta)", line, fit_depth(line; weights), oracle; weights)
        record("axis_query_$(delta)", axis, fit_depth(axis; weights), oracle; weights)
    end
    CSV.write(joinpath(outdir,"inputs.csv"), DataFrame(inputs))
    table = DataFrame(observations)
    CSV.write(joinpath(outdir,"observations.csv"), table)
    summary = combine(groupby(table,:example), :raw_error => maximum => :raw_error,
        :percentile_error => maximum => :percentile_error)
    CSV.write(joinpath(outdir,"summary.csv"), summary)
    show(stdout, MIME"text/plain"(), summary); println()
    println("Independent exact premises verified; recorded differences are investigative findings.")
end

main(isempty(ARGS) ? joinpath(@__DIR__, "results") : only(ARGS))
