using Random, Statistics, StatsBase

const SCENARIOS = (
    :uniform, :gaussian, :heavy_tails, :lognormal, :positive_dependence,
    :negative_dependence, :parabola, :noisy_ring, :two_moons,
    :unequal_clusters, :bridge_cluster, :dependence_violations,
    :ordinal, :zero_inflated, :duplicate_cloud, :constant_cloud,
    :constant_feature, :collinear, :nearly_collinear, :cartesian_grid,
    :parity, :mixed_units, :large_offset, :tiny_units, :huge_units,
    :many_features, :two_entities, :three_entities,
)

"""Seeded numerical stress fixtures; planted labels are diagnostic, not an oracle."""
function dataset(scenario, seed; n=63)
    rng = MersenneTwister(seed)
    truth = falses(n)
    x = randn(rng,n)
    M = if scenario === :uniform
        rand(rng,n,4)
    elseif scenario === :gaussian
        randn(rng,n,4)
    elseif scenario === :heavy_tails
        tan.(pi .* (rand(rng,n,4) .- 0.5))
    elseif scenario === :lognormal
        exp.(2 .* randn(rng,n,4))
    elseif scenario === :positive_dependence
        x .+ 0.05 .* randn(rng,n,4)
    elseif scenario === :negative_dependence
        hcat(x,-x .+ 0.05randn(rng,n),randn(rng,n))
    elseif scenario === :parabola
        hcat(x,x.^2 .+ 0.01randn(rng,n),sin.(x) .+ 0.01randn(rng,n))
    elseif scenario === :noisy_ring
        angle=2pi .* rand(rng,n); radius=1 .+ 0.01randn(rng,n)
        hcat(radius.*cos.(angle),radius.*sin.(angle))
    elseif scenario === :two_moons
        angle=pi .* rand(rng,n); second=collect(1:n) .> n÷2
        hcat(cos.(angle) .+ second, sin.(angle).*(1 .- 2second) .+ 0.2second) .+
            0.01randn(rng,n,2)
    elseif scenario in (:unequal_clusters,:bridge_cluster)
        m=max(3,n÷10); a=n-m-n÷4; b=n-a-m
        truth[end-m+1:end] .= true
        vcat(-2 .+ 0.05randn(rng,a,4), 2 .+ 0.5randn(rng,b,4),
             (scenario === :bridge_cluster ? 0. : 6.) .+ 0.02randn(rng,m,4))
    elseif scenario === :dependence_violations
        A=x .+ 0.1randn(rng,n,4); m=max(3,n÷10)
        A[end-m+1:end,[2,4]] .*= -1; truth[end-m+1:end] .= true; A
    elseif scenario === :ordinal
        Float64.(rand(rng,0:4,n,4))
    elseif scenario === :zero_inflated
        randn(rng,n,4) .* (rand(rng,n,4) .> 0.85)
    elseif scenario === :duplicate_cloud
        centers=randn(rng,5,3); centers[rand(rng,1:5,n),:]
    elseif scenario === :constant_cloud
        fill(3.7,n,3)
    elseif scenario === :constant_feature
        hcat(randn(rng,n,2),fill(3.7,n))
    elseif scenario === :collinear
        hcat(x,2x,-3x)
    elseif scenario === :nearly_collinear
        hcat(x,x .+ 1e-10randn(rng,n),-x .+ 1e-10randn(rng,n))
    elseif scenario === :cartesian_grid
        side=round(Int,cbrt(n)); grid=collect(range(-1.,1.;length=side))
        reduce(vcat,([a b c] for a in grid for b in grid for c in grid))
    elseif scenario === :parity
        common=[0. 0. 0.;0. 1. 1.;1. 0. 1.;1. 1. 0.]
        rare=[0. 0. 1.;0. 1. 0.;1. 0. 0.;1. 1. 1.]
        A=vcat(repeat(common;inner=(max(2,n÷4),1)),rare)
        truth=vcat(falses(size(A,1)-4),trues(4)); A
    elseif scenario === :mixed_units
        randn(rng,n,3) .* [1e-150 1. 1e150]
    elseif scenario === :large_offset
        1e14 .+ randn(rng,n,3)
    elseif scenario === :tiny_units
        1e-150 .* randn(rng,n,3)
    elseif scenario === :huge_units
        1e150 .* randn(rng,n,3)
    elseif scenario === :many_features
        randn(rng,24,12)
    elseif scenario === :two_entities
        randn(rng,2,3)
    elseif scenario === :three_entities
        randn(rng,3,3)
    else
        error("unknown scenario: $scenario")
    end
    length(truth)==size(M,1) || (truth=falses(size(M,1)))
    weights=exp.(range(-3.,3.;length=size(M,1)))
    shuffle!(rng,weights)
    (matrix=M, truth=truth, weights=weights,
     names=Symbol.("f" .* string.(1:size(M,2))),
     labels=[isodd(i) ? "A" : "B" for i in axes(M,1)])
end
