# Independent covariance algebra and exact fixtures for numerical depth repairs.
@testset "Exact depth consistency" begin
    function depth_oracle(X, weights=ones(Int,size(X,1)))
        R = Rational{BigInt}.(X)
        p = Rational{BigInt}.(weights)
        p ./= sum(p)
        Z = R .- reshape(vec(p' * R), 1, :)
        C = Z' * Diagonal(p) * Z
        [1 / (1 + dot(Z[i,:], C \ Z[i,:])) for i in axes(Z,1)]
    end
    function exact_depth_fit(X; transform=:none, weights=nothing, rtol=nothing)
        fit_typicality(build_pairs(EntitySet(X, [:x,:y]); transform);
            method=:depth, depth_rtol=rtol,
            weights=weights === nothing ? nothing : ReferenceWeights(weights))
    end

    triangle = [0. 0.;1. 0.;0. 1.]
    @test depth_oracle(triangle) == fill(1//3,3)
    for transform in (:none,:rank,:zscore)
        r = exact_depth_fit(triangle; transform)
        @test vec(r.raw_scores) == fill(Float64(1//3),3)
        @test vec(r.pair_scores) == fill(.5,3)
        augmented = exact_depth_fit(vcat(triangle,triangle[1:1,:]);
            transform, weights=[1.,1.,1.,0.])
        @test augmented.raw_scores[1:3,:] == r.raw_scores
        @test augmented.pair_scores[1:3,:] == r.pair_scores
    end
    weighted = exact_depth_fit(triangle;weights=[1.,2.,3.])
    @test vec(weighted.raw_scores) ≈ Float64.(depth_oracle(triangle,[1,2,3])) rtol=2eps()

    # Two exact depth groups separated by approximately 1e-13 must not be
    # merged while repairing the triangle's exact equality.
    near = [-1. -1.;1. -1.;1. 1.;-1. 1. + 2.0^-40]
    oracle = depth_oracle(near)
    @test length(unique(oracle)) == 2
    @test 0 < maximum(oracle)-minimum(oracle) < 1//10^12
    r = exact_depth_fit(near)
    @test vec(r.raw_scores) == Float64.(oracle)
    @test vec(r.pair_scores) == [.25,.75,.75,.25]

    offset = [1.421875 -0.1875;0.5625 2.078125;0.421875 -1.65625;
              1.625 0.3125;0.171875 -0.71875]
    shifted = offset .+ 1e14
    @test shifted .- 1e14 == offset
    oracle = Float64.(depth_oracle(offset))
    for X in (offset, shifted, -offset, offset[:,[2,1]],
              offset .* 2.0^-500, offset .* 2.0^500)
        @test vec(exact_depth_fit(X).raw_scores) == oracle
    end
    order = [5,2,4,1,3]
    @test vec(exact_depth_fit(shifted[order,:]).raw_scores)[invperm(order)] == oracle

    for delta in (2.0^-30,2.0^-24,2.0^-20,2.0^-16)
        line = [-1. -1.;0. 0.;1. 1.;0. delta]
        shear = [1. -1.;0. 1.]
        expected = [.4,1.,.4,0.]
        for X in (line,line*shear,-line,line[:,[2,1]])
            @test vec(exact_depth_fit(X;weights=[1.,1.,1.,0.]).raw_scores) == expected
        end
    end

    # Exact full-rank data remains full rank even when a double precision
    # eigenvalue threshold would discard its second direction.
    thin = [0. 0.;1. 0.;0. 2.0^-500]
    @test vec(exact_depth_fit(thin).raw_scores) == fill(Float64(1//3),3)
    @test_throws ArgumentError exact_depth_fit(thin;rtol=1e-12)
    extreme = [-floatmax(Float64) 0.;floatmax(Float64) 0.;0. floatmax(Float64)]
    @test Platypus._maha_depth(extreme[:,1],extreme[:,2],ones(3)) == fill(Float64(1//3),3)
    point = [2. 3.;2. 3.;2. nextfloat(3.);nextfloat(2.) 3.]
    @test vec(exact_depth_fit(point;weights=[1.,1.,0.,0.]).raw_scores) == [1.,1.,0.,0.]
    raw, metadata = Platypus._maha_depth(triangle[:,1],triangle[:,2],ones(3);
        return_metadata=true)
    @test metadata.effective_rank == 2
    @test metadata.rank_policy == :exact_covariance_rank
    @test metadata.support_policy == :exact_affine_span
    @test !metadata.rank_truncation
    @test_throws ArgumentError Platypus._maha_depth([0.,1.],[0.,1.],ones(2);rtol=0)
end
