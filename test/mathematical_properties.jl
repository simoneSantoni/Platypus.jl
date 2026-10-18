using LinearAlgebra
using Statistics

# Public-API properties and small independent mathematical oracles. These
# checks validate numerical definitions, not agreement with human judgment.
@testset "Mathematical properties of public estimators" begin
    pairfit(X; transform=:rank, kwargs...) =
        fit_typicality(build_pairs(EntitySet(X, Symbol.(["x$j" for j in axes(X,2)]));
                                   transform=transform); kwargs...)

    @testset "LOF distance ties are independent of row order" begin
        X = [0. 0.; 0. 1.; 0. 2.; 1. 0.; 2. 1.]
        original = pairfit(X; transform=:none, method=:lof, k=2)
        for permutation in ([5,4,3,2,1], [3,1,5,2,4]), backend in (:brute,:kdtree)
            permuted = pairfit(X[permutation,:]; transform=:none,
                               method=:lof, k=2, backend=backend)
            @test permuted.raw_scores[invperm(permutation),:] == original.raw_scores
            @test permuted.pair_scores[invperm(permutation),:] == original.pair_scores
        end
    end

    @testset "No evidence of distinction gives midpercentile one half" begin
        for method in TYPICALITY_METHODS
            result = pairfit(fill(7.0,8,3); method=method)
            @test result.pair_scores == fill(0.5,8,3)
            @test result.atypicality == fill(0.5,8)
            @test result.dispersion == zeros(8)
        end
    end

    @testset "Relative membership and zero-weight queries" begin
        # Regression for the original absolute-epsilon weight-scaling defect.
        X = [.875 .125; .375 .875; .125 .375; .625 .625]
        weights = [.1,1.,10.,100.]
        for method in (:kde,:depth), transform in (:rank,:none,:zscore)
            original = pairfit(X; transform=transform, method=method,
                               weights=ReferenceWeights(weights))
            for scale in (1e-200,1e-14,1e14,1e200)
                scaled = pairfit(X; transform=transform, method=method,
                                 weights=ReferenceWeights(scale .* weights))
                @test scaled.raw_scores ≈ original.raw_scores rtol=1e-12
                @test scaled.pair_scores ≈ original.pair_scores rtol=1e-12
            end
            extended = pairfit(vcat(X,[-1e5 1e5; .4 .4]); transform=transform,
                method=method, weights=ReferenceWeights(vcat(weights,0.,0.)))
            @test extended.raw_scores[1:4,:] ≈ original.raw_scores rtol=1e-12
            @test extended.pair_scores[1:4,:] ≈ original.pair_scores rtol=1e-12
        end
        ps = build_pairs(EntitySet(X,[:x,:y]))
        for method in (:lof,:knn)
            @test_throws ArgumentError fit_typicality(ps; method=method,
                weights=ReferenceWeights(weights))
        end
        # Floating summation of unequal tied masses must not make row order a
        # model parameter. Weight normalization used to fail this fixture.
        Xties = hcat(Float64.(mod.(0:9,3)), Float64.(mod.(0:9,4)))
        unequal = [1.;fill(1e-16,9)]
        permutation = collect(10:-1:1)
        for method in (:kde,:depth), transform in (:rank,:zscore)
            original = pairfit(Xties;method=method,transform=transform,
                               weights=ReferenceWeights(unequal))
            reversed = pairfit(Xties[permutation,:];method=method,transform=transform,
                               weights=ReferenceWeights(unequal[permutation]))
            @test reversed.raw_scores[invperm(permutation),:] == original.raw_scores
            @test reversed.pair_scores[invperm(permutation),:] == original.pair_scores
        end
    end

    @testset "Mahalanobis depth agrees with direct covariance algebra" begin
        X = [-2. 1.; -1. -2.; 1. 3.; 4. -1.; 0. 0.]
        weights = [1.,2.,3.,4.,0.]
        p = weights ./ sum(weights)
        mu = vec(sum(X .* p; dims=1))
        centered = X .- permutedims(mu)
        covariance = centered' * Diagonal(p) * centered
        expected = [1/(1+dot(centered[i,:],covariance\centered[i,:]))
                    for i in axes(X,1)]
        original = pairfit(X; transform=:none, method=:depth,
                           weights=ReferenceWeights(weights))
        @test vec(original.raw_scores) ≈ expected rtol=1e-12
        # Changes of coordinate units, reflection, and a nonsingular shear.
        for A in (Diagonal([1e-100,1e100]), Diagonal([-1e100,1e-100]),
                  [2. 1.; -1. 3.])
            changed = pairfit(X*A; transform=:none, method=:depth,
                              weights=ReferenceWeights(weights))
            @test vec(changed.raw_scores) ≈ expected rtol=1e-12
            @test changed.pair_scores ≈ original.pair_scores rtol=1e-12
        end
        translated = pairfit(X .+ [100. -23.]; transform=:none, method=:depth,
                             weights=ReferenceWeights(weights))
        @test vec(translated.raw_scores) ≈ expected rtol=1e-12
    end

    @testset "Singular depth uses affine support, not a dimensionful ridge" begin
        # Reference mean=(0,0), variance of the intrinsic coordinate=2/3.
        X = [-1. -2.; 0. 0.; 1. 2.; 0. 1.; 2. 4.]
        weights = ReferenceWeights([1.,1.,1.,0.,0.])
        expected = [.4,1.,.4,0.,1/7]
        for A in (Diagonal([1.,1.]), Diagonal([1e-100,1e100]))
            result = pairfit(X*A; transform=:none,method=:depth,weights=weights)
            @test vec(result.raw_scores) ≈ expected rtol=1e-12 atol=1e-15
        end
        # A point mass has depth one on support, zero for off-support queries.
        result = pairfit([2. 3.;2. 3.;2. 4.;1. 3.]; transform=:none,
                         method=:depth,weights=ReferenceWeights([1.,1.,0.,0.]))
        @test vec(result.raw_scores) == [1.,1.,0.,0.]
        # Non-power-of-two normalized weights must still identify a point mass.
        Xconstant = vcat(fill(7.,10,2),[7. 8.])
        result = pairfit(Xconstant;transform=:none,method=:depth,
                         weights=ReferenceWeights(vcat([1.;fill(1e-16,9)],0.)))
        @test vec(result.raw_scores) == vcat(ones(10),0.)
    end

    @testset "Raw values and normalized summaries have distinct contracts" begin
        X = [0. 1. 2.;1. 3. 0.;4. 2. 1.;2. 0. 4.;3. 4. 3.]
        p = fill(1/5,5)
        for method in TYPICALITY_METHODS
            result = pairfit(X;method=method)
            # Independent quadratic mid-CDF oracle; intentionally not sorting.
            oracle = [sum(p[j] * (result.raw_scores[j,c] < result.raw_scores[i,c] ?
                        1.0 : result.raw_scores[j,c] == result.raw_scores[i,c] ? 0.5 : 0.0)
                        for j in 1:5) for i in 1:5,c in 1:3]
            @test result.pair_scores ≈ oracle atol=1e-15
            @test result.conventionality ≈ [median(oracle[i,:]) for i in 1:5]
            @test result.novelty ≈ [quantile(oracle[i,:],.1) for i in 1:5]
            @test result.atypicality ≈ 1 .- result.conventionality
            @test scoretable(result;raw=true).typicality == vec(permutedims(result.raw_scores))
        end
        @test_throws ArgumentError pairfit(X;normalize=false)
    end

    @testset "Gaussian independence and exact small beta copula" begin
        # Normal-score cross products cancel exactly by sign symmetry.
        ps = build_pairs(EntitySet(hcat(1:4,[2,4,1,3]),[:x,:y]))
        gaussian = fit_copula(ps;family=:gaussian)
        @test vec(gaussian.raw_scores) ≈ zeros(4) atol=1e-14
        @test gaussian.pair_scores == fill(.5,4,1)
        # For two concordant ranks, c(u,v)=2*((1-u)*(1-v)+u*v).
        two = build_pairs(EntitySet([1. 1.;2. 2.],[:x,:y]))
        beta = fit_copula(two;family=:beta)
        @test exp.(beta.raw_scores) ≈ fill(1.25,2,1) rtol=1e-14
        @test beta.pair_scores == fill(.5,2,1)
        @test_throws ArgumentError fit_copula(two;family=:gaussian)
        @test all(isfinite,fit_copula(two;family=:gaussian,regularization=.1).raw_scores)
        @test_throws ArgumentError fit_copula(two;family=:beta,
                                            weights=ReferenceWeights([1.,2.]))
        tied = build_pairs(EntitySet([0. 0.;0. 1.;1. 2.],[:x,:y]))
        @test_throws ArgumentError fit_copula(tied;family=:beta)
    end

    @testset "Bounded KDE does not reproduce Gaussian boundary depression" begin
        m = 21
        X = [((i-.5)/m) for i in 1:m for j in 1:m]
        Y = [((j-.5)/m) for i in 1:m for j in 1:m]
        result = pairfit(hcat(X,Y);method=:kde,bandwidth=.3)
        # Quadrature of a constant product density is near one throughout.
        # Finite-grid quadrature and LOO omission preclude exact equality.
        density = exp.(result.raw_scores)
        @test maximum(abs.(density .- 1)) < .06
        @test minimum(density)/maximum(density) > .94
    end

    @testset "Explicit higher-order limitation of pairwise summaries" begin
        # Binary illustration: each odd-parity joint cell is 100x rarer, but
        # all bivariate cells contain 101 entities. This is a limitation test,
        # not a claim that tied binary data satisfy continuous copula theory.
        common = [0. 0. 0.;0. 1. 1.;1. 0. 1.;1. 1. 0.]
        rare = [0. 0. 1.;0. 1. 0.;1. 0. 0.;1. 1. 1.]
        X = vcat(repeat(common;inner=(100,1)),rare)
        for method in (:knn,:lof)
            result = pairfit(X;method=method,k=20)
            @test result.pair_scores == fill(.5,404,3)
            @test result.atypicality == fill(.5,404)
        end
    end
end
