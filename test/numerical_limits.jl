@testset "Finite numerical limits" begin
    @testset "Outcome projection scales finite responses" begin
        es = EntitySet([sin(i*j)+i/20 for i in 1:20,j in 1:3],[:a,:b,:c])
        r = fit_typicality(build_pairs(es))
        for magnitude in (0.0,1e-300,1e300,1e308)
            fit = fit_outcome(r,fill(magnitude,20);predictors=())
            @test only(fit.coefficients) ≈ magnitude
            @test fit.fitted ≈ fill(magnitude,20)
            @test all(isfinite,fit.residuals)
            @test isnan(fit.r2)
        end
        response = collect(range(-0.75,0.75;length=20))
        baseline = fit_outcome(r,response;predictors=(:conventionality,))
        large = fit_outcome(r,response .* 1e308;predictors=(:conventionality,))
        @test large.coefficients ./ 1e308 ≈ baseline.coefficients
        @test large.fitted ./ 1e308 ≈ baseline.fitted
        @test large.residuals ./ 1e308 ≈ baseline.residuals
        @test large.r2 ≈ baseline.r2
        # With no intercept and bounded predictors, this response requires an
        # unrepresentable slope. Return a deliberate error, never an Inf fit.
        @test_throws ArgumentError fit_outcome(r,fill(floatmax(Float64),20);
            predictors=(:conventionality,),intercept=false)
    end

    @testset "Gaussian KDE uses representable standardized differences" begin
        es = EntitySet([-1e308 -1.;1e308 1.],[:x,:y])
        ps = build_pairs(es;transform=:none)
        fitted = fit_typicality(ps;method=:kde,bandwidth=(1e308,1.0))
        # Each leave-one-out density is one Gaussian product kernel evaluated
        # two bandwidths from its center in both coordinates.
        expected = -4 - log(1e308) - log(2pi)
        @test fitted.raw_scores ≈ fill(expected,2,1)
        @test fitted.pair_scores == fill(0.5,2,1)
        ordinary = fit_typicality(build_pairs(EntitySet([-1. -1.;1. 1.],[:x,:y]);
            transform=:none);method=:kde,bandwidth=(1.0,1.0))
        @test fitted.raw_scores .+ log(1e308) ≈ ordinary.raw_scores
        for invalid in (big"1e400",big"1e-400")
            @test_throws ArgumentError fit_typicality(ps;method=:kde,bandwidth=(invalid,1.0))
            constants = build_pairs(EntitySet(ones(3,2),[:x,:y]);transform=:none)
            @test_throws ArgumentError fit_typicality(constants;method=:kde,bandwidth=invalid)
        end
    end
end
