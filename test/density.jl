using Test

@testset "Weighted density and copula numerical oracles" begin
    @testset "Reflected kernel conserves support mass" begin
        # Trapezoid integration is independent of the analytic implementation.
        grid = collect(range(0.0, 1.0; length=20001))
        for center in (0.0, 0.02, 0.5, 0.98, 1.0), h in (0.05, 0.3, 1.0)
            y = [Platypus._reflected_kernel(x, center, h) for x in grid]
            mass = (sum(y) - (first(y)+last(y))/2) / (length(grid)-1)
            @test mass ≈ 1.0 atol=1e-6
        end
        @test_throws ArgumentError Platypus._reflected_kernel(0.5, 0.5, 1.1)
    end

    @testset "KDE membership and invariance" begin
        u = [0.1, 0.2, 0.45, 0.6, 0.9]
        v = [0.8, 0.1, 0.6, 0.3, 0.4]
        w = [1.0, 3.0, 2.0, 1.0, 0.5]
        for support in (:unit, :real)
            score, pars = Platypus._density_scores(u,v,w; support=support)
            perm = [5,3,1,4,2]
            permuted, _ = Platypus._density_scores(u[perm],v[perm],w[perm]; support=support)
            @test permuted[invperm(perm)] == score
            for factor in (1e-200, 1e200)
                scaled, _ = Platypus._density_scores(u,v,w .* factor; support=support)
                @test scaled ≈ score rtol=2e-14
            end
            extra = support === :unit ? 1.0 : 1e100
            augmented, augmented_pars = Platypus._density_scores(
                vcat(u,extra),vcat(v,extra),vcat(w,0.0); support=support)
            @test augmented[1:5] == score
            @test augmented_pars.bandwidth == pars.bandwidth
            logs, _ = Platypus._density_scores(u,v,w; support=support,logscale=true)
            @test exp.(logs) ≈ score
        end
        a, _ = Platypus._density_scores(u,v,ones(5))
        b, _ = Platypus._density_scores(u,v,fill(7.0,5))
        @test a == b
        @test_throws ArgumentError Platypus._density_scores(u,v,[1,0,0,0,0])
        @test_throws ArgumentError Platypus._density_scores(u,v,[1,1,1,1,NaN])
        @test_throws ArgumentError Platypus._density_scores(u,v,w; bandwidth=0)
        @test_throws ArgumentError Platypus._density_scores(u,v,w; bandwidth=(0.2,1.2))
        @test_throws ArgumentError Platypus._density_scores(u,v,w; support=:unknown)
    end

    @testset "KDE analytic values and independent-grid boundary" begin
        # Two real-space observations: each LOO density is the other's kernel.
        s, p = Platypus._density_scores([0.,1.],[0.,1.],[1.,1.];
                                      support=:real,bandwidth=(1.0,2.0))
        expected = exp(-0.5*(1+0.25))/(4pi)
        @test all(x -> isapprox(x,expected; rtol=2e-14), s)
        @test p.bandwidth == (1.0,2.0)
        @test p.raw_scale == :density
        logscore, p = Platypus._density_scores([0.,1.],[0.,1.],[1.,1.];
                                      support=:real,bandwidth=(1.0,2.0),logscale=true)
        @test p.raw_scale == :log_density
        @test exp.(logscore) ≈ s
        # Uniform product grid: no corner depression relative to the center.
        t = (collect(1:41) .- 0.5) ./ 41
        u, v = repeat(t,41), repeat(t; inner=41)
        d, _ = Platypus._density_scores(u,v,ones(length(u));bandwidth=0.2)
        @test 0.94 < d[1]/d[841] < 1.06
        @test all(x -> 0.9 < x < 1.1, d)
    end

    @testset "Constant coordinates have explicit reduced-support semantics" begin
        d, pars = Platypus._density_scores(fill(0.5,3),fill(0.5,3),ones(3))
        @test d == ones(3)
        @test pars.constant == (true,true)
        d, _ = Platypus._density_scores([0.5,0.5,0.9],[0.5,0.5,0.5],[1.,1.,0.])
        @test d == [1.,1.,0.]
        d, pars = Platypus._density_scores([0.5,0.5],[0.,1.],ones(2);
                                          support=:real,bandwidth=1.0)
        @test d ≈ fill(exp(-0.5)/sqrt(2pi),2)
        @test pars.constant == (true,false)
    end

    @testset "Gaussian copula analytic identities and fit" begin
        for z1 in (-2.,0.,1.), z2 in (-1.,0.,3.)
            @test Platypus._gaussian_copula_logdensity(z1,z2,0.0) == 0.0
            for rho in (-0.7,0.4)
                s = Platypus._gaussian_copula_logdensity(z1,z2,rho)
                @test s ≈ Platypus._gaussian_copula_logdensity(z2,z1,rho)
                @test s ≈ Platypus._gaussian_copula_logdensity(-z1,-z2,rho)
                @test s ≈ Platypus._gaussian_copula_logdensity(z1,-z2,-rho)
            end
        end
        @test exp(Platypus._gaussian_copula_logdensity(0.,0.,0.5)) ≈ 2/sqrt(3)
        @test Platypus._gaussian_copula_logdensity(1.,-1.,0.5) ≈ -0.5log(0.75)-1.0
        # A symmetric product grid has fitted rho=0 exactly.
        u = repeat([0.125,0.375,0.625,0.875],4)
        v = repeat([0.125,0.375,0.625,0.875];inner=4)
        scores, pars = Platypus._copula_scores(u,v,ones(16))
        @test abs(pars.rho) < 1e-15
        @test all(x -> abs(x) < 1e-14,scores)
        # Strong dependence requires explicit regularization at singularity.
        x = [0.1,0.3,0.7,0.9]
        @test_throws ArgumentError Platypus._copula_scores(x,x,ones(4))
        scores, pars = Platypus._copula_scores(x,x,ones(4);regularization=0.05)
        @test pars.rho ≈ 0.95
        @test all(isfinite,scores)
        @test_throws ArgumentError Platypus._copula_scores(vcat(x,0.),vcat(x,1.),[1.,1.,1.,1.,0.];regularization=0.05)
        scores2, pars2 = Platypus._copula_scores(vcat(x,0.),vcat(x,1.),[1.,1.,1.,1.,0.];regularization=0.05,query_clip=0.01)
        @test scores2[1:4] == scores
        @test pars2.rho == pars.rho
        @test isfinite(scores2[end])
        @test pars2.query_clip == 0.01
        @test_throws ArgumentError Platypus._copula_scores(fill(0.5,4),x,ones(4))
        @test_throws ArgumentError Platypus._copula_scores(x,x,ones(4);regularization=1)
    end

    @testset "Empirical beta copula exact small-n oracle" begin
        # n=2 concordant density: 2*((1-u)*(1-v)+u*v).
        u = [0.25,0.75]
        scores, pars = Platypus._copula_scores(u,u,ones(2);family=:beta)
        @test exp.(scores) ≈ [1.25,1.25]
        @test pars.in_sample
        scores2, _ = Platypus._copula_scores(u,reverse(u),ones(2);family=:beta)
        @test exp.(scores2) ≈ [1.25,1.25]
        # Marginal integration of this analytic density gives exactly one.
        grid = collect(range(0.,1.;length=1001))
        for x in (0.,0.25,0.75,1.)
            y = [2*((1-x)*(1-t)+x*t) for t in grid]
            @test (sum(y)-(first(y)+last(y))/2)/1000 ≈ 1.0
        end
        @test_throws ArgumentError Platypus._copula_scores([0.5,0.5],u,ones(2);family=:beta)
        @test_throws ArgumentError Platypus._copula_scores(u,u,[1.,2.];family=:beta)
        @test_throws ArgumentError Platypus._copula_scores([0.2,0.5,0.8],[0.2,0.5,0.8],[1.,1.,0.];family=:beta)
    end
end
