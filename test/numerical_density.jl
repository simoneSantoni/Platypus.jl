# Independent identities and invariances for numerical density repairs.
using Test, Platypus, Random

@testset "Numerical density identities" begin
    for weights in ([1.,2.,3.], [1.,2.,3.].*1e150, [1.,2.,3.].*1e-150)
        logdensity, _ = Platypus._density_scores(ones(3),ones(3),weights;logscale=true)
        density, _ = Platypus._density_scores(ones(3),ones(3),weights)
        @test logdensity == zeros(3)
        @test density == ones(3)
    end
    logdensity, _ = Platypus._density_scores([1.,1.,1.,2.],ones(4),[1.,2.,3.,0.];support=:real,logscale=true)
    @test logdensity == [0.,0.,0.,-Inf]
    varying, _ = Platypus._density_scores([0.,1.,3.],ones(3),[1.,2.,3.];support=:real,logscale=true)
    @test length(unique(varying)) > 1

    u = Rational{BigInt}.([1//6,1//2,5//6])
    beta, _ = Platypus._copula_scores(u,u,ones(3);family=:beta)
    # Differentiate the exact Bernstein polynomials independently.
    basis(n,r,x) = n*binomial(big(n-1),r-1)*x^(r-1)*(1-x)^(n-r)
    oracle = [sum(basis(3,j,x)^2 for j in 1:3)/3 for x in u]
    @test oracle == [121//72,9//8,121//72]
    @test beta ≈ log.(Float64.(oracle))
    @test beta[1] == beta[3]
    @test beta[1] != beta[2]
    reflected, _ = Platypus._copula_scores(1 .-u,1 .-u,ones(3);family=:beta)
    @test reflected == beta

    # Rank mixtures with nontrivial permutations, including reflection in
    # only one margin; the Bernstein identity does not require a symmetric fit.
    for perm in ([2,4,1,3],[4,3,2,1],[1,3,4,2])
        x=Rational{BigInt}.((1:2:7).//8); y=x[perm]
        a,_=Platypus._copula_scores(x,y,ones(4);family=:beta)
        b,_=Platypus._copula_scores(1 .-x,y,ones(4);family=:beta)
        @test a == b
        exact=[sum(basis(4,j,x[i])*basis(4,perm[j],y[i]) for j in 1:4)/4 for i in 1:4]
        @test a ≈ log.(Float64.(exact))
    end
    # An independent exact-polynomial sweep also checks accidental equalities,
    # not just visually symmetric examples.
    for n in 3:9, seed in 1:10
        x=Rational{BigInt}.((1:2:2n-1).//(2n))
        perm=randperm(MersenneTwister(seed),n); y=x[perm]
        actual,_=Platypus._copula_scores(x,y,ones(n);family=:beta)
        expected=[sum(basis(n,j,x[i])*basis(n,perm[j],y[i]) for j in 1:n)/n for i in 1:n]
        @test all(sign(actual[i]-actual[j]) == sign(expected[i]-expected[j])
                  for i in 1:n for j in 1:i-1)
    end
    # A representable, genuinely unequal nearby value must not be rounded away.
    shifted=copy(u); shifted[3] += 1//big(2)^40
    close,_=Platypus._copula_scores(shifted,shifted,ones(3);family=:beta)
    exactclose=[sum(basis(3,j,x)^2 for j in 1:3)/3 for x in shifted]
    @test exactclose[1] != exactclose[3]
    @test close[1] != close[3]

    grid=[(x,y) for x in u for y in u]
    x=first.(grid); y=last.(grid)
    for weights in (ones(9), [a*b for a in (1,2,3) for b in (4,5,6)])
        a,params=Platypus._copula_scores(x,y,weights;regularization=.05)
        @test params.rho == 0
        @test a == zeros(9)
        b,aug=Platypus._copula_scores(vcat(x,x[1]),vcat(y,y[1]),vcat(weights,0);regularization=.05)
        @test aug.rho == 0
        @test b == zeros(10)
    end
    # Public fits must retain original weight provenance: individually rounded
    # normalized product weights need not themselves factor exactly.
    matrix=hcat(Float64.(x),Float64.(y))
    for weights in (ones(9), [a*b for a in (1,2,3) for b in (4,5,6)])
        pairs=build_pairs(EntitySet(matrix,[:x,:y]))
        result=fit_copula(pairs;weights=ReferenceWeights(weights),regularization=.05)
        @test result.raw_scores == zeros(9,1)
        @test result.pair_scores == fill(.5,9,1)
        augmented=build_pairs(EntitySet(vcat(matrix,matrix[1:1,:]),[:x,:y]))
        result=fit_copula(augmented;weights=ReferenceWeights(vcat(weights,0)),regularization=.05)
        @test result.raw_scores == zeros(10,1)
        @test result.pair_scores == fill(.5,10,1)
    end
    beta_public=fit_copula(build_pairs(EntitySet([1. 1.;2. 2.;3. 3.],[:x,:y]));family=:beta)
    @test vec(beta_public.pair_scores) == [2/3,1/6,2/3]
    for transform in (:rank,:none,:zscore)
        constant=fit_typicality(build_pairs(EntitySet(ones(3,2),[:x,:y]);transform);
            method=:kde,weights=ReferenceWeights([1,2,3]))
        @test constant.pair_scores == fill(.5,3,1)
    end
    # Weak but nonzero dependence is preserved by exact covariance accumulation.
    weights=ones(9); weights[end] += 2.0^-40
    _,weak=Platypus._copula_scores(x,y,weights;regularization=.05)
    @test weak.rho != 0

    # Caller BigFloat settings must not change numerical results.
    baseline,_=Platypus._density_scores(u,u,ones(3);logscale=true)
    correlation,_=Platypus._copula_scores(u,u,ones(3);regularization=.05)
    setprecision(BigFloat,16) do
        reduced,_=Platypus._density_scores(u,u,ones(3);logscale=true)
        reduced_correlation,_=Platypus._copula_scores(u,u,ones(3);regularization=.05)
        @test reduced == baseline
        @test reduced_correlation == correlation
        tasks=[Threads.@spawn(first(Platypus._density_scores(u,u,ones(3);logscale=true))) for _ in 1:8]
        @test all(fetch(task) == baseline for task in tasks)
        tasks=[Threads.@spawn(first(Platypus._copula_scores(u,u,ones(3);regularization=.05))) for _ in 1:8]
        @test all(fetch(task) == correlation for task in tasks)
    end
    for exponent in (-2000,-1000,-100,0,100,1000,2000)
        q=exponent>=0 ? (big(1)<<exponent)//big(1) : big(1)//(big(1)<<(-exponent))
        @test Platypus._density_sqrt_rational(q) == 2.0^(exponent÷2)
    end
    for support in (:unit,:real)
        x=support===:unit ? u : [1e14,1e14+1,1e14+3]
        y=support===:unit ? u[[2,3,1]] : [1e14+2,1e14-2,1e14]
        a,_=Platypus._density_scores(x,y,ones(3);support,logscale=true)
        b,_=Platypus._density_scores(y,x,ones(3);support,logscale=true)
        @test a == b
        rx=support===:unit ? 1 .-x : -x
        ry=support===:unit ? 1 .-y : -y
        c,_=Platypus._density_scores(rx,ry,ones(3);support,logscale=true)
        @test a == c
        perm=[3,1,2]
        d,_=Platypus._density_scores(x[perm],y[perm],ones(3);support,logscale=true)
        @test a == d[invperm(perm)]
        e,_=Platypus._density_scores(vcat(x,x[1]),vcat(y,y[1]),[1,1,1,0];support,logscale=true)
        @test a == e[1:3]
    end
end

@testset "KDE ordering survives dimensional log rounding" begin
    # Four represented points suffice: changing physical units makes two
    # unequal log densities round to the same Float64 value.
    matrix=[-2.0 0.;0 0;1e-14 0;1 0]
    baseline=fit_typicality(build_pairs(EntitySet(matrix,[:x,:y]);transform=:none);method=:kde)
    h=baseline.configuration.pair_parameters[1].bandwidth[1]
    oracle=setprecision(BigFloat,512) do
        x=BigFloat.(matrix[:,1]); bandwidth=BigFloat(h)
        [sum(exp(-((x[i]-x[j])/bandwidth)^2/2) for j in 1:4 if i!=j)/3 for i in 1:4]
    end
    @test oracle[3] > oracle[2]
    expected=[sum(oracle[j]<oracle[i] ? 1. : oracle[j]==oracle[i] ? .5 : 0. for j in 1:4)/4 for i in 1:4]
    @test vec(baseline.pair_scores) == expected
    @test baseline.normalization_scores[3] > baseline.normalization_scores[2]
    for exponent in (-200,200)
        changed=fit_typicality(build_pairs(EntitySet(matrix.*2.0^exponent,[:x,:y]);transform=:none);method=:kde)
        @test changed.raw_scores[2] == changed.raw_scores[3]
        @test changed.normalization_scores == baseline.normalization_scores
        @test changed.pair_scores == baseline.pair_scores
        @test changed.raw_scores .+ exponent*log(2) ≈ baseline.raw_scores
    end
    # Zero-reference tail queries have all their direct exponential kernels
    # underflow; shifted dimensionless log masses must still order those tails.
    x=[0.,1.,100.,200.]; w=[1.,1.,0.,0.]
    raw,pars,order=Platypus._density_scores(x,zeros(4),w;
        support=:real,bandwidth=1.,logscale=true,return_ordering=true)
    @test all(isfinite,raw)
    @test order[4] < order[3] < order[1]
    @test pars.normalization_basis == :scale_free_kernel_order
    for exponent in (-200,200)
        _,_,changed=Platypus._density_scores(x.*2.0^exponent,zeros(4),w;
            support=:real,bandwidth=2.0^exponent,logscale=true,return_ordering=true)
        @test changed == order
    end
end

@testset "Weighted mixed-underflow KDE comparison" begin
    matrix=[0. 0.;100. 0.;38. 0.;139. 0.]
    for weights in ([1e-300,1.,0.,1e-200],[1e-20,1.,0.,1e-10])
        oracle,percentiles=setprecision(BigFloat,2048) do
            q=BigFloat.(weights); x=BigFloat.(matrix[:,1])
            densities=[sum(q[j]*exp(-(x[i]-x[j])^2/2) for j in 1:4 if j!=i)/
                sum(q[j] for j in 1:4 if j!=i) for i in 1:4]
            cdf=[sum(q[j]*(densities[j]<densities[i] ? 1 : densities[j]==densities[i] ? .5 : 0)
                for j in 1:4)/sum(q) for i in 1:4]
            densities,Float64.(cdf)
        end
        @test oracle[1] < oracle[3] < oracle[2] < oracle[4]
        baseline=nothing
        for exponent in (0,-200,200)
            scale=2.0^exponent
            result=fit_typicality(build_pairs(EntitySet(matrix.*scale,[:x,:y]);transform=:none);
                method=:kde,bandwidth=scale,weights=ReferenceWeights(weights))
            @test vec(result.pair_scores) == percentiles
            @test all(sign(result.normalization_scores[i]-result.normalization_scores[j]) ==
                sign(oracle[i]-oracle[j]) for i in 1:4 for j in 1:4)
            if baseline===nothing
                baseline=result.normalization_scores
            else
                @test result.normalization_scores == baseline
            end
        end
    end
    # A represented nonzero nearby kernel times tiny weight can be dominated
    # within the same row by an initially underflowed, distant heavy kernel.
    x=[0.,80.,38.,40.]; weights=[1e-300,1.,0.,0.]
    raw,_,order=Platypus._density_scores(x,zeros(4),weights;
        support=:real,bandwidth=1.,logscale=true,return_ordering=true)
    @test exp(-38.0^2/2)>0 && exp(-42.0^2/2)==0
    @test order == [1.,1.,2.,3.]
    @test raw[3] < raw[4]
    # Preserve terms beyond Float64's direct exponential range when exact
    # weights make them relevant; the half-exponential is squared rationally.
    @test exp(-746.) == 0
    @test Platypus._density_shifted_exp(-746.) > 0

    # Directed intervals must refine genuinely close comparisons, without
    # changing the caller's precision or assigning approximate equality.
    t=big(1)//(big(1)<<400); q=1+t; target=t-t^2/2
    @test Platypus._density_logratio_compare(q,target) == 1
    @test Platypus._density_logratio_compare(inv(q),-target) == -1
    @test Platypus._density_logratio_compare(big(1)//big(1),big(0)//big(1)) == 0
    masses=Rational{BigInt}.([0,1//7,1,3,17]); offsets=[0.,-1000.,-1.,0.,-2.]
    cmp(i,j)=Platypus._density_kernel_compare(masses[i],offsets[i],masses[j],offsets[j])
    @test all(cmp(i,j)==-cmp(j,i) for i in 1:5 for j in 1:5)
    @test all(!(cmp(i,j)<0 && cmp(j,k)<0) || cmp(i,k)<0 for i in 1:5 for j in 1:5 for k in 1:5)
    expected=Platypus._density_kernel_order(masses,offsets)
    setprecision(BigFloat,16) do
        tasks=[Threads.@spawn(Platypus._density_kernel_order(masses,offsets)) for _ in 1:8]
        @test all(fetch(task)==expected for task in tasks)
        @test precision(BigFloat)==16
    end
end
