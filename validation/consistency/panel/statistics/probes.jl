# Review-only probes; no package code is modified.
# julia --startup-file=no --project=. validation/consistency/panel/statistics/probes.jl
using Platypus, CSV, DataFrames, Distributions, Statistics, SHA

const ROOT = normpath(joinpath(@__DIR__, "../../../.."))
const OUT = @__DIR__

# Surgical diagnostic: reproduce density.jl in an isolated module, changing only
# multiplication/addition grouping to make coordinate exchange commutative.
# This is causal evidence, not a proposed complete numerical repair.
module SymmetricPrototype
using Distributions
const source_path = normpath(joinpath(@__DIR__, "../../../../src/density.jl"))
const source = replace(read(source_path, String),
    "log(p[j]) + lu + lv" => "log(p[j]) + (lu + lv)",
    "p[i]*(zu[i]-mu)*(zv[i]-mv)" => "p[i]*((zu[i]-mu)*(zv[i]-mv))",
    "2rho*z1*z2" => "2rho*(z1*z2)")
include_string(@__MODULE__, source, "isolated_symmetric_density_prototype.jl")
end

function fit(M; family=:gaussian, weights=nothing, regularization=.05)
    fit_copula(build_pairs(EntitySet(M, [:x,:y])); family,
        weights=weights === nothing ? nothing : ReferenceWeights(weights),
        regularization=family === :beta ? 0.0 : regularization)
end

function main()
    open(joinpath(OUT,"environment.txt"),"w") do io
        println(io,"Julia ",VERSION)
        for path in ("src/density.jl","src/pairs.jl","src/estimation.jl")
            println(io,path," SHA256 ",bytes2hex(sha256(read(joinpath(ROOT,path)))))
        end
    end
    evidence=NamedTuple[]

    # Analytic zero-dimensional KDE support: density is exactly 1 for any
    # positive weights, including every leave-one-out reference.
    M=ones(3,2); weights=[1.,2.,3.]
    r=fit_typicality(build_pairs(EntitySet(M,[:x,:y]));method=:kde,
                     weights=ReferenceWeights(weights))
    for i in 1:3
        push!(evidence,(case="constant_kde",entity=i,oracle_raw=0.,
            package_raw=r.raw_scores[i,1],oracle_score=.5,package_score=r.pair_scores[i,1]))
    end

    # Independent empirical-beta density oracle. Integer-shape beta densities
    # are polynomials with rational coefficients; no package density helper,
    # logpdf, or floating quantile is needed.
    n=3; u=[BigInt(2i-1)//BigInt(2n) for i in 1:n]
    beta_pdf(x,r)=n*binomial(BigInt(n-1),r-1)*x^(r-1)*(1-x)^(n-r)
    exact=[sum(beta_pdf(x,j)^2 for j in 1:n)/n for x in u]
    @assert exact[1] == exact[3] && exact[1] > exact[2]
    M=hcat(collect(1.:3.),collect(1.:3.)); r=fit(M;family=:beta); reflected=fit(-M;family=:beta)
    oracle_scores=[2/3,1/6,2/3]
    for i in 1:3
        push!(evidence,(case="beta_collinear",entity=i,oracle_raw=log(Float64(exact[i])),
            package_raw=r.raw_scores[i,1],oracle_score=oracle_scores[i],package_score=r.pair_scores[i,1]))
    end
    CSV.write(joinpath(OUT,"beta_exact.csv"),DataFrame(entity=1:3,
        exact_density=string.(exact),raw=r.raw_scores[:,1],reflected_raw=reflected.raw_scores[:,1],
        score=r.pair_scores[:,1],reflected_score=reflected.pair_scores[:,1],oracle_score=oracle_scores))

    # Uniform 5x5 Cartesian support is exactly an independent empirical
    # distribution. Any deterministic marginal map preserves factorization;
    # Gaussian normal-score covariance is therefore exactly zero. Appending
    # zero-mass copies must preserve it.
    M=reduce(vcat,[reshape([Float64(x),Float64(y)],1,2) for x in 1:5 for y in 1:5])
    augmented=vcat(M,M[1:2,:]); w=vcat(ones(25),zeros(2)); grid_r=fit(M); grid_a=fit(augmented;weights=w)
    columns=Platypus._transform_columns(augmented,:rank,w)
    zu=quantile.(Normal(),columns[1][1:25]); zv=quantile.(Normal(),columns[2][1:25])
    # Treat the rounded normal scores themselves as exact rational numbers.
    # Covariance still vanishes exactly: this isolates reduction arithmetic
    # from approximation of the normal quantile.
    qu=Rational{BigInt}.(zu); qv=Rational{BigInt}.(zv)
    exact_cov=sum((qu.-sum(qu)/25).*(qv.-sum(qv)/25))/25
    @assert exact_cov == 0
    for i in 1:25
        push!(evidence,(case="gaussian_independent_grid_zero_queries",entity=i,oracle_raw=0.,
            package_raw=grid_a.raw_scores[i,1],oracle_score=.5,package_score=grid_a.pair_scores[i,1]))
    end
    CSV.write(joinpath(OUT,"gaussian_grid.csv"),DataFrame(entity=1:25,x=M[:,1],y=M[:,2],
        baseline_raw=grid_r.raw_scores[:,1],augmented_raw=grid_a.raw_scores[1:25,1],
        baseline_score=grid_r.pair_scores[:,1],augmented_score=grid_a.pair_scores[1:25,1]))
    CSV.write(joinpath(OUT,"oracles.csv"),DataFrame(evidence))

    # Reproduce every existing KDE/Gaussian feature-exchange flag and check the
    # arithmetic-only isolated prototype on exactly the same coordinates.
    checks=CSV.read(joinpath(ROOT,"validation/consistency/results/checks.csv"),DataFrame)
    selected=filter(x->x.status=="violation" && x.check=="feature_permutation" &&
        x.method in ("kde","gaussian","gaussian_regularized"),checks)
    results=NamedTuple[]
    for row in eachrow(selected)
        d=CSV.read(joinpath(ROOT,"validation/consistency/results/datasets",row.dataset*".csv"),DataFrame)
        names=filter(x->occursin(r"^f\d+$",x),DataFrames.names(d))
        X=Matrix{Float64}(d[:,names]); n,p=size(X); w=fill(1/n,n)
        cols=Platypus._transform_columns(X,Symbol(row.transform))
        baseline_raw=0.; baseline_score=0.; prototype_raw=0.; prototype_score=0.
        for j in 1:p-1, k in j+1:p
            if row.method=="kde"
                kw=(support=row.transform=="rank" ? :unit : :real,logscale=true)
                f=Platypus._density_scores; g=SymmetricPrototype._density_scores
            else
                kw=(regularization=row.method=="gaussian" ? 0.0 : .05,)
                f=Platypus._copula_scores; g=SymmetricPrototype._copula_scores
            end
            a=f(cols[j],cols[k],w;kw...)[1]; b=f(cols[k],cols[j],w;kw...)[1]
            c=g(cols[j],cols[k],w;kw...)[1]; d=g(cols[k],cols[j],w;kw...)[1]
            delta(x,y)=maximum(a==b ? 0. : abs(a-b) for (a,b) in zip(x,y))
            baseline_raw=max(baseline_raw,delta(a,b)); prototype_raw=max(prototype_raw,delta(c,d))
            baseline_score=max(baseline_score,maximum(abs.(Platypus._score_percentiles(a,w).-Platypus._score_percentiles(b,w))))
            prototype_score=max(prototype_score,maximum(abs.(Platypus._score_percentiles(c,w).-Platypus._score_percentiles(d,w))))
        end
        push!(results,(dataset=row.dataset,method=row.method,transform=row.transform,
            baseline_raw,baseline_score,prototype_raw,prototype_score))
    end
    CSV.write(joinpath(OUT,"feature_exchange_attribution.csv"),DataFrame(results))
    open(joinpath(OUT,"summary.txt"),"w") do io
        println(io,"Constant KDE raw: ",evidence[1:3])
        println(io,"Beta exact densities: ",exact)
        println(io,"Gaussian grid exact rational covariance of rounded normal scores: ",exact_cov)
        println(io,"Gaussian grid rho baseline: ",grid_r.configuration.pair_parameters[1].rho)
        println(io,"Gaussian grid rho with zero-weight copies: ",grid_a.configuration.pair_parameters[1].rho)
        println(io,"Gaussian grid max score change: ",maximum(abs.(grid_r.pair_scores.-grid_a.pair_scores[1:25,:])))
        println(io,"Existing feature-exchange flagged configurations: ",length(results))
        println(io,"Original feature-exchange score violations: ",count(x->x.baseline_score>1e-10,results))
        println(io,"Prototype remaining feature-exchange score violations: ",count(x->x.prototype_score>1e-10,results))
        println(io,"Prototype maximum raw feature-exchange difference: ",maximum(x.prototype_raw for x in results))
    end
    print(read(joinpath(OUT,"summary.txt"),String))
end
main()
