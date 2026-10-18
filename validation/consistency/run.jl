# julia --threads=2 --project validation/consistency/run.jl [output_directory]
# Investigative harness: record violations rather than stop on the first one.
module ConsistencyStudy
using Platypus, CSV, DataFrames, Random, Statistics, StatsBase, LinearAlgebra, SHA
include("datasets.jl")

const METHODS=(:lof,:knn,:kde,:depth,:gaussian,:gaussian_regularized,:beta,:graph)
const SCORE_TOL=1e-10
const RAW_TOL=1e-10

function fit(M,names,method; transform=:rank, weights=nothing, labels=nothing,
             ids=collect(1:size(M,1)), backend=:kdtree, parallel=false, k=0)
    es=EntitySet(M,names;ids=ids,labels=labels)
    w=weights===nothing ? nothing : ReferenceWeights(weights)
    method===:graph && return fit_graph_typicality(es;transform,weights=w,backend,k)
    ps=build_pairs(es;transform)
    if method in (:gaussian,:gaussian_regularized,:beta)
        return fit_copula(ps;family=method===:beta ? :beta : :gaussian,weights=w,
                          regularization=method===:gaussian_regularized ? 0.05 : 0.0)
    end
    fit_typicality(ps;method,weights=w,backend,parallel,k)
end

function aligned(r, names; rows=collect(1:r.n_entities))
    if r isa GraphTypicalityResult
        return (raw=reshape(r.stationary[rows],:,1), scores=reshape(r.typicality[rows],:,1),
                normalization=reshape(r.stationary[rows],:,1),
                summaries=reshape(r.atypicality[rows],:,1))
    end
    keys=[Set(pair) for pair in r.pair_names]
    columns=[findfirst(==(Set((names[a],names[b]))),keys)
             for a in 1:length(names)-1 for b in a+1:length(names)]
    (raw=r.raw_scores[rows,columns],scores=r.pair_scores[rows,columns],
     normalization=r.normalization_scores[rows,columns],
     summaries=hcat(r.atypicality[rows],r.conventionality[rows],r.novelty[rows],r.dispersion[rows]))
end

function error_metrics(a,b)
    size(a)==size(b) || return (absolute=Inf,relative=Inf)
    absolute=0.; relative=0.
    for (x,y) in zip(a,b)
        x==y && continue # includes equal signed infinities for zero log-density
        (!isfinite(x) || !isfinite(y)) && return (absolute=Inf,relative=Inf)
        delta=abs(x-y)
        absolute=max(absolute,delta)
        relative=max(relative,delta/max(1.,abs(x),abs(y)))
    end
    (absolute=absolute,relative=relative)
end

function expected_rejection(method,M,error)
    error isa ArgumentError || return false
    msg=sprint(showerror,error)
    if method===:beta
        return any(length(unique(M[:,j]))<size(M,1) for j in axes(M,2)) &&
               occursin("untied",msg)
    elseif method in (:gaussian,:gaussian_regularized)
        return (any(length(unique(M[:,j]))==1 for j in axes(M,2)) && occursin("nonconstant",msg)) ||
               (method===:gaussian && occursin("singular Gaussian",msg))
    end
    false
end

function main(outdir=joinpath(@__DIR__,"results"); seeds=(17,101,2026))
    mkpath(joinpath(outdir,"datasets"))
    checks=NamedTuple[]; baselines=NamedTuple[]; catalog=NamedTuple[]; detail=NamedTuple[]
    values=NamedTuple[]; sensitivity=NamedTuple[]
    function record(key,method,transform,check,status;raw_abs=NaN,raw_rel=NaN,
                    pair_delta=NaN,summary_delta=NaN,note="")
        push!(checks,(dataset=key,method=string(method),transform=string(transform),
            check=check,status=status,raw_max_abs=raw_abs,raw_max_relative=raw_rel,
            pair_max_delta=pair_delta,summary_max_delta=summary_delta,note=note))
    end
    for scenario in SCENARIOS, (s,seed) in enumerate(seeds)
        d=dataset(scenario,seed;n=(32,63,96)[mod1(s,3)])
        M,names=d.matrix,d.names; n,p=size(M); key="$(scenario)_$(seed)"
        table=DataFrame(M,names)
        insertcols!(table,1,:entity_id=>collect(1:n),:planted=>d.truth,
                    :label=>d.labels,:reference_weight=>d.weights)
        CSV.write(joinpath(outdir,"datasets",key*".csv"),table)
        push!(catalog,(dataset=key,scenario=string(scenario),seed=seed,n=n,features=p,
            duplicate_rows=n-length(unique(eachrow(M))),
            constant_features=count(j->length(unique(M[:,j]))==1,1:p),
            min_value=minimum(M),max_value=maximum(M)))
        permutation=randperm(MersenneTwister(seed+31),n)
        columns=reverse(collect(1:p))
        for transform in (:rank,:none,:zscore), method in METHODS
            transform!==:rank && method in (:gaussian,:gaussian_regularized,:beta) && continue
            r=try
                fit(M,names,method;transform)
            catch e
                status=expected_rejection(method,M,e) ? "expected_rejection" : "unexpected_error"
                push!(baselines,(dataset=key,method=string(method),transform=string(transform),
                    status=status,n_unique_scores=0,min_atypicality=NaN,max_atypicality=NaN,
                    note=sprint(showerror,e)))
                continue
            end
            a=aligned(r,names)
            push!(baselines,(dataset=key,method=string(method),transform=string(transform),
                status="ok",n_unique_scores=length(unique(r.atypicality)),
                min_atypicality=minimum(r.atypicality),max_atypicality=maximum(r.atypicality),note=""))
            for i in 1:n
                push!(values,(dataset=key,method=string(method),transform=string(transform),
                    entity_id=i,planted=d.truth[i],atypicality=r.atypicality[i]))
            end
            valid=all(isfinite,a.scores) && all(x->0<=x<=1,a.scores) &&
                  all(isfinite,a.summaries) && all(x->0<=x<=1,a.summaries)
            record(key,method,transform,"finite_bounded_outputs",valid ? "pass" : "violation")
            w=fill(1/n,n)
            # KDE ordering can retain distinctions lost in physical Float64
            # log densities. Validate the declared, publicly exposed ordering
            # representation; its mathematical accuracy has separate oracles.
            oracle=[sum(w[j]*(a.normalization[j,c]<a.normalization[i,c] ? 1. : a.normalization[j,c]==a.normalization[i,c] ? .5 : 0.)
                        for j in 1:n) for i in 1:n,c in axes(a.normalization,2)]
            delta=error_metrics(oracle,a.scores).absolute
            record(key,method,transform,"independent_midcdf_oracle",delta<=SCORE_TOL ? "pass" : "violation";pair_delta=delta)
            if scenario===:constant_cloud
                delta=maximum(abs.(a.scores .- .5))
                record(key,method,transform,"constant_cloud_ties",delta<=SCORE_TOL ? "pass" : "violation";pair_delta=delta)
            end
            if method===:graph
                transition=Matrix(r.adjacency)
                for i in 1:n
                    r.degree[i]>0 ? (transition[i,:]./=r.degree[i]) : (transition[i,i]=1)
                end
                delta=maximum(abs.(transition'*r.stationary .- r.stationary))
                valid=delta<=SCORE_TOL && isapprox(sum(r.stationary),1.;atol=SCORE_TOL)
                record(key,method,transform,"stationarity",valid ? "pass" : "violation";raw_abs=delta)
            end
            function compare(check, thunk; rows=collect(1:n), raw_adjust=identity)
                try
                    other=thunk(); b=aligned(other,names;rows)
                    raw=error_metrics(a.raw,raw_adjust(b.raw))
                    pair=error_metrics(a.scores,b.scores).absolute
                    summary=error_metrics(a.summaries,b.summaries).absolute
                    status=raw.relative<=RAW_TOL && pair<=SCORE_TOL && summary<=SCORE_TOL ? "pass" : "violation"
                    record(key,method,transform,check,status;raw_abs=raw.absolute,raw_rel=raw.relative,pair_delta=pair,summary_delta=summary)
                    if status=="violation"
                        i,c=Tuple(argmax(abs.(a.scores .- b.scores)))
                        push!(detail,(dataset=key,method=string(method),transform=string(transform),check=check,
                            entity_id=i,pair_index=c,raw_before=a.raw[i,c],raw_after=raw_adjust(b.raw)[i,c],
                            score_before=a.scores[i,c],score_after=b.scores[i,c],
                            atypicality_before=a.summaries[i,1],atypicality_after=b.summaries[i,1]))
                    end
                catch e
                    record(key,method,transform,check,"unexpected_error";note=sprint(showerror,e))
                end
            end
            compare("row_permutation",()->fit(M[permutation,:],names,method;transform);rows=invperm(permutation))
            compare("feature_permutation",()->fit(M[:,columns],names[columns],method;transform))
            compare("label_and_id_change",()->fit(M,names,method;transform,
                labels=reverse(d.labels),ids=collect(1001:1000+n)))
            if transform===:rank
                recoded=hcat([tiedrank(M[:,j]).^3 for j in 1:p]...)
                compare("monotone_recode",()->fit(recoded,names,method;transform))
            end
            compare("coordinate_reflection",()->fit(-M,names,method;transform))
            if method in (:lof,:knn,:graph)
                compare("exact_backend",()->fit(M,names,method;transform,backend=:brute))
                for k in unique([1,n-1])
                    try
                        z=fit(M,names,method;transform,k)
                        push!(sensitivity,(dataset=key,method=string(method),transform=string(transform),
                            change="k=$k",max_atypicality_delta=maximum(abs.(r.atypicality.-z.atypicality))))
                    catch e
                        record(key,method,transform,"k=$k","unexpected_error";note=sprint(showerror,e))
                    end
                end
            end
            if method in (:lof,:knn,:kde,:depth)
                compare("serial_threaded",()->fit(M,names,method;transform,parallel=true))
            end
            if transform===:none
                for exponent in (-200,200)
                    scale=2.0^exponent
                    adjust = method===:knn ? (B->B./scale) : method===:kde ?
                        (B->B .+ reshape([count(!,par.constant)*log(scale) for par in r.configuration.pair_parameters],1,:)) : identity
                    compare("common_units_2^$exponent",()->fit(M.*scale,names,method;transform);raw_adjust=adjust)
                end
            end
            if method in (:kde,:depth,:gaussian_regularized)
                # Uniform baseline + excluded rows exercises a different
                # preprocessing branch from genuinely nonuniform references.
                extras=M[1:2,:]
                compare("append_zero_weight_uniform",()->fit(vcat(M,extras),names,method;
                    transform,weights=vcat(ones(n),0.,0.));rows=collect(1:n))
                try
                    weighted=fit(M,names,method;transform,weights=d.weights)
                    original=a
                    a=aligned(weighted,names)
                    for scale in (1e-150,1e150)
                        compare("weight_rescale_$scale",()->fit(M,names,method;transform,weights=d.weights.*scale))
                    end
                    compare("append_zero_weight_nonuniform",()->fit(vcat(M,extras),names,method;
                        transform,weights=vcat(d.weights,0.,0.));rows=collect(1:n))
                    a=original
                catch e
                    record(key,method,transform,"weighted_fit","unexpected_error";note=sprint(showerror,e))
                end
            end
        end
        println("Completed ",key," (",n," × ",p,")")
        flush(stdout)
    end
    CSV.write(joinpath(outdir,"datasets.csv"),DataFrame(catalog))
    CSV.write(joinpath(outdir,"baselines.csv"),DataFrame(baselines))
    CSV.write(joinpath(outdir,"checks.csv"),DataFrame(checks))
    details=isempty(detail) ? DataFrame(dataset=String[],method=String[],transform=String[],
        check=String[],entity_id=Int[],pair_index=Int[],raw_before=Float64[],raw_after=Float64[],
        score_before=Float64[],score_after=Float64[],atypicality_before=Float64[],atypicality_after=Float64[]) : DataFrame(detail)
    CSV.write(joinpath(outdir,"counterexamples.csv"),details)
    CSV.write(joinpath(outdir,"scores.csv"),DataFrame(values))
    CSV.write(joinpath(outdir,"sensitivity.csv"),DataFrame(sensitivity))
    counts=combine(groupby(DataFrame(checks),[:method,:transform,:check,:status]),nrow=>:count)
    CSV.write(joinpath(outdir,"summary.csv"),counts)
    open(joinpath(outdir,"environment.txt"),"w") do io
        println(io,"Julia ",VERSION,"; Platypus ",pkgversion(Platypus),"; threads=",Threads.nthreads())
        println(io,"Seeds: ",join(seeds,", "),"; raw relative tolerance=",RAW_TOL,"; score absolute tolerance=",SCORE_TOL)
        println(io,"Raw error denominator=max(1,abs(before),abs(after)); equal signed infinities match.")
        println(io,"No rank recovery threshold is treated as a correctness oracle; k sensitivity is descriptive.")
        println(io,"Mid-CDF oracle uses exposed normalization_scores; raw log-density unit checks remain unchanged.")
        source=join(read.(sort(readdir(joinpath(@__DIR__,"..","..","src");join=true)),String),"\n")
        println(io,"Source SHA256: ",bytes2hex(sha256(source)))
        println(io,"Baseline errors are expected only for documented beta ties or Gaussian singular/constant fits.")
    end
    println("Checks: ",length(checks),"; violations: ",count(x->x.status=="violation",checks),
            "; unexpected comparison errors: ",count(x->x.status=="unexpected_error",checks))
    println("Baseline unexpected errors: ",count(x->x.status=="unexpected_error",baselines))
    outdir
end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && ConsistencyStudy.main(isempty(ARGS) ?
    joinpath(@__DIR__,"results") : only(ARGS))
