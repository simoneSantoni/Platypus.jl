# Run from an instantiated package environment:
# julia --project validation/experiments.jl [output_directory]
# Small, fixed-design numerical comparison; not a general anomaly benchmark.
using Platypus, Random, Statistics, StatsBase, CSV, DataFrames, Printf

const SEEDS = (101, 202, 303)
const METHODS = (:lof, :knn, :kde, :depth, :gaussian_copula, :beta_copula)

function fit_method(M, method; transform=:rank, backend=:kdtree, provenance=(;))
    ps = build_pairs(EntitySet(M, Symbol.("f" .* string.(1:size(M,2)))); transform)
    method === :gaussian_copula && return fit_copula(ps; family=:gaussian, provenance)
    method === :beta_copula && return fit_copula(ps; family=:beta, provenance)
    fit_typicality(ps; method, backend, provenance)
end

# AUROC gives ties half credit. AP groups tied score thresholds, avoiding an
# arbitrary row-order tie break. Precision at the known anomaly count averages
# uniformly across ties at the selection boundary.
function recovery(scores, truth)
    pos, neg = findall(truth), findall(.!truth)
    auc = sum(scores[i] > scores[j] ? 1.0 : scores[i] == scores[j] ? 0.5 : 0.0
              for i in pos, j in neg) / (length(pos)*length(neg))
    order = sortperm(scores; rev=true)
    tp = 0
    ap = 0.0
    i = 1
    while i <= length(order)
        j = i
        while j < length(order) && scores[order[j+1]] == scores[order[i]]
            j += 1
        end
        block_tp = count(truth[order[i:j]])
        tp += block_tp
        ap += block_tp/length(pos) * tp/j
        i = j+1
    end
    cutoff = scores[order[length(pos)]]
    above, tied = scores .> cutoff, scores .== cutoff
    expected_tp = count(above .& truth) +
        (length(pos)-count(above))*count(tied .& truth)/count(tied)
    (auc=auc, average_precision=ap, precision_at_m=expected_tp/length(pos))
end

function sample_case(scenario, seed; n=160, m=12)
    rng = MersenneTwister(seed)
    truth = [falses(n-m); trues(m)]
    if scenario === :independent
        # Labels contain no signal, assigned independently of all coordinates.
        M = rand(rng,n,4)
    elseif scenario === :dependence || scenario === :irrelevant_features
        z = randn(rng,n)
        M = z .+ 0.12 .* randn(rng,n,4)
        # Reverse two of four features: four of six pairs violate dependence,
        # so the median, not only a selected pair, can register the violation.
        M[end-m+1:end,[2,4]] .*= -1
        scenario === :irrelevant_features && (M = hcat(M,randn(rng,n,4)))
    elseif scenario === :unequal_clusters
        # Unequal cluster mass and dispersion; bridge points are planted
        # anomalies. Density and global centrality need not agree here.
        n1 = 108
        M = vcat(-2 .+ 0.08 .* randn(rng,n1,4),
                  2 .+ 0.55 .* randn(rng,n-m-n1,4),
                       0.12 .* randn(rng,m,4))
    elseif scenario === :higher_order
        # U,V,(U+V) mod 1 have independent pairwise margins but a joint
        # constraint. Replacing the third coordinate violates only that
        # higher-order constraint; bivariate densities retain no contrast.
        U,V = rand(rng,n),rand(rng,n)
        W = mod.(U .+ V,1)
        W[end-m+1:end] = rand(rng,m)
        M = hcat(U,V,W)
    else
        error("unknown scenario")
    end
    M,truth
end

function main(outdir=joinpath(@__DIR__,"results"))
    mkpath(outdir)
    rows = NamedTuple[]
    stability = NamedTuple[]
    # Warm each implementation; wall times below remain indicative only.
    warm = rand(MersenneTwister(9),30,4)
    for method in METHODS
        fit_method(warm,method)
    end
    for scenario in (:independent,:dependence,:unequal_clusters,:irrelevant_features,:higher_order), seed in SEEDS
        M,truth = sample_case(scenario,seed)
        for method in METHODS
            timed = @timed fit_method(M,method;
                provenance=(dataset="synthetic-v1",cohort=string(scenario),seed=seed))
            r = timed.value
            for (summary,values) in ((:atypicality,r.atypicality),(:tail_unusualness,1 .- r.novelty))
                metrics = recovery(values,truth)
                push!(rows,merge((scenario=string(scenario),seed=seed,method=string(method),
                    summary=string(summary),n=size(M,1),features=size(M,2),
                    prevalence=mean(truth)),metrics,(seconds=timed.time,bytes=timed.bytes,)))
            end
        end
    end
    for seed in SEEDS
        base,_ = sample_case(:dependence,seed)
        extended,_ = sample_case(:irrelevant_features,seed)
        permutation = randperm(MersenneTwister(seed+1),size(base,1))
        for method in METHODS
            original = fit_method(base,method)
            permuted = fit_method(base[permutation,:],method)
            extra = fit_method(extended,method)
            aligned = permuted.atypicality[invperm(permutation)]
            brute_error = method in (:lof,:knn) ?
                maximum(abs.(original.raw_scores .- fit_method(base,method;backend=:brute).raw_scores)) : missing
            push!(stability,(seed=seed,method=string(method),
                row_permutation_max_error=maximum(abs.(original.atypicality .- aligned)),
                row_permutation_spearman=corspearman(original.atypicality,aligned),
                irrelevant_feature_spearman=corspearman(original.atypicality,extra.atypicality),
                exact_backend_max_raw_error=brute_error))
        end
    end
    results = DataFrame(rows)
    CSV.write(joinpath(outdir,"recovery.csv"),results)
    CSV.write(joinpath(outdir,"stability.csv"),DataFrame(stability))
    summary = combine(groupby(results,[:scenario,:method,:summary]),
        :auc=>mean=>:mean_auc,:auc=>minimum=>:min_auc,:auc=>maximum=>:max_auc,
        :average_precision=>mean=>:mean_average_precision,
        :precision_at_m=>mean=>:mean_precision_at_m)
    CSV.write(joinpath(outdir,"summary.csv"),summary)
    open(joinpath(outdir,"environment.txt"),"w") do io
        println(io,"Julia ",VERSION,"; Platypus ",pkgversion(Platypus))
        println(io,"Seeds: ",join(SEEDS,", "),"; Threads: ",Threads.nthreads())
        println(io,"CPU: ",Sys.CPU_NAME,"; architecture: ",Sys.ARCH)
        println(io,"n=160; anomalies=12; rank transform; default k and bandwidth; no tuning")
        println(io,"Negative controls and higher-order cases are retained regardless of performance.")
        println(io,"Scores are evaluated within each endogenous collection; no held-out prediction claim.")
        println(io,"Times exclude one warm-up per method but are not a controlled HPC benchmark.")
    end
    show(stdout,MIME("text/plain"),summary[summary.summary .== "atypicality",:])
    println("\nResults written to ",abspath(outdir))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(isempty(ARGS) ? joinpath(@__DIR__,"results") : first(ARGS))
