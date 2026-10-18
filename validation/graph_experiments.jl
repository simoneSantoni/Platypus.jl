# Supplemental full-feature connectivity comparison. Reuses the exact same
# generators and metrics without rerunning the pairwise/copula experiments.
# julia --project validation/graph_experiments.jl [output_directory]
include(joinpath(@__DIR__,"experiments.jl"))

function fit_graph_case(M; provenance=(;))
    es = EntitySet(M,Symbol.("f" .* string.(1:size(M,2))))
    fit_graph_typicality(es;provenance)
end

function graph_main(outdir=joinpath(@__DIR__,"results"))
    mkpath(outdir)
    rows = NamedTuple[]
    stability = NamedTuple[]
    fit_graph_case(rand(MersenneTwister(9),30,4))
    for scenario in (:independent,:dependence,:unequal_clusters,:irrelevant_features,:higher_order), seed in SEEDS
        M,truth = sample_case(scenario,seed)
        timed = @timed fit_graph_case(M;
            provenance=(dataset="synthetic-v1",cohort=string(scenario),seed=seed))
        graph = timed.value
        metadata = analysis_metadata(graph)
        metrics = recovery(graph.atypicality,truth)
        push!(rows,merge((scenario=string(scenario),seed=seed,method="self_tuning_graph",
            target="full_feature_connectivity",summary="atypicality",n=size(M,1),
            features=size(M,2),prevalence=mean(truth)),metrics,
            (k=metadata.k,n_edges=metadata.n_edges,n_components=metadata.n_components,
             seconds=timed.time,bytes=timed.bytes)))
    end
    for seed in SEEDS
        M,_ = sample_case(:dependence,seed)
        extra,_ = sample_case(:irrelevant_features,seed)
        permutation = randperm(MersenneTwister(seed+1),size(M,1))
        graph = fit_graph_case(M)
        permuted = fit_graph_case(M[permutation,:])
        extended = fit_graph_case(extra)
        aligned = permuted.atypicality[invperm(permutation)]
        push!(stability,(seed=seed,method="self_tuning_graph",
            target="full_feature_connectivity",
            row_permutation_max_error=maximum(abs.(graph.atypicality .- aligned)),
            row_permutation_spearman=corspearman(graph.atypicality,aligned),
            irrelevant_feature_spearman=corspearman(graph.atypicality,extended.atypicality)))
    end
    results = DataFrame(rows)
    CSV.write(joinpath(outdir,"graph_recovery.csv"),results)
    CSV.write(joinpath(outdir,"graph_stability.csv"),DataFrame(stability))
    summary = combine(groupby(results,[:scenario,:method,:target,:summary]),
        :auc=>mean=>:mean_auc,:auc=>minimum=>:min_auc,:auc=>maximum=>:max_auc,
        :average_precision=>mean=>:mean_average_precision,
        :precision_at_m=>mean=>:mean_precision_at_m)
    CSV.write(joinpath(outdir,"graph_summary.csv"),summary)
    open(joinpath(outdir,"graph_environment.txt"),"w") do io
        println(io,"Julia ",VERSION,"; Platypus ",pkgversion(Platypus))
        println(io,"Seeds: ",join(SEEDS,", "),"; Threads: ",Threads.nthreads())
        println(io,"CPU: ",Sys.CPU_NAME,"; architecture: ",Sys.ARCH)
        println(io,"n=160; anomalies=12; rank transform; default k; exact kdtree; no tuning")
        println(io,"Target: full-feature self-tuning graph connectivity; no novelty summary.")
        println(io,"Same generators and metrics as experiments.jl; original results were not rerun.")
        println(io,"Endogenous in-sample comparison; timings indicative after a four-feature warm-up.")
        println(io,"First fits at three/eight features can include dimension-specific compilation.")
    end
    show(stdout,MIME("text/plain"),summary)
    println("\nGraph supplement written to ",abspath(outdir))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && graph_main(isempty(ARGS) ? joinpath(@__DIR__,"results") : first(ARGS))
