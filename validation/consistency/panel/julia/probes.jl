# Review-only reproductions: preserve current behavior, including violations.
using Platypus
using Random
using Statistics
using CSV
using DataFrames
using SHA

const OUT = @__DIR__
const REPO = normpath(joinpath(OUT,"../../../.."))

function environment_record()
    open(joinpath(OUT,"environment.txt"),"w") do io
        println(io,"Julia ",VERSION)
        println(io,"Threads ",Threads.nthreads())
        for path in sort(readdir(joinpath(REPO,"src");join=true))
            endswith(path,".jl") || continue
            println(io,relpath(path,REPO)," sha256=",bytes2hex(sha256(read(path))))
        end
        println(io,"Project.toml sha256=",bytes2hex(sha256(read(joinpath(REPO,"Project.toml")))))
    end
end

fit_for(ps,method) = method in (:gaussian,:beta) ? fit_copula(ps;family=method) :
    fit_typicality(ps;method=method)

function all_pairs_probe()
    rng = MersenneTwister(20261018)
    X = Float64.(hcat(collect(1:12),randperm(rng,12),randperm(rng,12)))
    names = [:a,:b,:c]
    CSV.write(joinpath(OUT,"pair_fixture.csv"),DataFrame(entity_id=1:12,a=X[:,1],b=X[:,2],c=X[:,3]))
    rows = NamedTuple[]
    scores = NamedTuple[]
    for method in (:knn,:lof,:kde,:depth,:gaussian,:beta)
        baseline_ps = build_pairs(EntitySet(X,names))
        baseline = fit_for(baseline_ps,method)
        duplicate = build_pairs(EntitySet(X,names))
        duplicate.pairs[2] = duplicate.pairs[1]
        validation_accepted = Platypus._validate_pairset(duplicate) === duplicate
        corrupted = fit_for(duplicate,method)
        # No numerical tolerance: surviving and repeated pair scores are exact
        # copies of baseline columns. Only the selected collection changed.
        @assert corrupted.raw_scores == baseline.raw_scores[:,[1,1,3]]
        @assert corrupted.pair_scores == baseline.pair_scores[:,[1,1,3]]
        expected = [median(baseline.pair_scores[i,[1,1,3]]) for i in 1:12]
        @assert corrupted.conventionality == expected
        changed = count(!iszero,corrupted.atypicality .- baseline.atypicality)
        delta = maximum(abs.(corrupted.atypicality .- baseline.atypicality))
        missing = PairSet(baseline_ps.pairs[[1,3]],copy(names),copy(baseline_ps.entity_ids),
            copy(baseline_ps.labels),12,2,baseline_ps.transform,copy(baseline_ps.source_matrix))
        missing_accepted = Platypus._validate_pairset(missing) === missing
        incomplete = fit_for(missing,method)
        reordered = build_pairs(EntitySet(X,names))
        reverse!(reordered.pairs)
        reordered_fit = fit_for(reordered,method)
        @assert reordered_fit.raw_scores[:,end:-1:1] == baseline.raw_scores
        @assert reordered_fit.atypicality == baseline.atypicality
        push!(rows,(method=string(method),duplicate_accepted=validation_accepted,
            unique_pairs=length(unique(corrupted.pair_names)),reported_pairs=corrupted.n_pairs,
            changed_entities=changed,max_atypicality_change=delta,
            max_conventionality_change=maximum(abs.(corrupted.conventionality.-baseline.conventionality)),
            max_novelty_change=maximum(abs.(corrupted.novelty.-baseline.novelty)),
            max_dispersion_change=maximum(abs.(corrupted.dispersion.-baseline.dispersion)),
            missing_pair_accepted=missing_accepted,missing_pair_reported_count=incomplete.n_pairs,
            missing_pair_max_atypicality_change=maximum(abs.(incomplete.atypicality.-baseline.atypicality)),
            pair_order_control_passed=true))
        for i in 1:12
            push!(scores,(method=string(method),entity_id=i,
                baseline_atypicality=baseline.atypicality[i],duplicate_atypicality=corrupted.atypicality[i],
                incomplete_atypicality=incomplete.atypicality[i],
                baseline_ab=baseline.pair_scores[i,1],baseline_ac=baseline.pair_scores[i,2],
                baseline_bc=baseline.pair_scores[i,3]))
        end
    end
    CSV.write(joinpath(OUT,"pair_validation.csv"),DataFrame(rows))
    CSV.write(joinpath(OUT,"pair_scores.csv"),DataFrame(scores))
    return rows
end

function snapshot_controls()
    X = Float64[0 2 3;1 3 2;2 1 0;3 0 1]
    es = EntitySet(X,[:a,:b,:c])
    ps = build_pairs(es)
    provenance = (dataset="review-fixture",notes=["original"])
    result = fit_typicality(ps;method=:knn,provenance=provenance)
    baseline = deepcopy(result)
    X[1,1] = 123
    @assert es.matrix[1,1] == 0
    es.matrix[1,1] = 99
    es.entity_ids[1] = 500
    @assert ps.source_matrix[1,1] == 0
    @assert ps.entity_ids[1] == 1
    snapshot = es[2]
    snapshot.features[1] = 900
    @assert es.matrix[2,1] == 1
    ps.entity_ids[1] = 700
    provenance.notes[1] = "edited"
    @assert result.entity_ids == baseline.entity_ids
    @assert analysis_metadata(result).provenance.notes == ["original"]
    table = summarytable(result)
    table.atypicality[1] = -99
    @assert result.atypicality == baseline.atypicality
    metadata = analysis_metadata(result)
    metadata.reference_weights[1] = 0
    @assert result.configuration.reference_weights == baseline.configuration.reference_weights
    return (input_matrix_copied=true,entity_accessor_detached=true,
        pair_source_detached=true,result_metadata_detached=true,
        summary_table_detached=true,metadata_accessor_detached=true)
end

function provenance_probe()
    X = Float64[1 2;2 4;3 1;4 3]
    ps = build_pairs(EntitySet(X,[:a,:b]))
    malformed = PairSet(ps.pairs,copy(ps.feature_names),copy(ps.entity_ids),copy(ps.labels),
        ps.n_entities,ps.n_pairs,:rank,nothing)
    accepted = Platypus._validate_pairset(malformed) === malformed
    caught = try
        fit_copula(malformed)
        "no error"
    catch error
        sprint(showerror,error)
    end
    legacy = PairSet(ps.pairs,copy(ps.feature_names),copy(ps.entity_ids),copy(ps.labels),
        ps.n_entities,ps.n_pairs)
    legacy_fit = fit_typicality(legacy;method=:knn)
    @assert legacy_fit.pair_scores == fit_typicality(ps;method=:knn).pair_scores
    rejected = try
        fit_typicality(legacy;method=:depth,weights=ReferenceWeights([1,2,3,4]))
        false
    catch error
        error isa ArgumentError
    end
    @assert rejected
    return (malformed_state_accepted=accepted,copula_error=caught,
        legacy_uniform_knn_control=true,legacy_nonuniform_rejected=true)
end

function backend_probe()
    X = [1e-100 0.;0. 1e-100;1. 2.;2. 1.;1e100 1.;1. 1e100]
    es = EntitySet(X,[:a,:b])
    ps = build_pairs(es;transform=:none)
    search = Platypus._neighbor_setup(ps[1].u,ps[1].v,2,:kdtree)
    fallback = search.tree === nothing
    fit = fit_typicality(ps;method=:knn,k=2,backend=:kdtree)
    brute = fit_typicality(ps;method=:knn,k=2,backend=:brute)
    @assert fit.raw_scores == brute.raw_scores
    @assert fit.pair_scores == brute.pair_scores
    graph = fit_graph_typicality(es;k=2,transform=:none,backend=:kdtree)
    ordinary = Platypus._neighbor_setup([0.,1.,2.],[2.,0.,1.],1,:kdtree)
    @assert ordinary.tree !== nothing
    return (requested="kdtree",effective_fallback=fallback ? "brute" : "kdtree",
        recorded_pair_backend=string(analysis_metadata(fit).backend),
        recorded_graph_backend=string(analysis_metadata(graph).backend),
        exact_agreement_with_brute=true,ordinary_tree_control=true)
end

function metadata_binding_probe()
    X = Float64[0 2;1 3;2 1;3 0]
    es = EntitySet(X,[:a,:b];ids=[30,10,40,20])
    result = fit_typicality(build_pairs(es);method=:depth,weights=ReferenceWeights([1,2,3,4]))
    metadata = analysis_metadata(result)
    table = sort(summarytable(result),:entity_id)
    original_binding = Dict(zip(result.entity_ids,metadata.reference_weights))
    expected = [original_binding[id] for id in table.entity_id]
    positional = metadata.reference_weights
    mismatch = count(!iszero,expected .- positional)
    @assert mismatch > 0
    CSV.write(joinpath(OUT,"metadata_binding.csv"),DataFrame(entity_id=table.entity_id,
        expected_reference_weight=expected,sidecar_position_weight=positional))
    outcome = fit_outcome(result,Float64[1,2,3,4];predictors=())
    outcome_metadata = analysis_metadata(outcome)
    graph = fit_graph_typicality(es;k=2)
    return (metadata_has_entity_ids=haskey(metadata,:entity_ids),
        sorted_table_ids=table.entity_id,original_ids=copy(result.entity_ids),
        incorrect_positional_matches=mismatch,max_weight_mapping_error=maximum(abs.(expected.-positional)),
        outcome_nested_metadata_has_entity_ids=haskey(outcome_metadata.measurement,:entity_ids),
        outcome_object_preserves_ids=outcome.entity_ids == result.entity_ids,
        graph_only_uniform_weights=all(==(first(graph.configuration.reference_weights)),graph.configuration.reference_weights))
end

function main()
    environment_record()
    pairs = all_pairs_probe()
    controls = snapshot_controls()
    provenance = provenance_probe()
    backend = backend_probe()
    metadata = metadata_binding_probe()
    open(joinpath(OUT,"summary.txt"),"w") do io
        println(io,"PAIR COLLECTION VALIDATION")
        for row in pairs
            println(io,row)
        end
        println(io,"SNAPSHOT CONTROLS: ",controls)
        println(io,"PROVENANCE CONTRACT: ",provenance)
        println(io,"BACKEND METADATA: ",backend)
        println(io,"SORTED EXPORT METADATA: ",metadata)
    end
    println(read(joinpath(OUT,"summary.txt"),String))
end
main()
