using Platypus, Test, DataFrames

# Run from the repository root with --project=.
println("Sociology acceptance controls")
println("Julia ", VERSION)
println("Scope: identity/reference binding and preserved operational contracts")

@testset "Sociology final identity and reference controls" begin
    X = [0.0 3.0; 1.0 0.0; 2.0 2.0; 4.0 1.0; 1.5 1.5]
    ids = [30, 10, 40, 20, 99]
    w = [1.0, 2.0, 3.0, 4.0, 0.0]
    es = EntitySet(X, [:x, :y]; ids=ids, labels=["a", "b", "a", "b", "query"])
    ps = build_pairs(es)
    fits = [p -> fit_typicality(p; method=:kde, weights=ReferenceWeights(w)),
            p -> fit_typicality(p; method=:depth, weights=ReferenceWeights(w)),
            p -> fit_copula(p; weights=ReferenceWeights(w))]
    for fit in fits
        result = fit(ps)
        metadata = analysis_metadata(result)
        exported = sort(summarytable(result), :entity_id)
        # A join does not promise to retain the left table's order. Sort its
        # output explicitly before checking the independently known ID map.
        joined = sort(leftjoin(exported,
            DataFrame(entity_id=metadata.entity_ids,
                      original_mass=metadata.source_reference_weights,
                      normalized_mass=metadata.reference_weights);
            on=:entity_id), :entity_id)
        @test joined.original_mass == [2.0, 4.0, 1.0, 3.0, 0.0]
        @test joined.normalized_mass == [0.2, 0.4, 0.1, 0.3, 0.0]
        @test joined.entity_id[end] == 99 && joined.normalized_mass[end] == 0
        @test result.atypicality == 1 .- result.conventionality
        @test all(result.novelty .<= result.conventionality)
        @test result.pair_names == [(:x, :y)]
        metadata.entity_ids[1] = -1
        metadata.source_reference_weights[1] = -1
        @test analysis_metadata(result).entity_ids == ids
        @test analysis_metadata(result).source_reference_weights == w
        altered = build_pairs(EntitySet(X, [:x, :y];
            ids=reverse(ids), labels=fill("relabelled", 5)))
        @test fit(altered).raw_scores == result.raw_scores
        @test fit(altered).pair_scores == result.pair_scores
    end

    # Positive-mass splitting changes the row-level leave-one-out estimand.
    # Fix bandwidth to isolate kernel exclusion from bandwidth adaptation.
    A = [0.0 0.0; 1.0 1.0; 2.0 0.0]
    a = fit_typicality(build_pairs(EntitySet(A, [:x, :y]); transform=:none);
        method=:kde, weights=ReferenceWeights([2.0, 1.0, 1.0]), bandwidth=1.0)
    b = fit_typicality(build_pairs(EntitySet(A[[1, 1, 2, 3], :], [:x, :y]);
        transform=:none); method=:kde, bandwidth=1.0)
    K(x, y) = exp(-(x*x + y*y)/2)/(2pi)
    @test a.raw_scores[1, 1] ≈ log((K(1, 1) + K(2, 0))/2)
    @test b.raw_scores[1, 1] ≈ log((K(0, 0) + K(1, 1) + K(2, 0))/3)
    @test a.raw_scores[1, 1] != b.raw_scores[1, 1]
end
