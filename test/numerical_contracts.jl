@testset "Structural and numerical provenance contracts" begin
    X = [1.0 8.0 3.0; 2.0 3.0 8.0; 3.0 10.0 5.0; 4.0 1.0 11.0;
         5.0 7.0 2.0; 6.0 12.0 7.0; 7.0 2.0 10.0; 8.0 9.0 1.0;
         9.0 4.0 12.0; 10.0 11.0 4.0; 11.0 6.0 9.0; 12.0 5.0 6.0]
    es = EntitySet(X, [:a, :b, :c])
    fitters = Function[ps -> fit_typicality(ps; method=method, k=method in (:knn, :lof) ? 2 : 0)
               for method in (:knn, :lof, :kde, :depth)]
    append!(fitters, [ps -> fit_copula(ps; family=family) for family in (:gaussian, :beta)])

    @testset "Complete unique feature-pair coverage" begin
        for fit in fitters
            baseline = fit(build_pairs(es))
            duplicate = build_pairs(es)
            duplicate.pairs[2] = duplicate.pairs[1]
            @test_throws ArgumentError fit(duplicate)
            original = build_pairs(es)
            missing = PairSet(original.pairs[1:2], original.feature_names,
                original.entity_ids, original.labels, original.n_entities, 2,
                original.transform, original.source_matrix)
            @test_throws ArgumentError fit(missing)
            reordered = build_pairs(es)
            reverse!(reordered.pairs)
            result = fit(reordered)
            @test result.raw_scores == reverse(baseline.raw_scores; dims=2)
            @test result.normalization_scores == reverse(baseline.normalization_scores; dims=2)
            @test result.pair_scores == reverse(baseline.pair_scores; dims=2)
            @test result.conventionality == baseline.conventionality
            @test result.novelty == baseline.novelty
            @test result.dispersion == baseline.dispersion
        end
    end

    @testset "Inspectable normalization values and exports" begin
        ps = build_pairs(EntitySet(X[1:4, 1:2], [:a, :b]); transform=:none)
        raw = fill(-300.0, 4, 1)
        keys = reshape([2.0, 1.0, 3.0, 2.0], :, 1)
        source_weights = [1.0, 2.0, 3.0, 4.0]
        result = Platypus._summarize(ps, raw, :kde, 0, 0.1,
            source_weights ./ 10, (;); score_weights=source_weights,
            normalization_scores=keys, normalization_basis=:scale_free_kernel_order)
        @test result.raw_scores == raw
        @test result.normalization_scores == keys
        @test result.pair_scores[:, 1] == [0.45, 0.1, 0.85, 0.45]
        @test scoretable(result; raw=true).typicality == vec(raw)
        @test scoretable(result; normalization=true).typicality == vec(keys)
        @test scoretable(result).typicality == vec(result.pair_scores)
        @test_throws ArgumentError scoretable(result; raw=true, normalization=true)
        @test analysis_metadata(result).normalization_basis === :scale_free_kernel_order
        keys[1] = -1
        @test result.normalization_scores[1] == 2
        @test_throws DimensionMismatch Platypus._summarize(ps, raw, :kde,
            0, 0.1, source_weights ./ 10, (;); normalization_scores=zeros(3, 1))
        for method in (:knn, :lof, :depth)
            fitted = fit_typicality(build_pairs(es); method=method)
            @test fitted.normalization_scores == fitted.raw_scores
            @test fitted.normalization_scores !== fitted.raw_scores
            @test analysis_metadata(fitted).normalization_basis === :raw_score
        end
        for family in (:gaussian, :beta)
            fitted = fit_copula(build_pairs(es); family=family)
            @test fitted.normalization_scores == fitted.raw_scores
            @test analysis_metadata(fitted).normalization_basis === :raw_score
        end
        rank_kde = fit_typicality(build_pairs(es); method=:kde)
        @test rank_kde.normalization_scores == rank_kde.raw_scores
        @test analysis_metadata(rank_kde).normalization_basis === :raw_score

        # A common physical-unit shift collapses two distinct Float64 log
        # densities. Their inspectable ordering must still determine the CDF.
        close = [-2.0 0; 0 0; 1e-14 0; 1 0]
        real_kde(M; weights=nothing) = fit_typicality(
            build_pairs(EntitySet(M, [:x, :y]); transform=:none);
            method=:kde, weights=weights)
        baseline = real_kde(close)
        changed_units = real_kde(ldexp.(close, 200))
        @test baseline.raw_scores[2] < baseline.raw_scores[3]
        @test changed_units.raw_scores[2] == changed_units.raw_scores[3]
        @test changed_units.normalization_scores[2] < changed_units.normalization_scores[3]
        @test changed_units.pair_scores == baseline.pair_scores
        @test changed_units.normalization_scores == baseline.normalization_scores
        @test analysis_metadata(changed_units).normalization_basis === :scale_free_kernel_order
        @test changed_units.raw_scores ≈ baseline.raw_scores .- 200log(2.0)
        # Ordinal label magnitudes need not remain identical when an additional
        # query creates a new ordering level, but existing percentiles must.
        extra_query = real_kde(vcat(close, [10.0 0.0]);
            weights=ReferenceWeights([1.0, 1.0, 1.0, 1.0, 0.0]))
        @test extra_query.pair_scores[1:4, :] == baseline.pair_scores
    end

    @testset "Transform provenance is coherent" begin
        ps = build_pairs(es)
        known_without_source = PairSet(ps.pairs, ps.feature_names, ps.entity_ids,
            ps.labels, ps.n_entities, ps.n_pairs, :rank, nothing)
        unknown_with_source = PairSet(ps.pairs, ps.feature_names, ps.entity_ids,
            ps.labels, ps.n_entities, ps.n_pairs, :unknown, ps.source_matrix)
        for fit in fitters
            @test_throws ArgumentError fit(known_without_source)
            @test_throws ArgumentError fit(unknown_with_source)
        end
        legacy = PairSet(ps.pairs, ps.feature_names, ps.entity_ids, ps.labels,
            ps.n_entities, ps.n_pairs)
        @test legacy.transform === :unknown
        @test all(isfinite, fit_typicality(legacy; method=:depth).raw_scores)
        @test_throws ArgumentError fit_typicality(legacy; method=:depth,
            weights=ReferenceWeights(collect(1.0:12.0)))
    end

    @testset "Exact source mass survives normalization" begin
        source = [0.8594636272788749, 0.9564597692887471, 0.9916107558439862,
                  0.41124825032947365, 0.47321488616625174, 0.045419171599359265]
        normalize(w) = Platypus._validated_weights(ReferenceWeights(w), length(w))
        baseline = normalize(source)
        permutation = [6, 2, 5, 1, 4, 3]
        @test normalize(source[permutation]) == baseline[permutation]
        @test normalize(vcat(source, zeros(257)))[1:6] == baseline
        @test all(iszero, normalize(vcat(source, zeros(257)))[7:end])
        @test normalize([1e308, 1e308]) == [0.5, 0.5]
        @test normalize([4.0, 2.0, 1.0]) == [4/7, 2/7, 1/7]
        @test normalize(ldexp.([4.0, 2.0, 1.0], 900)) == normalize([4.0, 2.0, 1.0])
        tiny = nextfloat(0.0)
        @test normalize([tiny, 2tiny, 3tiny]) == [1/6, 1/3, 1/2]
        @test_throws ArgumentError normalize([tiny, floatmax(Float64)])
        @test_throws ArgumentError Platypus._validated_weights(ReferenceWeights([1.0]), 1; min_positive=0)
    end

    @testset "Exported mass is bound to entity identifiers" begin
        ids = [30, 10, 40, 20]
        small = EntitySet([0.0 3.0; 1.0 0.0; 2.0 2.0; 4.0 1.0], [:x, :y]; ids=ids)
        rw = ReferenceWeights([1.0, 2.0, 3.0, 4.0])
        result = fit_typicality(build_pairs(small); method=:depth, weights=rw)
        metadata = analysis_metadata(result)
        @test metadata.entity_ids == ids
        @test metadata.source_reference_weights == rw.weights
        exported = sort(summarytable(result), :entity_id)
        masses = Dict(zip(metadata.entity_ids, metadata.reference_weights))
        @test [masses[id] for id in exported.entity_id] == [0.2, 0.4, 0.1, 0.3]
        metadata.entity_ids[1] = -1
        metadata.source_reference_weights[1] = -1
        @test analysis_metadata(result).entity_ids == ids
        @test analysis_metadata(result).source_reference_weights == rw.weights
        outcome = fit_outcome(result, [2.0, 1.0, 4.0, 3.0]; predictors=())
        @test analysis_metadata(outcome).measurement.entity_ids == ids
        graph = fit_graph_typicality(small; k=2)
        @test analysis_metadata(graph).entity_ids == ids
        for method in (:knn, :lof)
            backend = analysis_metadata(fit_typicality(build_pairs(small); method=method, k=2, backend=:kdtree))
            @test backend.requested_backend === :kdtree
            @test backend.effective_backend === :brute
            @test backend.backend === backend.effective_backend
        end
    end

    @testset "Entity selection retains ordinary Boolean-mask semantics" begin
        selected = es[isodd.(1:length(es))]
        @test [e.id for e in selected] == [1, 3, 5, 7, 9, 11]
        selected[1].features[1] = -100
        @test es.matrix[1, 1] == 1
        @test [e.id for e in es[[3, 1, 3]]] == [3, 1, 3]
        @test_throws BoundsError es[trues(length(es) - 1)]
    end

    @testset "Finite square-root results and ambient precision" begin
        standardized = Platypus._weighted_zscore([0.0, 1.0, 1e200], [1.0, 1.0, 0.0])
        @test all(isfinite, standardized)
        @test standardized[1:2] == [-sqrt(0.5), sqrt(0.5)]
        @test standardized[3] ≈ sqrt(2.0) * 1e200
        weights = ReferenceWeights(collect(1.0:12.0))
        baseline_rank = Platypus._transform_columns(X, :rank, weights.weights)
        baseline_zscore = Platypus._transform_columns(X, :zscore, weights.weights)
        baseline = fit_typicality(build_pairs(es); method=:depth, weights=weights)
        # BigFloat precision is deliberately varied only in this serial test.
        # No fit implementation may change process-global MPFR settings.
        setprecision(BigFloat, 16) do
            @test Platypus._transform_columns(X, :rank, weights.weights) == baseline_rank
            @test Platypus._transform_columns(X, :zscore, weights.weights) == baseline_zscore
            low_precision = fit_typicality(build_pairs(es); method=:depth, weights=weights)
            @test low_precision.raw_scores == baseline.raw_scores
            @test low_precision.pair_scores == baseline.pair_scores
            @test analysis_metadata(low_precision).reference_weights == analysis_metadata(baseline).reference_weights
        end
    end
end
