using Platypus
using Test
using DataFrames

@testset "Platypus.jl" begin

    @testset "Entity and EntitySet" begin
        e1 = Entity(id=1, features=[0.2, 0.8, 0.5])
        @test e1.id == 1
        @test features(e1) == [0.2, 0.8, 0.5]
        @test label(e1) === nothing

        e2 = Entity(id=2, features=[0.4, 0.1, 0.9], label="rock")
        @test label(e2) == "rock"

        # An entity with no features is rejected
        @test_throws ArgumentError Entity(id=3, features=Float64[])

        es = EntitySet([e1, e2], [:a, :b, :c])
        @test n_entities(es) == 2
        @test n_features(es) == 3
        @test feature_names(es) == [:a, :b, :c]
        @test es.matrix[2, :] == [0.4, 0.1, 0.9]
        @test es.entity_ids == [1, 2]

        # Mismatched feature count is rejected
        bad = Entity(id=9, features=[1.0, 2.0])
        @test_throws DimensionMismatch EntitySet([e1, bad], [:a, :b, :c])

        # Matrix constructor
        es2 = EntitySet([0.1 0.2; 0.3 0.4; 0.5 0.6], [:x, :y])
        @test n_entities(es2) == 3
        @test es2.entity_ids == [1, 2, 3]
    end

    @testset "Data Loading" begin
        df = DataFrame(
            id    = [10, 20, 30],
            genre = ["jazz", "rock", "jazz"],
            danceability = [0.3, 0.7, 0.5],
            acousticness = [0.9, 0.2, 0.6],
        )
        es = load_entities(df; id_col=:id, label_col=:genre)
        @test n_entities(es) == 3
        @test n_features(es) == 2                       # genre/id excluded
        @test feature_names(es) == [:danceability, :acousticness]
        @test es.entity_ids == [10, 20, 30]
        @test es.labels[2] == "rock"

        # Explicit feature columns
        es2 = load_entities(df; feature_cols=[:danceability], id_col=:id)
        @test n_features(es2) == 1

        # Missing feature values are rejected
        df_missing = DataFrame(id=[1, 2], f=[0.5, missing])
        @test_throws ArgumentError load_entities(df_missing; id_col=:id)
    end

    @testset "build_pairs" begin
        es = EntitySet([0.1 0.2 0.9;
                        0.4 0.5 0.6;
                        0.7 0.8 0.3;
                        0.2 0.1 0.5], [:a, :b, :c])
        ps = build_pairs(es)
        @test n_pairs(ps) == 3                          # 3*(3-1)/2
        @test n_entities(ps) == 4

        # :rank transform puts every coordinate in (0, 1)
        for fp in ps
            @test all(0 .< fp.u .< 1)
            @test all(0 .< fp.v .< 1)
        end

        # Transforms are accepted; unknown ones are rejected
        @test n_pairs(build_pairs(es; transform=:zscore)) == 3
        @test n_pairs(build_pairs(es; transform=:none)) == 3
        @test_throws ArgumentError build_pairs(es; transform=:nonsense)

        # Need at least two features
        es1 = EntitySet(reshape([0.1, 0.2, 0.3], 3, 1), [:only])
        @test_throws ArgumentError build_pairs(es1)
    end

    @testset "EDA" begin
        es = EntitySet([0.1 0.2 0.9;
                        0.4 0.5 0.6;
                        0.7 0.8 0.3;
                        0.2 0.1 0.5], [:a, :b, :c])

        fs = feature_summary(es)
        @test nrow(fs) == 3
        @test "median" in names(fs)

        R = rank_matrix(es)
        @test size(R) == (4, 3)
        @test all(0 .< R .< 1)

        pc = pair_correlations(es)
        @test nrow(pc) == 3
        @test "correlation" in names(pc)

        ps = build_pairs(es)
        pr = pair_redundancy(ps)
        @test nrow(pr) == 3
        @test all(0 .<= pr.weight .<= 1)

        ne = nearest_entities(es, 1; k=2)
        @test nrow(ne) == 2
        @test issorted(ne.distance)

        ces = compute_entity_statistics(es)
        @test nrow(ces) == 4
        @test "extremeness" in names(ces)

        pt = pair_table(ps)
        @test nrow(pt) == 3
    end

    @testset "fit_typicality: every method" begin
        # 30 entities loosely clustered, plus features 2 and 3
        M = Float64[(i % 7) / 7 + (j - 1) * 0.05 for i in 1:30, j in 1:4]
        es = EntitySet(M, [:f1, :f2, :f3, :f4])
        ps = build_pairs(es)

        for method in TYPICALITY_METHODS
            r = fit_typicality(ps; method=method)
            @test r.method == method
            @test r.n_entities == 30
            @test r.n_pairs == 6
            @test size(r.pair_scores) == (30, 6)
            @test length(atypicality(r)) == 30
            @test all(0 .<= conventionality(r) .<= 1)
            @test all(atypicality(r) .≈ 1 .- conventionality(r))

            st = summarytable(r)
            @test nrow(st) == 30
            @test "atypicality" in names(st)

            sc = scoretable(r)
            @test nrow(sc) == 30 * 6
            @test "typicality" in names(sc)
        end
    end

    @testset "fit_typicality: flags an off-diagonal outlier" begin
        # 39 entities on the diagonal of a 3-feature space; the 40th
        # violates the correlation (high f1, low f2) and so should be
        # the most atypical entity under a pairwise density measure.
        rows = Vector{Float64}[]
        for i in 1:39
            x = i / 40
            push!(rows, [x, x, x])
        end
        push!(rows, [0.95, 0.05, 0.50])              # the outlier
        M = permutedims(reduce(hcat, rows))
        es = EntitySet(M, [:f1, :f2, :f3])
        ps = build_pairs(es)

        r = fit_typicality(ps; method=:depth)
        @test argmax(atypicality(r)) == 40           # outlier is row 40
        @test atypicality(r)[40] > atypicality(r)[20]
    end

    @testset "fit_typicality: ReferenceWeights and errors" begin
        es = EntitySet([0.1 0.2 0.3;
                        0.4 0.5 0.6;
                        0.7 0.8 0.9;
                        0.2 0.9 0.1], [:a, :b, :c])
        ps = build_pairs(es)

        rw = ReferenceWeights([1.0, 1.0, 1.0, 0.0])
        r = fit_typicality(ps; method=:kde, weights=rw)
        @test r.n_entities == 4

        # Uniform-weight constructor
        @test length(ReferenceWeights(4)) == 4
        @test_throws ArgumentError ReferenceWeights([-1.0, 1.0])
        @test_throws ArgumentError ReferenceWeights([0.0, 0.0])

        # Bad arguments to fit_typicality
        @test_throws ArgumentError fit_typicality(ps; method=:bogus)
        @test_throws DimensionMismatch fit_typicality(
            ps; method=:kde, weights=ReferenceWeights([1.0, 1.0]))
    end

    @testset "Estimation stubs raise" begin
        es = EntitySet([0.1 0.2; 0.3 0.4; 0.5 0.6], [:a, :b])
        ps = build_pairs(es)
        r  = fit_typicality(ps; method=:knn)
        @test_throws ErrorException fit_graph_typicality(es)
        @test_throws ErrorException fit_outcome(r, [1.0, 2.0, 3.0])
        @test_throws ErrorException fit_copula(ps)
    end

    @testset "Visualization stubs raise" begin
        es = EntitySet([0.1 0.2 0.7;
                        0.3 0.4 0.5;
                        0.5 0.6 0.2;
                        0.9 0.1 0.4], [:a, :b, :c])
        ps = build_pairs(es)
        r  = fit_typicality(ps; method=:lof)

        @test_throws ErrorException plot_feature_distributions(es)
        @test_throws ErrorException plot_pair_cloud(ps, 1)
        @test_throws ErrorException plot_correlation_heatmap(es)
        @test_throws ErrorException plot_pair_redundancy(ps)
        @test_throws ErrorException plot_typicality_profile(r, 1)
        @test_throws ErrorException plot_atypicality_distribution(r)
        @test_throws ErrorException plot_novelty_conventionality(r)
        @test_throws ErrorException plot_outcome_fit(nothing)
        @test_throws ErrorException plot_pair_heatmap(r)
    end
end
