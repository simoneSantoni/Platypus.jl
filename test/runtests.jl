using Platypus
using Test
using DataFrames
using Random
using LinearAlgebra
using RecipesBase

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
        @test all(0 .<= pr.abs_correlation .<= 1)

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

    @testset "Additional implemented APIs" begin
        rng = MersenneTwister(77)
        es = EntitySet(rand(rng, 40, 4), [:a,:b,:c,:d])
        ps = build_pairs(es)
        r = fit_typicality(ps)
        @test fit_copula(ps).method == :gaussian_copula
        @test fit_copula(ps;family=:beta).method == :beta_copula
        y = 2 .+ 3 .* r.conventionality .- 4 .* r.novelty .+ 0.5 .* r.dispersion
        out = fit_outcome(r,y)
        @test out.coefficients ≈ [2,3,-4,0.5]
        @test out.fitted ≈ y
        @test out.r2 ≈ 1
        @test_throws ArgumentError fit_outcome(r,y;predictors=(:conventionality,:atypicality))
        @test_throws DimensionMismatch fit_outcome(r,y[1:3])
        @test_throws ArgumentError fit_outcome(r,fill(NaN,40))
        plots = [plot_feature_distributions(es), plot_pair_cloud(ps,1),
                 plot_pair_cloud(ps,(:a,:b)), plot_correlation_heatmap(es),
                 plot_pair_redundancy(ps), plot_typicality_profile(r,1),
                 plot_typicality_profile(r,:entity_id=>1), plot_atypicality_distribution(r),
                 plot_novelty_conventionality(r), plot_outcome_fit(out), plot_pair_heatmap(r)]
        for p in plots
            @test p isa TypicalityPlot
            recipe = RecipesBase.apply_recipe(Dict{Symbol,Any}(),p)
            @test length(recipe) == length(p.series)
            @test all(!isempty(x.args) for x in recipe)
        end
        @test_throws ArgumentError plot_pair_cloud(ps,(:a,:missing))
        @test_throws ArgumentError plot_typicality_profile(r,:entity_id=>999)
        @test_throws ArgumentError fit_typicality(ps;normalize=false)
        @test_throws ArgumentError fit_typicality(ps;k=-1)
        @test_throws ArgumentError fit_typicality(ps;k=40)
        @test_throws ArgumentError fit_typicality(ps;weights=ReferenceWeights(1:40))
        @test_throws ArgumentError fit_typicality(build_pairs(EntitySet([1. 2.],[:a,:b])))
        @test scoretable(r;raw=true).typicality == vec(permutedims(r.raw_scores))
        metadata = analysis_metadata(r)
        metadata.reference_weights[1] = 0
        @test r.configuration.reference_weights[1] > 0
        ps.entity_ids[1] = 100
        @test r.entity_ids[1] == 1
    end
end

include("data_integrity.jl")
include("neighbors.jl")
include("density.jl")
include("mathematical_properties.jl")
include("graph.jl")
include("workflow.jl")
include("numerical_limits.jl")
include("visualization.jl")
include("numerical_depth.jl")
include("numerical_density.jl")
include("numerical_neighbors.jl")
include("numerical_precision.jl")
include("numerical_contracts.jl")

@testset "Serial and threaded pair scoring agree" begin
    X = [sin(i*j)+i/20 for i in 1:40,j in 1:4]
    ps = build_pairs(EntitySet(X,[:a,:b,:c,:d]))
    for method in TYPICALITY_METHODS
        serial = fit_typicality(ps;method=method)
        parallel = fit_typicality(ps;method=method,parallel=true)
        @test serial.raw_scores == parallel.raw_scores
        @test serial.pair_scores == parallel.pair_scores
    end
end
