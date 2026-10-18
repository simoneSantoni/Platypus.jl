@testset "Data integrity and reference preprocessing" begin
    @testset "Finite values and dimensions" begin
        for bad in (NaN, Inf, -Inf)
            @test_throws ArgumentError Entity(id=1, features=[bad])
            @test_throws ArgumentError EntitySet([1.0 2.0; bad 3.0], [:x, :y])
            @test_throws ArgumentError load_entities(DataFrame(x=[1.0,bad], y=[2.0,3.0]))
            @test_throws ArgumentError ReferenceWeights([1.0, bad])
        end
        @test_throws ArgumentError ReferenceWeights(Float64[])
        @test_throws ArgumentError ReferenceWeights([big"1e-1000", big"1.0"])
        @test_throws ArgumentError ReferenceWeights([big"-1e-1000", big"1.0"])
        @test_throws ArgumentError EntitySet(ones(2,2), [:x,:x])
        @test_throws ArgumentError EntitySet(ones(2,2), [:x,:y]; ids=[1,1])
        for length in (1,3)
            @test_throws DimensionMismatch EntitySet(ones(2,2), [:x,:y]; ids=collect(1:length))
            @test_throws DimensionMismatch EntitySet(ones(2,2), [:x,:y]; labels=fill("a",length))
        end
        es = EntitySet([1.0 2.0;3.0 4.0], [:x,:y])
        es.matrix[1,1] = NaN
        @test_throws ArgumentError build_pairs(es)
        @test_throws ArgumentError feature_summary(es)
        @test_throws ArgumentError nearest_entities(es,1)
        rw = ReferenceWeights([1.0, 1.0])
        rw.weights[1] = NaN
        @test_throws ArgumentError Platypus._validated_weights(rw,2)
        rw.weights .= 0
        @test_throws ArgumentError Platypus._validated_weights(rw,2)
        @test Platypus._validated_weights(ReferenceWeights([1e308,1e308]),2) == [0.5,0.5]
    end

    @testset "Canonical matrix and detached snapshots" begin
        vector = [1.0,2.0]
        e = Entity(1, vector, nothing)
        vector[1] = 99
        @test features(e) == [1,2]
        input = [e,Entity(id=2,features=[3.0,4.0])]
        names = [:x,:y]
        es = EntitySet(input,names)
        input[1].features[1] = 88
        names[1] = :edited
        @test es.matrix[1,1] == 1
        @test feature_names(es) == [:x,:y]
        es.matrix[1,1] = 42
        @test features(es[1]) == [42,2]
        @test features(first(es)) == [42,2]
        @test features(es.entities[1]) == [42,2]
        snapshot = es[1]
        snapshot.features[1] = 0
        features(es[1])[1] = -1
        feature_names(es)[1] = :changed
        @test es.matrix[1,1] == 42
        @test feature_names(es) == [:x,:y]
        ps = build_pairs(es)
        es.matrix[1,1] = 10
        es.entity_ids[1] = 99
        es.feature_names[1] = :changed
        @test ps.source_matrix[1,1] == 42
        @test ps.entity_ids == [1,2]
        @test ps.feature_names == [:x,:y]
    end

    @testset "Shared coordinates and weighted transforms" begin
        es = EntitySet([1.0 2.0 4.0;2.0 3.0 6.0;3.0 4.0 8.0], [:x,:y,:z])
        ps = build_pairs(es)
        @test ps[1].u === ps[2].u
        @test ps[1].v === ps[3].u
        @test ps.transform === :rank
        @test ps.source_matrix == es.matrix
        @test ps.source_matrix !== es.matrix
        for transform in (:rank,:zscore,:none)
            plain = Platypus._transform_columns(es.matrix,transform)
            uniform = Platypus._transform_columns(es.matrix,transform,ones(3))
            @test all(isapprox.(plain,uniform))
        end
        # Zero-weight rows do not change the reference transform, even at tails.
        M = reshape([-100.0,1.0,2.0,2.0,4.0,100.0],:,1)
        w = [0.0,1.0,2.0,1.0,1.0,0.0]
        ranks = only(Platypus._transform_columns(M,:rank,w))
        @test ranks ≈ [0,0.1,0.5,0.5,0.9,1]
        kept = only(Platypus._transform_columns(M[2:5,:],:rank,w[2:5]))
        @test ranks[2:5] ≈ kept
        for transform in (:rank,:zscore)
            @test only(Platypus._transform_columns(M,transform,w))[2:5] ≈
                only(Platypus._transform_columns(M[2:5,:],transform,w[2:5]))
            @test only(Platypus._transform_columns(M,transform,w)) ≈
                only(Platypus._transform_columns(M,transform,w .* 1e100))
        end
        # Mutated coordinates must not bypass the estimator's validation helper.
        ps[1].u[1] = Inf
        @test_throws ArgumentError Platypus._validate_pairset(ps)
        ps = build_pairs(es)
        ps[1].u[1] = 0.123
        @test_throws ArgumentError Platypus._validate_pairset(ps)
        ps = build_pairs(es; transform=:none)
        ps.source_matrix[1,1] = 999
        @test_throws ArgumentError Platypus._validate_pairset(ps)
    end

    @testset "EDA edge cases" begin
        es = EntitySet([1.0 2.0;1.0 3.0;1.0 4.0],[:constant,:varying])
        diagnostic = pair_redundancy(build_pairs(es))
        @test isnan(only(diagnostic.abs_correlation))
        @test !("weight" in names(diagnostic))
        @test_throws ArgumentError nearest_entities(es,1;k=-1)
        @test nrow(nearest_entities(es,1;k=0)) == 0
        @test build_pairs(EntitySet([1.0 2.0],[:x,:y]);transform=:zscore)[1].u == [0.0]
    end
end
