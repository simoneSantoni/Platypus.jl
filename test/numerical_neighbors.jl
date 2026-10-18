# These oracles use exact input coordinates, not the production distance kernel.
function rational_neighbor_oracle(matrix, k)
    q = Rational{BigInt}.(matrix)
    n, p = size(matrix)
    squared = [sum((q[i,c] - q[j,c])^2 for c in 1:p) for i in 1:n, j in 1:n]
    boundary = [sort([squared[i,j] for j in 1:n if j != i])[k] for i in 1:n]
    membership = [i != j && squared[i,j] <= boundary[i] for i in 1:n, j in 1:n]
    adjacency = setprecision(BigFloat, 256) do
        a = zeros(BigFloat, n, n)
        for i in 1:n, j in i+1:n
            (membership[i,j] || membership[j,i]) || continue
            a[i,j] = a[j,i] = squared[i,j] == 0 ? BigFloat(1) :
                boundary[i] == 0 || boundary[j] == 0 ? BigFloat(0) :
                exp(-BigFloat(squared[i,j]) / sqrt(BigFloat(boundary[i]) * BigFloat(boundary[j])))
        end
        a
    end
    return membership, Float64.(adjacency)
end

function exact_search_membership(matrix, k, backend)
    search = Platypus._neighbor_setup_points(permutedims(matrix), backend)
    radii = Platypus._neighbor_kdistances(search, k)
    result = falses(size(matrix, 1), size(matrix, 1))
    for i in axes(matrix, 1)
        result[i, Platypus._neighbor_indices!(Int[], search, i, radii[i])] .= true
    end
    return result
end

function rational_lof_oracle(matrix, k)
    q = Rational{BigInt}.(matrix)
    n, p = size(matrix)
    squared = [sum((q[i,c] - q[j,c])^2 for c in 1:p) for i in 1:n, j in 1:n]
    boundary = [sort([squared[i,j] for j in 1:n if j != i])[k] for i in 1:n]
    neighbors = [[j for j in 1:n if i != j && squared[i,j] <= boundary[i]] for i in 1:n]
    return setprecision(BigFloat, 256) do
        distances = sqrt.(BigFloat.(squared))
        radii = sqrt.(BigFloat.(boundary))
        reach = [sum(max(radii[j], distances[i,j]) for j in neighbors[i]) /
                 length(neighbors[i]) for i in 1:n]
        Float64.([iszero(reach[i]) ? BigFloat(1) :
            sum(reach[i] / reach[j] for j in neighbors[i]) / length(neighbors[i]) for i in 1:n])
    end
end

@testset "Certified distance predicates and graph symmetry" begin
    integer = [0. 0. 3.; -2. 1. 1.; -1. -2. 2.; -1. 2. 1.;
               3. 3. 1.; 2. -1. -2.; 1. 0. 0.; -3. 1. 1.]
    ranks = [1. 2.; 2. 4.; 3. 5.; 4. 3.; 5. 1.]
    # Distinct exact distances with the SAME rounded norm must remain distinct.
    almost_tied = [0. 0.; 1. 2.0^-27; 1. 0.]
    @test hypot(1., 2.0^-27) == hypot(1., 0.)
    ref, _ = rational_neighbor_oracle(almost_tied, 1)
    @test findall(ref[1,:]) == [3]
    # Frozen non-dyadic coordinates are deliberately NOT an ideal integer grid.
    frozen = [Float64[x, y, z] for x in (-1., -1/3, 1/3, 1.)
              for y in (-1., -1/3, 1/3, 1.) for z in (-1., -1/3, 1/3, 1.)]
    frozen = reduce(vcat, permutedims.(frozen))
    for (matrix, k) in ((integer, 2), (ranks, 2), (almost_tied, 1), (frozen, 8),
                        ([0. 0.; 1e-300 0.; 1. 0.; 1e300 0.], 1))
        expected, affinity = rational_neighbor_oracle(matrix, k)
        p = size(matrix, 2)
        names = [Symbol("f", c) for c in 1:p]
        baseline = fit_graph_typicality(EntitySet(matrix, names); k, transform=:none)
        @test Matrix(baseline.adjacency) ≈ affinity atol=2e-15 rtol=2e-15
        @test baseline.configuration.requested_backend == :kdtree
        @test baseline.configuration.effective_backend == :brute
        @test baseline.configuration.distance_policy == :filtered_exact
        @test baseline.configuration.entity_ids == baseline.entity_ids
        for backend in (:kdtree, :brute)
            @test exact_search_membership(matrix, k, backend) == expected
            @test exact_search_membership(-matrix, k, backend) == expected
            @test exact_search_membership(matrix[:,end:-1:1], k, backend) == expected
            rows = reverse(collect(axes(matrix, 1)))
            @test exact_search_membership(matrix[rows,:], k, backend)[rows,rows] == expected
            for transformed in (-matrix, matrix[:,end:-1:1])
                actual = fit_graph_typicality(EntitySet(transformed, names); k, transform=:none, backend)
                @test actual.adjacency == baseline.adjacency
                @test actual.stationary == baseline.stationary
                @test actual.typicality == baseline.typicality
            end
        end
    end
    @test integer[6,:] - integer[3,:] == [3., 1., -4.]
    @test exact_search_membership(integer, 2, :kdtree)[6,3]

    # Rank geometry keeps integer midranks until distance evaluation.
    rank_result = fit_graph_typicality(EntitySet(ranks, [:x,:y]); k=2)
    _, rank_affinity = rational_neighbor_oracle(ranks, 2)
    @test Matrix(rank_result.adjacency) ≈ rank_affinity
    @test rank_result.configuration.geometry == :integer_midrank
    for matrix in (-ranks, ranks[:,[2,1]])
        actual = fit_graph_typicality(EntitySet(matrix, [:x,:y]); k=2)
        @test actual.adjacency == rank_result.adjacency
        @test actual.typicality == rank_result.typicality
    end
    kresult = fit_typicality(build_pairs(EntitySet(ranks, [:x,:y])); method=:knn, k=2)
    expected = Platypus._neighbor_kdistance(ranks[:,1], ranks[:,2], 2) ./ size(ranks,1)
    @test vec(kresult.raw_scores) ≈ -expected
    lof = fit_typicality(build_pairs(EntitySet(ranks, [:x,:y])); method=:lof, k=2)
    for matrix in (-ranks, ranks[:,[2,1]])
        actual = fit_typicality(build_pairs(EntitySet(matrix, [:x,:y])); method=:lof, k=2)
        @test actual.raw_scores == lof.raw_scores
        @test actual.pair_scores == lof.pair_scores
    end

    # Decimal rescaling may genuinely change a binary-coordinate tie. Check
    # each represented dataset against its OWN mathematical geometry and LOF.
    original = [0. 0.; 1. 3.; 2. 1.; 4. 2.; 8. 1.]
    for scale in (1e-200, 1e-10, 8., 1e10, 1e200)
        matrix = original .* scale
        membership, _ = rational_neighbor_oracle(matrix, 2)
        @test exact_search_membership(matrix, 2, :brute) == membership
        @test exact_search_membership(matrix, 2, :kdtree) == membership
        @test Platypus._neighbor_lof(matrix[:,1], matrix[:,2], 2) ≈ rational_lof_oracle(matrix, 2)
    end

    # Exact sqrt rounding at ordinary, subnormal, midpoint and overflow boundaries.
    rationals = [Rational{BigInt}(nextfloat(0.0))^2,
                 Rational{BigInt}(floatmax(Float64))^2,
                 BigInt(2)//BigInt(1), BigInt(1)//(BigInt(1)<<2150)]
    for q in rationals
        expected = setprecision(BigFloat, 8192) do
            Float64(sqrt(BigFloat(q)))
        end
        @test Platypus._neighbor_sqrt(q) === expected
    end
    midpoint = (Rational{BigInt}(1.) + Rational{BigInt}(nextfloat(1.))) / 2
    @test Platypus._neighbor_sqrt(midpoint^2) === 1.
    @test Platypus._neighbor_sqrt(midpoint^2 + (BigInt(1)//(BigInt(1)<<200))) === nextfloat(1.)
    overflow = Rational{BigInt}(floatmax(Float64)) + (BigInt(1)<<970)
    @test isinf(Platypus._neighbor_sqrt(overflow^2))
end
