using SparseArrays

function graph_oracle(matrix, k)
    n = size(matrix, 1)
    distances = [sqrt(sum(abs2, matrix[i, :] - matrix[j, :])) for i in 1:n, j in 1:n]
    radii = [sort([distances[i, j] for j in 1:n if i != j])[k] for i in 1:n]
    adjacency = zeros(n, n)
    for i in 1:n, j in 1:n
        i == j && continue
        d = distances[i, j]
        if d <= radii[i] || d <= radii[j]
            adjacency[i, j] = d == 0 ? 1.0 :
                (radii[i] == 0 || radii[j] == 0 ? 0.0 : exp(-d^2 / (radii[i] * radii[j])))
        end
    end
    return adjacency
end

function graph_transition(result)
    n = result.n_entities
    transition = Matrix(result.adjacency)
    for i in 1:n
        if result.degree[i] > 0
            transition[i, :] ./= result.degree[i]
        else
            transition[i, i] = 1.0
        end
    end
    return transition
end

@testset "Full-feature exact graph" begin
    fixtures = [
        [0. 0.; 0. 1.; 1. 0.; 1. 1.],
        [0. 0.; 1. 0.; 10. 0.; 11. 0.],
        [0. 0.; 0. 0.; 0. 0.; 1. 0.],
        zeros(6, 3),
        [sin(1.7i + j) + i / 40 for i in 1:17, j in 1:4],
    ]
    for matrix in fixtures
        n, p = size(matrix)
        names = [Symbol("f", j) for j in 1:p]
        for k in unique([1, min(2, n - 1), n - 1])
            es = EntitySet(matrix, names)
            result = fit_graph_typicality(es; k=k, transform=:none)
            brute = fit_graph_typicality(es; k=k, transform=:none, backend=:brute)
            @test Matrix(result.adjacency) ≈ graph_oracle(matrix, k)
            @test result.adjacency == transpose(result.adjacency)
            @test result.adjacency == brute.adjacency
            @test result.typicality == brute.typicality
            @test sum(result.stationary) ≈ 1.0
            @test all(>=(0), result.stationary)
            @test transpose(graph_transition(result)) * result.stationary ≈ result.stationary
            @test result.atypicality == 1 .- result.typicality
            @test nrow(summarytable(result)) == n
            @test scoretable(result; raw=true).typicality == result.stationary
            permutation = reverse(collect(1:n))
            permuted = fit_graph_typicality(EntitySet(matrix[permutation, :], names);
                                            k=k, transform=:none)
            @test permuted.stationary[invperm(permutation)] == result.stationary
            @test permuted.typicality[invperm(permutation)] == result.typicality
            @test permuted.adjacency[invperm(permutation), invperm(permutation)] == result.adjacency
        end
    end

    constant = fit_graph_typicality(EntitySet(zeros(6, 3), [:x, :y, :z]); k=1)
    @test all(==(1 / 6), constant.stationary)
    @test constant.typicality ≈ fill(0.5, 6)
    @test constant.configuration.n_edges == 15 # all boundary ties, not six edges
    @test constant.configuration.stationary_unique
    @test_throws ArgumentError fit_graph_typicality(EntitySet(zeros(6, 3), [:x, :y, :z]); k=1, max_edges=14)

    one_feature = fit_graph_typicality(EntitySet(reshape([0., 1., 2.], 3, 1), [:x]);
                                       k=1, transform=:none)
    @test one_feature.stationary == [0.25, 0.5, 0.25]

    isolated = fit_graph_typicality(EntitySet([0. 0.; 0. 0.; 0. 0.; 1. 0.], [:x, :y]); k=2)
    @test isolated.stationary == [1 / 3, 1 / 3, 1 / 3, 0.0]
    @test isolated.degree == [2., 2., 2., 0.]
    @test isolated.configuration.n_components == 2
    @test !isolated.configuration.stationary_unique

    # Pairwise projections have identical multiplicities, but the complete
    # three-way configurations have a 10:1 multiplicity contrast.
    rows = [[Float64(x), Float64(y), Float64(z)] for x in 0:1 for y in 0:1 for z in 0:1
            for repeat in 1:(iseven(x + y + z) ? 10 : 1)]
    matrix = reduce(vcat, permutedims.(rows))
    full = fit_graph_typicality(EntitySet(matrix, [:x, :y, :z]); k=5)
    rare = [isodd(Int(sum(row))) for row in rows]
    @test all(iszero, full.stationary[rare])
    @test all(>(0), full.stationary[.!rare])
    for (a, b) in ((1, 2), (1, 3), (2, 3))
        @test Platypus._neighbor_lof(matrix[:, a], matrix[:, b], 5) == ones(size(matrix, 1))
    end

    matrix = [0. 1. 3.; 1. 4. 0.; 2. 2. 1.; 3. 0. 2.]
    es = EntitySet(matrix, [:x, :y, :z])
    reference = fit_graph_typicality(es; k=2)
    @test fit_graph_typicality(EntitySet(matrix .^ 3, [:x, :y, :z]); k=2).typicality == reference.typicality
    @test fit_graph_typicality(es; k=2, weights=ReferenceWeights(fill(1e-200, 4))).stationary == reference.stationary
    @test_throws ArgumentError fit_graph_typicality(es; weights=ReferenceWeights([1., 1., 1., 2.]))
    @test_throws ArgumentError fit_graph_typicality(es; k=-1)
    @test_throws ArgumentError fit_graph_typicality(es; k=4)
    @test_throws ArgumentError fit_graph_typicality(es; max_edges=0)
    @test_throws ArgumentError fit_graph_typicality(EntitySet(zeros(1, 2), [:x, :y]))
    ids = copy(reference.entity_ids)
    es.entity_ids[1] = 100
    @test reference.entity_ids == ids
    values = atypicality(reference)
    values[1] = -999
    @test reference.atypicality[1] >= 0
end
