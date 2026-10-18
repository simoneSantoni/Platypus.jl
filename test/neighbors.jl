# Independent small dense oracle: clarity takes precedence over performance here.
function neighbor_oracle(u, v, k)
    n = length(u)
    d = [hypot(u[i] - u[j], v[i] - v[j]) for i in 1:n, j in 1:n]
    radii = [sort([d[i, j] for j in 1:n if i != j])[k] for i in 1:n]
    neighbors = [[j for j in 1:n if i != j && d[i, j] <= radii[i]] for i in 1:n]
    reach = [sum(max(radii[j], d[i, j]) for j in neighbors[i]) / length(neighbors[i]) for i in 1:n]
    lof = [iszero(reach[i]) ? 1.0 :
           sum(reach[i] / reach[j] for j in neighbors[i]) / length(neighbors[i]) for i in 1:n]
    return radii, lof
end

@testset "Exact neighbor kernels" begin
    fixtures = [
        ([0., 0., 0., 1., 2.], [0., 1., 2., 0., 1.]),
        (Float64[x for x in -2:2 for y in -2:2], Float64[y for x in -2:2 for y in -2:2]),
        ([sin(1.7i) + i / 40 for i in 1:37], [cos(2.3i) for i in 1:37]),
        ([0., 0., 0., 1., 2., 2.], [0., 0., 0., 0., 1., 1.]),
        (zeros(8), zeros(8)),
        ([0., 1.], [0., 0.]),
    ]
    for (u, v) in fixtures, k in unique([1, min(2, length(u) - 1), length(u) - 1])
        radii, lof = neighbor_oracle(u, v, k)
        actual_radii = Platypus._neighbor_kdistance(u, v, k)
        actual_lof = Platypus._neighbor_lof(u, v, k)
        @test actual_radii ≈ radii
        @test actual_lof ≈ lof
        @test actual_radii == Platypus._neighbor_kdistance(u, v, k; backend=:brute)
        @test actual_lof == Platypus._neighbor_lof(u, v, k; backend=:brute)
        for perm in (reverse(collect(eachindex(u))), sortperm(sin.(eachindex(u))))
            @test Platypus._neighbor_kdistance(u[perm], v[perm], k)[invperm(perm)] == actual_radii
            @test Platypus._neighbor_lof(u[perm], v[perm], k)[invperm(perm)] == actual_lof
        end
    end

    # Exact binary changes of units preserve the represented geometry. Decimal
    # scaling may split a boundary tie through input rounding; those cases are
    # checked against independent exact-input oracles in numerical_neighbors.jl.
    u, v = [0., 1., 2., 4., 8.], [0., 3., 1., 2., 1.]
    radii, lof = neighbor_oracle(u, v, 2)
    for scale in (2.0^-660, 2.0^-30, 8., 2.0^30, 2.0^660)
        @test Platypus._neighbor_kdistance(u .* scale, v .* scale, 2) ./ scale ≈ radii
        @test Platypus._neighbor_lof(u .* scale, v .* scale, 2) ≈ lof
    end

    # Extremely heterogeneous magnitudes take the robust linear-scan path.
    u = [0., 1e-300, 1., 1e300]
    @test Platypus._neighbor_kdistance(u, zeros(4), 1) == [1e-300, 1e-300, 1., 1e300]
    @test !any(isnan, Platypus._neighbor_lof(u, zeros(4), 1))

    # Explicit singular extension, with no added length scale.
    @test Platypus._neighbor_lof([0., 0., 0., 1.], zeros(4), 2) == [1., 1., 1., Inf]
    @test Platypus._neighbor_kdistance([0., 0., 0., 1.], zeros(4), 2) == [0., 0., 0., 1.]
    @test_throws ArgumentError Platypus._neighbor_lof([0.], [0.], 1)
    @test_throws ArgumentError Platypus._neighbor_lof([0., NaN], zeros(2), 1)
    @test_throws ArgumentError Platypus._neighbor_kdistance([0., Inf], zeros(2), 1)
    @test_throws ArgumentError Platypus._neighbor_lof(zeros(3), zeros(3), 0)
    @test_throws ArgumentError Platypus._neighbor_lof(zeros(3), zeros(3), 3)
    @test_throws ArgumentError Platypus._neighbor_lof(zeros(3), zeros(3), 1; backend=:unknown)
    @test_throws DimensionMismatch Platypus._neighbor_lof(zeros(3), zeros(2), 1)
    @test_throws ArgumentError Platypus._neighbor_kdistance([-1e308, 1e308], zeros(2), 1)
end
