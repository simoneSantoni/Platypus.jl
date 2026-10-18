using Test, Random

@testset "Exact binary rounding independent of caller precision" begin
    rng = MersenneTwister(541)
    values = [0.0, nextfloat(0.0), floatmin(Float64), prevfloat(floatmin(Float64)),
              1.0, prevfloat(1.0), nextfloat(1.0), floatmax(Float64)]
    append!(values, [reinterpret(Float64, rand(rng, UInt64) & 0x7fefffffffffffff) for _ in 1:100])
    for precision in (16, 53, 256)
        setprecision(BigFloat, precision) do
            for x in values
                q = Rational{BigInt}(x)
                @test Platypus._rational_float(q) === x
                @test Platypus._rational_float(-q) == -x
                @test Platypus._rational_sqrt(q*q) === x
                x == floatmax(Float64) && continue
                y = nextfloat(x)
                midpoint = (q + Rational{BigInt}(y))/2
                expected = iseven(reinterpret(UInt64,x)) ? x : y
                @test Platypus._rational_float(midpoint) === expected
                @test Platypus._rational_sqrt(midpoint^2) === expected
                epsilon = (Rational{BigInt}(y)-q)/8
                @test Platypus._rational_float(midpoint-epsilon) === x
                @test Platypus._rational_float(midpoint+epsilon) === y
            end
            overflow = Rational{BigInt}(floatmax(Float64)) + (big(1)<<970)
            @test Platypus._rational_float(overflow) == Inf
            @test Platypus._rational_float(overflow-1) == floatmax(Float64)
            @test Platypus._rational_sqrt(overflow^2) == Inf
            @test Platypus._rational_sqrt((overflow-1)^2) == floatmax(Float64)
        end
    end

    X = [sin(i*j)+i/8 for i in 1:12, j in 1:3]
    es = EntitySet(X, [:a,:b,:c])
    baseline = fit_graph_typicality(es; k=3)
    setprecision(BigFloat,16) do
        actual = fit_graph_typicality(es; k=3)
        @test actual.adjacency == baseline.adjacency
        @test actual.typicality == baseline.typicality
    end
end
