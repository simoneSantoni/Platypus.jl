# Run with: julia --project benchmark/neighbors.jl
# Single pair, uniform cloud. Allocations are cumulative, not peak resident RAM.
using Platypus
using Random
using Printf

function measure_neighbors(n, k, backend)
    rng = MersenneTwister(20261018 + n)
    u, v = rand(rng, n), rand(rng, n)
    Platypus._neighbor_lof(u, v, k; backend=backend) # compilation warm-up
    GC.gc()
    run = @timed Platypus._neighbor_lof(u, v, k; backend=backend)
    return run.time, run.bytes
end

println("Julia ", VERSION, "; threads=", Threads.nthreads(), "; CPU=", Sys.CPU_NAME)
println("n,k,backend,seconds,allocated_bytes,dense_distance_matrix_bytes")
for n in (250, 1000, 10000), backend in (:kdtree, :brute)
    backend == :brute && n > 1000 && continue
    k = min(20, n - 1)
    seconds, bytes = measure_neighbors(n, k, backend)
    @printf("%d,%d,%s,%.6f,%d,%d\n", n, k, backend, seconds, bytes, 8n^2)
end
