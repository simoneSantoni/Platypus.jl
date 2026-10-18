# Run with: julia --project benchmark/graph.jl
# Cumulative allocation is not peak resident RAM. Warm compilation first.
using Platypus
using Random

function measure_graph(n=1000, p=8, k=20)
    rng = MersenneTwister(20261018)
    entities = EntitySet(rand(rng, n, p), [Symbol("x", j) for j in 1:p])
    fit_graph_typicality(entities; k=k)
    GC.gc()
    trial = @timed fit_graph_typicality(entities; k=k)
    return (n=n, p=p, k=k, seconds=trial.time, allocated_bytes=trial.bytes,
            result_bytes=Base.summarysize(trial.value),
            edges=trial.value.configuration.n_edges,
            components=trial.value.configuration.n_components)
end

println("Julia ", VERSION, "; threads=", Threads.nthreads(), "; CPU=", Sys.CPU_NAME)
println(measure_graph())
