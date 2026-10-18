# Focused downstream consequence audit; does not change package estimators or
# original consistency outputs. Run with the repository Julia project.
using Platypus, CSV, DataFrames, Statistics, SHA

const OUT = @__DIR__
const ROOT = normpath(joinpath(OUT,"../../../.."))

function summarize_effect(case, before, after; interpretation)
    length(before) == length(after) || error("unaligned observations")
    n = length(before)
    reversals = ties_split = ties_created = 0
    for i in 1:n-1, j in i+1:n
        a, b = sign(before[i]-before[j]), sign(after[i]-after[j])
        reversals += a*b < 0
        ties_split += a == 0 && b != 0
        ties_created += a != 0 && b == 0
    end
    # A deliberately illustrative operational cutoff, not a calibrated anomaly
    # decision. Count transitions in either direction rather than a false rate.
    cutoff = 0.75
    (case=case,n=n,interpretation=interpretation,
     changed_entities=count(abs.(before.-after) .> 1e-10),
     max_headline_delta=maximum(abs.(before.-after)),
     strict_order_reversals=reversals,computed_ties_split=ties_split,
     computed_ties_created=ties_created,
     selected_before=count(before .>= cutoff),selected_after=count(after .>= cutoff),
     cutoff_switches=count((before .>= cutoff) .!= (after .>= cutoff)))
end

function main()
    effects = NamedTuple[]
    scores = NamedTuple[]
    function record(case,before,after;interpretation)
        push!(effects,summarize_effect(case,before,after;interpretation))
        for i in eachindex(before)
            push!(scores,(case=case,entity_id=i,atypicality_before=before[i],
                atypicality_after=after[i],delta=after[i]-before[i]))
        end
    end
    M = [1. 2.;2. 4.;3. 5.;4. 3.;5. 1.]
    a = fit_typicality(build_pairs(EntitySet(M,[:x,:y]));method=:lof,k=2)
    b = fit_typicality(build_pairs(EntitySet(-M,[:x,:y]));method=:lof,k=2)
    record("lof_rank_reflection",a.atypicality,b.atypicality;
        interpretation="same distance geometry; numerical inconsistency")

    fixture = joinpath(ROOT,"validation/consistency/panel/hpc/integer_feature_fixture.csv")
    M = Matrix(CSV.read(fixture,DataFrame))
    a = fit_graph_typicality(EntitySet(M,[:x,:y,:z]);transform=:none,k=2)
    b = fit_graph_typicality(EntitySet(M[:,3:-1:1],[:z,:y,:x]);transform=:none,k=2)
    record("graph_feature_order",a.atypicality,b.atypicality;
        interpretation="same Euclidean feature geometry; numerical inconsistency")

    M = reduce(vcat,([Float64(x) Float64(y)] for x in 1:5 for y in 1:5))
    a = fit_copula(build_pairs(EntitySet(M,[:x,:y]));regularization=.05)
    augmented = vcat(M,M[1:2,:])
    b = fit_copula(build_pairs(EntitySet(augmented,[:x,:y]));regularization=.05,
        weights=ReferenceWeights(vcat(ones(25),0.,0.)))
    record("gaussian_zero_weight_append",a.atypicality,b.atypicality[1:25];
        interpretation="unchanged positive reference mass; numerical inconsistency")
    c = fit_copula(build_pairs(EntitySet(augmented,[:x,:y]));regularization=.05)
    record("gaussian_positive_weight_append",a.atypicality,c.atypicality[1:25];
        interpretation="changed empirical reference; descriptive sensitivity, not invariant")

    triangle = [0. 0.;1. 0.;0. 1.]
    r = fit_typicality(build_pairs(EntitySet(triangle,[:x,:y]));method=:depth)
    record("depth_exact_triangle_tie",fill(.5,3),r.atypicality;
        interpretation="exact equal-depth oracle versus computed output")
    r = fit_typicality(build_pairs(EntitySet(ones(3,2),[:x,:y]));method=:kde,
        weights=ReferenceWeights([1.,2.,3.]))
    record("kde_exact_constant_tie",fill(.5,3),r.atypicality;
        interpretation="exact constant-support oracle versus computed output")

    CSV.write(joinpath(OUT,"effects.csv"),DataFrame(effects))
    CSV.write(joinpath(OUT,"entity_scores.csv"),DataFrame(scores))
    open(joinpath(OUT,"environment.txt"),"w") do io
        println(io,"Julia ",VERSION,"; Platypus ",pkgversion(Platypus),"; threads=",Threads.nthreads())
        source=join(read.(sort(readdir(joinpath(ROOT,"src");join=true)),String),"\n")
        println(io,"Source SHA256: ",bytes2hex(sha256(source)))
        println(io,"Illustrative headline cutoff >=0.75; no calibrated decision or population-rate claim.")
        println(io,"Pair reversal counts use computed headline order; they do not certify exact mathematical ties.")
        println(io,"Triangle/constant-support expectations derive from the separately documented rational/algebraic oracles.")
    end
    show(stdout,MIME("text/plain"),DataFrame(effects))
    println()
end

main()
