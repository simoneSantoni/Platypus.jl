# Small, independently checked fixtures, retained from the original failures.
# These record discrepancies; the package suite separately asserts repairs.
module ConsistencyReproducers
using Platypus, LinearAlgebra, CSV, DataFrames

function fit(X; method=:depth, transform=:none, k=0, weights=nothing)
    ps=build_pairs(EntitySet(X,[:x,:y]);transform)
    fit_typicality(ps;method,k,weights=weights===nothing ? nothing : ReferenceWeights(weights))
end

"""Exact population-covariance depth for a nonsingular rational input matrix."""
function exact_depth(X)
    R=Rational{BigInt}.(X); n=size(R,1)
    Z=R .- sum(R;dims=1)./n
    C=Z'*Z/n
    [1/(1+dot(Z[i,:],C\Z[i,:])) for i in axes(Z,1)]
end

"""Integer squared distances determine ties; BigFloat evaluates LOF afterward."""
function exact_rank_lof(X,k)
    @assert all(isinteger,X)
    n=size(X,1)
    @assert all(sort(X[:,c])==collect(1:n) for c in axes(X,2))
    ranks=Int.(X)
    sq=[sum(abs2,ranks[i,c]-ranks[j,c] for c in axes(ranks,2)) for i in 1:n,j in 1:n]
    radius=[sort([sq[i,j] for j in 1:n if j!=i])[k] for i in 1:n]
    neighbors=[[j for j in 1:n if j!=i && sq[i,j]<=radius[i]] for i in 1:n]
    raw=setprecision(BigFloat,256) do
        d=sqrt.(BigFloat.(sq)); r=sqrt.(BigFloat.(radius))
        reach=[sum(max(r[j],d[i,j]) for j in neighbors[i])/length(neighbors[i]) for i in 1:n]
        [inv(sum(reach[i]/reach[j] for j in neighbors[i])/length(neighbors[i])) for i in 1:n]
    end
    Float64.(raw),neighbors
end

function actual_neighbors(X,k)
    # Match the public rank-scoring path's exact doubled-midrank geometry.
    coordinates=Platypus._neighbor_rank_coordinates(X)
    search=Platypus._neighbor_setup_points(coordinates,:brute)
    radii=Platypus._neighbor_kdistances(search,k)
    [copy(Platypus._neighbor_indices!(Int[],search,i,radii[i])) for i in axes(X,1)]
end

function main(outdir=joinpath(@__DIR__,"results"))
    mkpath(outdir)
    rows=NamedTuple[]; inputs=NamedTuple[]; checks=NamedTuple[]
    function save(name,X,result,raw_expected;score_expected=fill(NaN,size(X,1)),weights=ones(size(X,1)))
        for i in axes(X,1)
            push!(inputs,(example=name,entity_id=i,x=X[i,1],y=X[i,2],weight=weights[i]))
            push!(rows,(example=name,entity_id=i,raw_expected=Float64(raw_expected[i]),
                raw_observed=result.raw_scores[i,1],score_expected=score_expected[i],
                score_observed=result.pair_scores[i,1],atypicality=result.atypicality[i]))
        end
    end
    function check(name,delta;tolerance=1e-10)
        push!(checks,(example=name,max_error=Float64(delta),tolerance=tolerance,
            status=delta<=tolerance ? "pass" : "violation"))
    end

    triangle=[0. 0.;1. 0.;0. 1.]
    expected=exact_depth(triangle)
    @assert all(==(1//3),expected)
    for transform in (:none,:rank)
        r=fit(triangle;transform)
        save("triangle_depth_$transform",triangle,r,expected;score_expected=fill(.5,3))
        check("triangle_depth_$(transform)_percentiles",maximum(abs.(r.pair_scores .- .5)))
    end

    constant=ones(3,2); weights=[1.,2.,3.]
    r=fit(constant;method=:kde,transform=:rank,weights)
    # Documented reduced-support KDE: density=1, hence log density=0.
    save("weighted_constant_kde",constant,r,zeros(3);score_expected=fill(.5,3),weights)
    check("weighted_constant_kde_percentiles",maximum(abs.(r.pair_scores .- .5)))

    ranks=[1. 2.;2. 4.;3. 5.;4. 3.;5. 1.]; k=2
    raw,neighbors=exact_rank_lof(ranks,k)
    reflected_raw,reflected_neighbors=exact_rank_lof(6 .- ranks,k)
    @assert raw==reflected_raw && neighbors==reflected_neighbors
    a=fit(ranks;method=:lof,transform=:rank,k)
    b=fit(-ranks;method=:lof,transform=:rank,k)
    save("lof_original",ranks,a,raw)
    save("lof_reflected",-ranks,b,raw)
    check("lof_original_raw_oracle",maximum(abs.(vec(a.raw_scores).-raw)))
    check("lof_reflected_raw_oracle",maximum(abs.(vec(b.raw_scores).-raw)))
    check("lof_reflection_percentiles",maximum(abs.(a.pair_scores.-b.pair_scores)))
    an,bn=actual_neighbors(ranks,k),actual_neighbors(-ranks,k)
    CSV.write(joinpath(outdir,"minimal_neighbors.csv"),DataFrame(
        entity_id=1:5,expected=join.(neighbors," "),original=join.(an," "),reflected=join.(bn," ")))

    offset=[1.421875 -0.1875;0.5625 2.078125;0.421875 -1.65625;1.625 0.3125;0.171875 -0.71875]
    shifted=offset .+ 1e14
    @assert shifted .- 1e14==offset # exact representability: no input information lost
    expected=exact_depth(offset)
    @assert exact_depth(shifted)==expected
    a,b=fit(offset),fit(shifted)
    save("depth_before_offset",offset,a,expected)
    save("depth_after_offset",shifted,b,expected)
    check("depth_before_offset_raw_oracle",maximum(abs.(vec(a.raw_scores).-Float64.(expected))))
    check("depth_after_offset_raw_oracle",maximum(abs.(vec(b.raw_scores).-Float64.(expected))))

    CSV.write(joinpath(outdir,"minimal_inputs.csv"),DataFrame(inputs))
    CSV.write(joinpath(outdir,"minimal_examples.csv"),DataFrame(rows))
    CSV.write(joinpath(outdir,"minimal_checks.csv"),DataFrame(checks))
    show(stdout,MIME"text/plain"(),DataFrame(checks));println()
    println("Exact oracle premises verified. Recorded violations: ",count(x->x.status=="violation",checks))
end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && ConsistencyReproducers.main(isempty(ARGS) ?
    joinpath(@__DIR__,"results") : only(ARGS))
