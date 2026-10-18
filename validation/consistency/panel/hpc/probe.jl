module HPCProbe
using Platypus, CSV, DataFrames, Random, StatsBase, LinearAlgebra

const ROOT=normpath(joinpath(@__DIR__,"../.."))

"Exact squared distances and membership; high-precision affinity evaluation."
function oracle(X,k)
    R=Rational{BigInt}.(X); n,p=size(R)
    sq=[sum((R[i,c]-R[j,c])^2 for c in 1:p) for i in 1:n,j in 1:n]
    radius=[sort([sq[i,j] for j in 1:n if j!=i])[k] for i in 1:n]
    membership=[i!=j && sq[i,j]<=radius[i] for i in 1:n,j in 1:n]
    A=setprecision(BigFloat,256) do
        a=zeros(BigFloat,n,n)
        for i in 1:n,j in i+1:n
            (membership[i,j] || membership[j,i]) || continue
            v=sq[i,j]==0 ? BigFloat(1) : radius[i]==0 || radius[j]==0 ? BigFloat(0) :
                exp(-BigFloat(sq[i,j])/sqrt(BigFloat(radius[i])*BigFloat(radius[j])))
            a[i,j]=a[j,i]=v
        end
        a
    end
    degree=vec(sum(A;dims=2)); total=sum(degree)
    mass=total==0 ? fill(BigFloat(1)/n,n) : degree/total
    (membership=membership, adjacency=Float64.(A), mass=Float64.(mass))
end

function inspect(X,k;transform=:rank,backend=:brute)
    columns=Platypus._transform_columns(X,transform)
    points=permutedims(hcat(columns...))
    s=Platypus._neighbor_setup_points(points,backend)
    radii=Platypus._neighbor_kdistances(s,k)
    membership=falses(size(X,1),size(X,1))
    for i in axes(X,1)
        membership[i,Platypus._neighbor_indices!(Int[],s,i,radii[i])].=true
    end
    graph=fit_graph_typicality(EntitySet(X,[Symbol("f$c") for c in axes(X,2)]);k,transform,backend)
    (membership=membership,graph=graph,points=permutedims(s.points),search=s)
end

function main()
    rows=NamedTuple[]
    interventions=NamedTuple[]
    function compare(label,a,b;order=collect(axes(a.membership,1)))
        aligned=invperm(order)
        push!(interventions,(intervention=label,
            membership_changes=count(a.membership.!=b.membership[aligned,aligned]),
            edge_changes=count((a.graph.adjacency.>0).!=(b.graph.adjacency[aligned,aligned].>0))÷2,
            max_mass_delta=maximum(abs.(a.graph.stationary.-b.graph.stationary[aligned])),
            max_typicality_delta=maximum(abs.(a.graph.typicality.-b.graph.typicality[aligned]))))
    end
    function record(label,X,k,reference;transform=:rank,order=collect(axes(X,1)))
        local a=inspect(X,k;transform,backend=:brute)
        local b=inspect(X,k;transform,backend=:kdtree)
        invorder=invperm(order)
        edge(a)=Matrix(a.graph.adjacency)[invorder,invorder]
        push!(rows,(fixture=label,transform=string(transform),n=size(X,1),p=size(X,2),k=k,
            directed_membership_errors=count(a.membership[invorder,invorder].!=reference.membership),
            directed_missing=count(reference.membership .& .!a.membership[invorder,invorder]),
            directed_extra=count(.!reference.membership .& a.membership[invorder,invorder]),
            undirected_edge_errors=count((edge(a).>0).!=(reference.adjacency.>0))÷2,
            max_affinity_error=maximum(abs.(edge(a).-reference.adjacency)),
            max_mass_error=maximum(abs.(a.graph.stationary[invorder].-reference.mass)),
            backend_membership_errors=count(a.membership.!=b.membership),
            backend_mass_error=maximum(abs.(a.graph.stationary.-b.graph.stationary))))
        a
    end

    X=[1. 2.;2. 4.;3. 5.;4. 3.;5. 1.]; ref=oracle(X,2)
    record("five_ranks_original",X,2,ref)
    record("five_ranks_reflected",-X,2,ref)
    record("five_ranks_integer_intervention",X,2,ref;transform=:none)

    # Dyadic integer lattice removes rank coordinate rounding: any remaining
    # feature-order error comes from distance accumulation, not preprocessing.
    rng=MersenneTwister(20261018)
    for trial in 1:10000
        X=Float64.(rand(rng,-3:3,8,3)); k=2
        a=inspect(X,k;transform=:none); b=inspect(X[:,3:-1:1],k;transform=:none)
        if (a.graph.adjacency.>0)!=(b.graph.adjacency.>0)
            ref=oracle(X,k)
            record("integer_feature_original",X,k,ref;transform=:none)
            record("integer_feature_reversed",X[:,3:-1:1],k,ref;transform=:none)
            CSV.write(joinpath(@__DIR__,"integer_feature_fixture.csv"),DataFrame(X,[:x,:y,:z]))
            CSV.write(joinpath(@__DIR__,"integer_feature_membership.csv"),DataFrame(
                entity=1:8,oracle=[join(findall(ref.membership[i,:])," ") for i in 1:8],
                original=[join(findall(a.membership[i,:])," ") for i in 1:8],
                reversed=[join(findall(b.membership[i,:])," ") for i in 1:8]))
            compare("integer_feature_reversal",a,b)
            CSV.write(joinpath(@__DIR__,"integer_feature_distances.csv"),DataFrame(
                entity=fill(6,7),neighbor=[j for j in 1:8 if j!=6],
                exact_squared=[sum(abs2,X[6,:].-X[j,:]) for j in 1:8 if j!=6],
                original_scaled=[Platypus._neighbor_distance(a.search,6,j) for j in 1:8 if j!=6],
                reversed_scaled=[Platypus._neighbor_distance(b.search,6,j) for j in 1:8 if j!=6]))
            break
        end
        trial==10000 && error("No feature permutation witness found")
    end

    # 4^3 Cartesian points. Exact integer grid is the mathematical rank geometry;
    # the saved non-dyadic raw values are used only by the zscore experiments.
    table=CSV.read(joinpath(ROOT,"results/datasets/cartesian_grid_101.csv"),DataFrame)
    X=Matrix{Float64}(table[:,r"^f"]); k=8
    lattice=hcat([tiedrank(X[:,c]) for c in axes(X,2)]...)
    ref=oracle(lattice,k)
    permutation=randperm(MersenneTwister(132),size(X,1))
    a=record("grid_original",X,k,ref;transform=:zscore)
    b=record("grid_permuted",X[permutation,:],k,ref;transform=:zscore,order=permutation)
    compare("grid_row_permutation",a,b;order=permutation)
    record("grid_rank_intervention",X,k,ref;transform=:rank)
    record("grid_integer_intervention",lattice,k,ref;transform=:none)
    record("grid_integer_zscore",lattice,k,ref;transform=:zscore)
    record("grid_integer_zscore_permuted",lattice[permutation,:],k,ref;transform=:zscore,order=permutation)
    # Freeze already transformed geometry: isolates row-order effects after preprocessing.
    freeze=oracle(a.points,k)
    frozen=record("grid_frozen_geometry",a.points,k,freeze;transform=:none)
    frozen_permuted=record("grid_frozen_geometry_permuted",a.points[permutation,:],k,freeze;transform=:none,order=permutation)
    compare("grid_frozen_row_permutation",frozen,frozen_permuted;order=permutation)
    @assert frozen.membership==frozen_permuted.membership[invperm(permutation),invperm(permutation)]
    @assert frozen.graph.stationary==frozen_permuted.graph.stationary[invperm(permutation)]
    coords=b.points[invperm(permutation),:]
    CSV.write(joinpath(@__DIR__,"grid_transform_delta.csv"),DataFrame(
        feature=1:3,max_coordinate_delta=[maximum(abs.(a.points[:,c].-coords[:,c])) for c in 1:3]))

    CSV.write(joinpath(@__DIR__,"oracle_checks.csv"),DataFrame(rows))
    CSV.write(joinpath(@__DIR__,"interventions.csv"),DataFrame(interventions))
    show(stdout,MIME"text/plain"(),DataFrame(rows));println()
    open(joinpath(@__DIR__,"environment.txt"),"w") do io
        println(io,"Julia ",VERSION,"; threads=",Threads.nthreads(),"; seed=20261018")
        println(io,"Exact Rational{BigInt} squared distances; 256-bit BigFloat affinities.")
        println(io,"No source estimator edits; oracle is quadratic and intended only for small fixtures.")
    end
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && HPCProbe.main()
