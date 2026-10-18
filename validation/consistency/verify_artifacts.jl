# Verify that saved inputs match the generators and that result tables agree.
# This checks audit integrity, not whether Platypus passed the audit.
using CSV, DataFrames, SHA, TOML, InteractiveUtils, LinearAlgebra
include("datasets.jl")

root=isempty(ARGS) ? joinpath(@__DIR__,"results") : only(ARGS)
catalog=CSV.read(joinpath(root,"datasets.csv"),DataFrame)
checks=CSV.read(joinpath(root,"checks.csv"),DataFrame)
baselines=CSV.read(joinpath(root,"baselines.csv"),DataFrame)
counterexamples=CSV.read(joinpath(root,"counterexamples.csv"),DataFrame)
summary=CSV.read(joinpath(root,"summary.csv"),DataFrame)
scores=CSV.read(joinpath(root,"scores.csv"),DataFrame)
sizes=Dict(17=>32,101=>63,2026=>96)
for row in eachrow(catalog)
    generated=dataset(Symbol(row.scenario),row.seed;n=sizes[row.seed])
    saved=CSV.read(joinpath(root,"datasets",row.dataset*".csv"),DataFrame)
    @assert Matrix(saved[:,generated.names])==generated.matrix
    @assert saved.reference_weight==generated.weights
    @assert saved.planted==generated.truth && saved.label==generated.labels
    @assert saved.entity_id==collect(1:row.n)
    @assert size(generated.matrix)==(row.n,row.features)
end
@assert Set(checks.dataset) ⊆ Set(catalog.dataset)
@assert Set(baselines.dataset)==Set(catalog.dataset)
@assert nrow(unique(baselines,[:dataset,:method,:transform]))==nrow(baselines)
@assert nrow(unique(checks,[:dataset,:method,:transform,:check]))==nrow(checks)
flags=filter(:status=>==("violation"),checks)
@assert nrow(counterexamples)==nrow(flags)
@assert Set(Tuple.(eachrow(counterexamples[:,[:dataset,:method,:transform,:check]])))==
    Set(Tuple.(eachrow(flags[:,[:dataset,:method,:transform,:check]])))
@assert sum(summary.count)==nrow(checks)
for row in eachrow(summary)
    @assert row.count==count(x->x.method==row.method && x.transform==row.transform &&
        x.check==row.check && x.status==row.status,eachrow(checks))
end
for group in groupby(scores,[:dataset,:method,:transform])
    n=only(catalog.n[catalog.dataset .== first(group.dataset)])
    @assert sort(group.entity_id)==collect(1:n)
end
@assert length(groupby(scores,[:dataset,:method,:transform]))==count(==("ok"),baselines.status)
source=join(read.(sort(readdir(joinpath(@__DIR__,"..","..","src");join=true)),String),"\n")
fingerprint=bytes2hex(sha256(source))
@assert occursin("Source SHA256: $fingerprint",read(joinpath(root,"environment.txt"),String))

open(joinpath(root,"provenance.txt"),"w") do io
    println(io,"Artifact integrity verified against deterministic generators and result tables.")
    println(io,"Datasets: ",nrow(catalog),"; checks: ",nrow(checks),"; recorded violations: ",nrow(flags))
    println(io,"Package source SHA256 matches original run: ",fingerprint)
    versioninfo(io)
    println(io,"BLAS configuration: ",BLAS.get_config())
    println(io,"BLAS threads: ",BLAS.get_num_threads())
    println(io,"\nScript SHA256:")
    for name in ("datasets.jl","run.jl","reproduce.jl","plot_results.py","compare_repairs.py","verify_artifacts.jl")
        println(io,bytes2hex(sha256(read(joinpath(@__DIR__,name)))),"  ",name)
    end
    println(io,"\nResolved dependency versions:")
    projectdir=dirname(Base.active_project())
    versioned=joinpath(projectdir,"Manifest-v$(VERSION.major).$(VERSION.minor).toml")
    manifest=isfile(versioned) ? versioned : joinpath(projectdir,"Manifest.toml")
    println(io,"Manifest SHA256: ",bytes2hex(sha256(read(manifest))))
    dependencies=TOML.parsefile(manifest)["deps"]
    for name in sort(collect(keys(dependencies))), entry in dependencies[name]
        println(io,name," ",get(entry,"version","unversioned")," uuid=",entry["uuid"],
            " tree=",get(entry,"git-tree-sha1","stdlib/path"))
    end
end
println("Verified ",nrow(catalog)," datasets, ",nrow(checks)," checks and ",nrow(scores),
    " baseline scores. Package source fingerprint unchanged; ",nrow(flags)," violations remain.")
