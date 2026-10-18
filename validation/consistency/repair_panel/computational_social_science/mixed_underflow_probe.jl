# Independent weighted Gaussian-mixture oracle for mixed represented-positive
# and all-underflow kernel masses. Records discrepancies; no package changes.
using Platypus, CSV, DataFrames, SHA
out=isempty(ARGS) ? joinpath(@__DIR__,"mixed_underflow") : only(ARGS)
mkpath(out)
X=[0. 0.;100. 0.;38. 0.;139. 0.]
rows=NamedTuple[]
for (case,w) in (("wide_weights",[1e-300,1.,0.,1e-200]),
                 ("moderate_weights",[1e-20,1.,0.,1e-10]))
    fit=fit_typicality(build_pairs(EntitySet(X,[:x,:y]);transform=:none);
        method=:kde,bandwidth=(1.,1.),weights=ReferenceWeights(w))
    oracle_logs,oracle_percentiles=setprecision(BigFloat,2048) do
        q=BigFloat.(w)
        logs=[log(sum(q[j]*exp(-(BigFloat(X[i,1])-BigFloat(X[j,1]))^2/2)
                      for j in 1:4 if j!=i)/sum(q[j] for j in 1:4 if j!=i))
              for i in 1:4]
        @assert logs[3] < logs[2] && logs[3] < logs[4]
        @assert logs[1] < logs[3]
        percentiles=[sum(q[j]*(logs[j]<logs[i] ? 1 : logs[j]==logs[i] ? .5 : 0)
                         for j in 1:4)/sum(q) for i in 1:4]
        (Float64.(logs),Float64.(percentiles))
    end
    for i in 1:4
        push!(rows,(case=case,entity_id=i,x=X[i,1],weight=w[i],
            oracle_log_kernel_mass=oracle_logs[i],raw_log_density=fit.raw_scores[i],
            normalization_score=fit.normalization_scores[i],
            oracle_percentile=oracle_percentiles[i],package_percentile=fit.pair_scores[i],
            percentile_error=fit.pair_scores[i]-oracle_percentiles[i]))
    end
end
CSV.write(joinpath(out,"observations.csv"),DataFrame(rows))
open(joinpath(out,"environment.txt"),"w") do io
    println(io,"Julia ",VERSION,"; Platypus ",pkgversion(Platypus),"; oracle precision 2048 bits")
    src=normpath(joinpath(@__DIR__,"../../../../src"))
    println(io,"Source SHA256: ",bytes2hex(sha256(join(read.(sort(readdir(src;join=true)),String),"\n"))))
    println(io,"Primary acceptance requirement: query3 is below reference2 by a well-separated log-density difference.")
    println(io,"This does not demand resolving every exponentially smaller distinction between reference2 and reference4.")
end
show(stdout,MIME("text/plain"),DataFrame(rows));println()
