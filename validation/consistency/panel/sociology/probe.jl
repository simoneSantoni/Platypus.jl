using Platypus, CSV, DataFrames

# Contract controls, not newly discovered numerical defects.
# Row splitting preserves weighted empirical mass but changes kernel-contribution
# leave-one-out exclusion. Fix bandwidth to isolate this distinction.
X = [0.0 0.0; 1.0 1.0; 2.0 0.0]
ids = [10, 20, 30]
w = [2.0, 1.0, 1.0]
a = fit_typicality(build_pairs(EntitySet(X, [:x, :y]; ids=ids); transform=:none);
    method=:kde, weights=ReferenceWeights(w), bandwidth=1.0)
Y = X[[1,1,2,3], :]
b = fit_typicality(build_pairs(EntitySet(Y, [:x, :y]; ids=[10,11,20,30]); transform=:none);
    method=:kde, weights=ReferenceWeights(ones(4)), bandwidth=1.0)
# For entity 10, original LOO density = [K(1,1)+K(2,0)]/2;
# split LOO density = [K(0,0)+K(1,1)+K(2,0)]/3.
kernel(x,y) = exp(-(x*x+y*y)/2)/(2pi)
expected_a = log((kernel(1,1)+kernel(2,0))/2)
expected_b = log((kernel(0,0)+kernel(1,1)+kernel(2,0))/3)
@assert isapprox(a.raw_scores[1,1], expected_a; atol=1e-14)
@assert isapprox(b.raw_scores[1,1], expected_b; atol=1e-14)
@assert isapprox(a.raw_scores[2,1], b.raw_scores[3,1]; atol=1e-14)
@assert isapprox(a.raw_scores[3,1], b.raw_scores[4,1]; atol=1e-14)
meta = analysis_metadata(a)
@assert :reference_weights in keys(meta)
# Record this artifact contract observation without treating the missing key as
# an intended behavior that future production tests should preserve.
CSV.write(joinpath(@__DIR__, "contract_controls.csv"), DataFrame(
    case=["original_weighted_entity", "split_duplicate_entity"],
    raw_log_density=[a.raw_scores[1,1],b.raw_scores[1,1]],
    independent_expected=[expected_a,expected_b],
    normalized_typicality=[a.pair_scores[1,1],b.pair_scores[1,1]],
    effective_sample_size=[meta.effective_sample_size,analysis_metadata(b).effective_sample_size]))
open(joinpath(@__DIR__, "metadata_observation.txt"), "w") do io
    println(io, "Julia ", VERSION)
    println(io, "result_entity_ids = ", a.entity_ids)
    println(io, "metadata_reference_weights = ", meta.reference_weights)
    println(io, "metadata_has_entity_ids = ", :entity_ids in keys(meta))
    println(io, "metadata_keys = ", keys(meta))
end
println("Contract controls passed; no numerical package defect is inferred from row splitting.")
