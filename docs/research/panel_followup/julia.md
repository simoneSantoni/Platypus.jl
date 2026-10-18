# Julia engineering panel: collection integrity and reproducibility

Reviewed 2026-10-18 on Julia 1.13.2, one thread. This is a focused review, not a
repair. No package source, tests, other panel reports, or original audit outputs
were changed. The current suite's 911 passing checks were reported by the
integration reviewer; this panel ran its own focused probes. Human-perception
validation is outside scope.

## Assessment

The numerical findings in the simulation, mathematics, statistics and HPC
reports remain unresolved. This review independently identifies an additional
structural defect: a public `PairSet` can contain a repeated pair and omit another
pair, pass validation, and silently change aggregate scores. It also confirms
three narrower provenance/input-contract issues. These are **not explanations
of the original 630 simulation flags**, whose harness did not mutate pair lists.

| Finding | Severity and scope |
|---|---|
| Duplicate/missing pair collections pass the all-pairs validator | High conditional consequence: silently changes scientific summaries after public collection mutation or manual construction |
| Exported metadata does not bind reference weights to entity IDs | Medium reproducibility defect when tables are reordered or filtered |
| Effective neighbor fallback is not recorded for pairwise fits | Medium provenance/performance ambiguity; no score defect demonstrated |
| Incoherent eight-argument `PairSet` provenance passes validation | Low-priority malformed-input diagnostic defect; ordinary `build_pairs` output is unaffected |

The existing snapshot repairs largely work. The controls verified separation of
caller matrices, entity snapshots, pair sources, result IDs, nested caller
provenance, exported DataFrames and metadata accessor copies.

## J1. Duplicate pairs silently become undeclared aggregation weights

Location: [`src/types.jl:251`](../../../src/types.jl), particularly lines 253 and
258–273; aggregation at [`src/estimation.jl:178`](../../../src/estimation.jl).
The documented contract is every unordered feature pair exactly once.
`_validate_pairset` checks stored length against `n_pairs`, each pair's indices,
names, dimensions and source coordinates. It checks neither uniqueness of
`(i,j)` nor completeness against the feature schema.

With three features the builder returns `(a,b)`, `(a,c)`, `(b,c)`. This single
operation preserves every check currently performed:

```julia
ps.pairs[2] = ps.pairs[1]
```

The validator accepts `(a,b)`, `(a,b)`, `(b,c)`, and the fitted result still
reports three pairs. Repeated `(a,b)` now occupies two of the three median
positions, effectively controlling conventionality; `(a,c)` contributes nothing.
This is a collection-selection defect, not floating-point ranking instability.

The deterministic 12-row integer fixture and every row's results are saved in
[`pair_fixture.csv`](../../../validation/consistency/panel/julia/pair_fixture.csv)
and [`pair_scores.csv`](../../../validation/consistency/panel/julia/pair_scores.csv).
All six pairwise configurations accepted the corrupted collection:

| Configuration | Entities with changed atypicality, out of 12 | Maximum absolute atypicality change |
|---|---:|---:|
| kNN | 7 | 0.416667 |
| LOF | 8 | 0.291667 |
| KDE | 7 | 0.583333 |
| Depth | 7 | 0.500000 |
| Gaussian copula | 7 | 0.583333 |
| Beta copula | 7 | 0.333333 |

For every configuration the probe asserts **exact equality**, without a
numerical tolerance, of both corrupted raw and normalized matrices to the
original matrices indexed by columns `[1,1,3]`. It independently computes the
changed median from those original columns. This isolates the cause from all
numerical defects discussed by the other panelists. Novelty and dispersion also
change; full measurements are in
[`pair_validation.csv`](../../../validation/consistency/panel/julia/pair_validation.csv).

A separate constructor probe supplies only `(a,b)` and `(b,c)` with a three-feature
schema and `n_pairs=2`; validation and all six fits accept that incomplete set.
By contrast, reversing the three distinct pair entries is an exact passing
control: aligned raw values and every entity's atypicality are unchanged.

Unmodified `build_pairs` does construct the correct collection. The defect is
conditional on subsequent public-array mutation or caller construction; it does
not establish that ordinary builder output corrupts itself. Direct mutations
are described as advanced operations, but an explicit validator must enforce
the collection contract it presents as valid.

**Acceptance criteria:** enforce one occurrence of every allowed unordered pair,
including `n_pairs == p*(p-1)/2`; reject both duplicate-with-missing and incomplete
collections at construction and/or fitting boundaries; retain valid arbitrary
pair ordering and existing identity/schema checks. Test the six methods against
these invalid objects and retain the reordered valid collection as a control.
If selected subsets or pair weighting are wanted, introduce an explicit distinct
contract rather than implicitly accepting them as an all-pairs `PairSet`.

## J2. Positional weights lose their identity binding in exported sidecars

Locations: [`src/estimation.jl:49`](../../../src/estimation.jl) and lines 188–194;
[`src/outcome.jl:94`](../../../src/outcome.jl). Documentation tells callers to save
`analysis_metadata(result)` alongside exported scores. The configuration includes
`reference_weights` but no associated `entity_ids`. IDs exist in the result and
summary table; the link is only their original row order.

The probe uses IDs `[30,10,40,20]` and normalized weights `[0.1,0.2,0.3,0.4]`.
After sorting the summary table by ID, the correct weights are
`[0.2,0.4,0.1,0.3]`. Pairing the exported positional sidecar with the sorted table
instead misassigns all four weights, by up to 0.2. Neither the sorted table nor
that metadata alone preserves the original ID-to-weight relation. This is a
reproducibility/export omission, **not an error in the fitted scores**, and the
original `TypicalityResult` still preserves the relation correctly.

[`metadata_binding.csv`](../../../validation/consistency/panel/julia/metadata_binding.csv)
records the true binding and the incompatible positional interpretation. Outcome
metadata inherits the same omission in its nested `measurement` configuration,
although `OutcomeResult.entity_ids` preserves the original order. Current graph
fits permit only uniform positive reference weights, so their metadata lacks
IDs too but does not currently have this individual-weight ambiguity.

**Acceptance criteria:** export immutable-by-copy identifier bindings with every
positional entity setting, for example matched `entity_ids`/`reference_weights`
vectors or keyed records. Keep the complete reference even when an exported
score table is filtered. Test sorted/subset tables, arbitrary unique IDs,
nonuniform weights, and nested outcome metadata. Apply the same principle to
pair-indexed settings if pair ordering is allowed to vary. Snapshot isolation
and identity binding are different properties; copying arrays establishes only
the first. This corroborates the sociology panel's provenance observation.

## J3. Requested and effective backend have different metadata semantics

Locations: [`src/neighbors.jl:35`](../../../src/neighbors.jl) through line 43;
[`src/estimation.jl:270`](../../../src/estimation.jl); compare
[`src/graph.jl:148`](../../../src/graph.jl).

For raw coordinates spanning `1e-100` to `1e100`, the neighbor setup deliberately
uses a linear scan even when `backend=:kdtree` was requested. The probe confirms
`search.tree === nothing`. The pairwise result records `backend=:kdtree`, while
the graph result records `backend=:brute` for the same dataset and request.
Explicit-brute and requested-KD-tree kNN raw/normalized scores match exactly; a
moderate-scale control actually constructs a tree.

Recording the requested option is legitimate. The defect is the **unlabeled
requested/effective ambiguity**, especially across two APIs whose metadata uses
the same key differently. It prevents the sidecar from identifying actual
algorithm execution and can obscure quadratic fallback cost. This is not a new
neighbor correctness counterexample or evidence that the fallback itself is bad.

**Acceptance criteria:** record requested backend separately from effective
backend, with the effective choice per pair when pairs can differ, and include
a fallback reason. Test an ordinary tree, extreme-range fallback, and explicit
brute request; require numerical agreement separately from accurate metadata.

## J4. Invalid transform/source combinations are accepted until incidental failure

Locations: [`src/types.jl:234`](../../../src/types.jl) and lines 265–274;
[`src/estimation.jl:173`](../../../src/estimation.jl) and line 304.

The default eight-field constructor accepts `transform=:rank` with
`source_matrix=nothing`. `_validate_pairset` accepts it; uniform `_fit_setup`
returns `columns=nothing`; `fit_copula` then raises
`MethodError: no method matching getindex(::Nothing, ::Int64)`.
This is a lower-priority malformed-input path, not a normal `build_pairs` failure.
The preserved six-argument constructor's coherent `:unknown`/`nothing` state
works for uniform kNN and correctly rejects nonuniform weighted fitting in the
probe.

**Acceptance criteria:** specify and validate permitted transform/source
combinations; reject missing copula provenance with a meaningful `ArgumentError`
before indexing. Keep the valid legacy constructor control. Julia inner
constructors are the appropriate place to establish cross-field invariants,
with later validation still required for mutable array contents. Merely typing
the fields does not express these relationships.
[Julia constructor documentation](https://docs.julialang.org/en/v1/manual/constructors/#Inner-Constructor-Methods)

## Architecture and testing implications

An immutable Julia `struct` can contain mutable vectors. Public mutable storage
is not itself a bug; the problem is an incomplete stated validation/snapshot
contract. A read-only wrapper, validated mutation API, or private owned storage
could reduce this surface, but none is required to fix J1. Do not destroy the
intended shared coordinate storage merely to prevent duplicate pair descriptors.
[Julia type documentation](https://docs.julialang.org/en/v1/manual/types/#Composite-Types)

`copy` is shallow and `deepcopy` recursively detaches nested values. The current
result/provenance copy controls pass; a wholesale claim that snapshots remain
broken would be unsupported.
[Julia copying documentation](https://docs.julialang.org/en/v1/base/base/#Base.copy)

No concurrent mutation was attempted. Sharing mutable inputs with fitting tasks
would require a separate synchronization contract; ordinary thread-count
agreement does not establish safety under caller mutation.
[Julia threading documentation](https://docs.julialang.org/en/v1/manual/multi-threading/#Communication-and-data-races-between-threads)

The existing data-integrity tests cover coordinate corruption, NaNs, schema/ID
validity and copy isolation, but not pair uniqueness/completeness or positional
metadata binding. The mathematical/HPC/statistical follow-ups also demonstrate
why toleranced raw equality and agreement between shared implementations can
miss discontinuities in neighborhoods and percentiles. These are distinct test
blindspots: add structural invariants, export round trips, and independent
mathematical controls rather than merely increasing the current pass count.

## Reproduction and evidence boundary

From the repository root:

```sh
JULIA_DEPOT_PATH=/tmp/platypus-review-julia-depot:/home/simon/.julia \
  julia --startup-file=no --project=. validation/consistency/panel/julia/probes.jl
```

Artifacts are contained in
[`validation/consistency/panel/julia/`](../../../validation/consistency/panel/julia/):
`probes.jl`, `run.log`, `summary.txt`, the fixture and result CSVs, and
`environment.txt` with Julia version and SHA-256 hashes of every source file and
Project.toml. The script asserts its causal identities and passing controls,
then records the current violations; its successful exit is not a claim those
violations are repaired. It overwrites only this panel's own diagnostic outputs.

Procedure: read the preserved simulation report and current mathematics/HPC/
statistics follow-ups; inspect constructors, validators, fitting and export
paths; revisit primary Julia documentation; construct pair-list interventions
with unchanged per-pair numerical scores; check copy controls, malformed
provenance and backend fallback; add the independently requested sorted-export
probe; rerun and preserve the observed outputs. No broad benchmark, numerical
algorithm replacement or perception validation was performed.
