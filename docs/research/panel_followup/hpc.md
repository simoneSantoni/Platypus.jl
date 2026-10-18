# Julia HPC follow-up: neighborhood and graph correctness

Review date: 2026-10-18. Julia 1.13.2, two threads. This is an investigation,
not a repair. Neither package source, existing tests, nor the original audit
outputs were changed. Human-perception validation is excluded.

## Conclusion and severity

The neighbor layer does not currently satisfy its advertised exact,
tie-inclusive Euclidean contract. There are three distinct numerical pathways:
rank-coordinate rounding changes intended rank-lattice ties; sequential `hypot`
changes distances under feature permutation; and row-dependent z-score moments
change geometry before search. Exact comparison of rounded distances can also
include genuinely unequal distances. These are high-priority correctness defects
because they alter discrete neighborhoods and graph edges, not merely the final
percentile of nearly equal raw scores.

KD-tree and brute agreement is insufficient: both call the same distance
function and final membership predicate. All 13 new oracle configurations agree
across backends, including the incorrect ones. No independent candidate-pruning
failure was observed; the existing radius guard has not thereby been proved
sufficient for every dimension and floating-point range.

## Independent oracle and interventions

The [probe](../../../validation/consistency/panel/hpc/probe.jl) computes squared
distances and kth-boundary membership with `Rational{BigInt}` independently of
package kernels. Affinities use 256-bit `BigFloat` evaluation of
`exp(-squared_distance/sqrt(squared_radius_i*squared_radius_j))`; degree and
stationary mass then follow directly. Only small fixtures use this quadratic
oracle. It is not proposed as a production implementation.

The [oracle results](../../../validation/consistency/panel/hpc/oracle_checks.csv)
and [intervention results](../../../validation/consistency/panel/hpc/interventions.csv)
separate topology, affinity, raw stationary mass, and normalized typicality.

### 1. Exact integer coordinates expose feature-order dependence

An eight-row, three-feature integer fixture, `transform=:none`, `k=2`, is saved
in [integer_feature_fixture.csv](../../../validation/consistency/panel/hpc/integer_feature_fixture.csv).
For entity 6, squared distances to entities 3 and 5 are exactly 26, while its
nearest neighbor 7 is at squared distance 6. Both boundary neighbors must be
included. Original feature order gives scaled distances
`2.5495097567963927` and `2.5495097567963922`; reversed order gives
`2.5495097567963922` for both. Original order therefore excludes entity 3.

The excluded directed relation is not restored by graph union: edge `(3,6)` is
missing. Reversing features restores the oracle topology, changes stationary
mass by **0.020025433368322146**, and changes typicality by **0.125**. The correct
edge affinity is approximately **0.214937**. No rank transform, non-dyadic input,
or covariance calculation is involved.

Root cause: [`src/neighbors.jl:46–54`](../../../src/neighbors.jl) applies `hypot`
successively in feature order for dimensions other than two. These operations
are not associative. The final `<=` predicate at lines 61–72 turns one-bit
distance differences into different graph structure. Sorting degree reductions
later cannot restore an edge that was never constructed.

### 2. Rank-lattice geometry and rounded coordinates are different contracts

The prior five-row LOF witness has two missing directed neighbors before
reflection and one afterward. This probe reproduces those differences with
exact integer squared-distance membership. Passing the integer ranks directly
with `transform=:none` restores all intended neighborhoods. The rank transform
`(tiedrank(x)-0.5)/n` introduces non-dyadic coordinate rounding; later exact
floating comparisons do not preserve its mathematical distance ties.

For this particular five-row fixture graph union restores every undirected
edge, so graph stationary mass remains within `5.56e-17` of the oracle. It is
incorrect to infer that every LOF neighborhood defect must change graph topology.
The eight-row fixture above demonstrates a graph defect separately.

### 3. Z-score preprocessing explains the observed grid row-order change

On saved `cartesian_grid_101`, z-score preprocessing followed by row permutation
changes one transformed feature by at most **4.440892098500626e-16** after
alignment. It changes **104 directed memberships**, **32 undirected edges**,
stationary mass by **0.0028561189599899638**, and typicality by **0.25**.

Freezing the original transformed coordinates, then fitting `transform=:none`
before and after row permutation, removes every membership, edge, mass and
typicality difference. The probe asserts the exact aligned membership and mass
identities. This localizes the permutation failure to `_zscore`'s ordinary
`mean`/`std` reductions at [`src/pairs.jl:71–79`](../../../src/pairs.jl), rather
than dictionary iteration or final degree reduction. Graph degree already uses
sorted `_stable_sum` calls at [`src/graph.jl:140–143`](../../../src/graph.jl).

There are two deliberately different grid oracles:

- The **ideal integer lattice** is the intended equally spaced grid/rank geometry.
  Against it, the original z-score fit omits 104 directed neighbors and 32 edges;
  the permuted fit, rank intervention, and raw integer intervention restore the
  ideal graph. Saved raw thirds are non-dyadic Float64 values, so those ideal-grid
  omission counts alone do not establish correctness relative to exact binary
  raw inputs. They are diagnostic comparisons, not interchangeable definitions.
- The **exact frozen Float64 geometry** converts every already-transformed value
  to its exact rational representation. Against this literal-input oracle, both
  row orders include **80 extra directed neighbors and 32 extra edges**; maximum
  stationary-mass error is **0.00115101**. Distinct exact squared distances have
  collapsed to equal floating distances. Freezing removes row-order dependence,
  but does not make the distance predicate mathematically exact.

The latter result and the integer fixture are independent, unambiguous failures
of exact-input Euclidean membership. An integer-grid z-score control additionally
misses four undirected edges against the ideal transformed geometry in both row
orders: deterministic preprocessing alone is necessary but insufficient.

## Numerical policy and algorithm screen

1. **Exact unweighted rank predicates are a strong fit.** Preserve doubled
   midranks as integers and compare sums of squared rank differences before
   taking square roots. Common rank translation and scale cancel from membership.
   Check accumulator overflow explicitly and choose wider arithmetic or fallback
   when the dimension/sample-size bound requires it. Weighted ranks need their
   own declared semantics; their cumulative masses are not this integer lattice.
2. **Filtered/adaptive exact comparisons suit arbitrary Float64 coordinates.**
   Evaluate a fast squared-distance comparison with a justified error bound;
   resolve ambiguous comparisons with an exact representation or sufficient
   higher precision. Preserve input information through preprocessing. This
   follows the robust-predicate design principle described by
   [Shewchuk](https://www.cs.cmu.edu/~quake/robust.html); his orientation/incircle
   routines are not a drop-in squared-distance implementation. An adaptive
   distance predicate and safe candidate bounds must be designed and tested.
3. **Deterministic scaled norm accumulation improves reproducibility** but is
   not an exact-membership proof. Sorting components may fix feature-order
   dependence at additional cost; compensated or reproducible summation can
   improve accuracy. Neither guarantees equality/order of boundary distances.
4. **Canonical preprocessing reductions** can remove row-order dependence.
   Stable centering is also needed for large offsets. A deterministic but biased
   moment estimate is not sufficient. Uniform and weighted paths should share
   a stated numerical contract before threading or SIMD optimization.
5. **Do not add a blanket epsilon.** It changes neighborhoods, can merge real
   near-neighbors, and does not define transitive score ties. A distinct policy
   for uncertain raw-score equality is still needed after topology is repaired;
   this review does not resolve the percentile layer.
6. **Retain exact search as the default.** Approximate retrieval would change
   the target graph and LOF estimator and would not repair these defects. The
   KD-tree can remain a candidate accelerator only if its pruning cannot omit
   candidates required by the chosen exact predicate. The current `8dimensions`
   successor-radius guard is a heuristic until its bound is justified.

The original [LOF definitions 3–4](https://cs.fit.edu/~pkc/ml/related/breunig00sigmod.pdf)
require all kth-distance ties. The
[Zelnik-Manor–Perona paper, equations 1–2](https://papers.nips.cc/paper/2004/file/40173ea48d9567f1f393b20c855bb40b-Paper.pdf)
supplies local-scale affinity; Platypus's sparse union graph and degree-percentile
score are additional package choices, not the paper's complete spectral
clustering algorithm. The official
[NearestNeighbors.jl documentation](https://github.com/KristofferC/NearestNeighbors.jl)
describes search APIs but does not establish exact real-arithmetic predicates.

## HPC review and acceptance criteria

No new broad benchmarks were run: correctness is the blocker. Neighbor kernels
avoid a dense distance matrix, but brute search remains quadratic, and sorted
neighbor reductions add sorting cost. Tie-inclusive graph construction may be
dense even for small `k`. `max_edges` bounds retained undirected edges, not work
already spent computing radii, candidate buffers, dictionary overhead or total
live memory. Pair threading multiplies temporary working storage. These are
existing engineering limits, not newly established bugs.

Require the following before claiming this layer corrected:

- Exact membership agreement on the five-row rank and eight-row integer fixtures,
  their reflections and permutations, with both backends independently checked.
- Exact-input and intended-transform oracles explicitly distinguished. Include
  near-but-unequal distances to prevent a tolerance rule from passing by merging
  genuinely distinct observations.
- Grid permutation invariance through preprocessing, graph topology and raw mass,
  plus the frozen-coordinate exact-membership fixture.
- High-precision affinity/mass agreement after topology agrees; separately verify
  normalized-score ties and stationary equations against the intended graph.
- Evidence or conservative fallback for candidate-pruning bounds, dimensional
  accumulation, overflow, underflow and rank accumulator limits.
- Then benchmark warmed kernels, allocations and peak memory across dimensions,
  boundary-tie multiplicities and thread counts. Julia's
  [performance guidance](https://docs.julialang.org/en/v1.12/manual/performance-tips/)
  distinguishes compilation from steady execution and warns that `@fastmath`
  permits changes to floating evaluation; it should not conceal a failed
  numerical contract.

The earlier HPC report's description of brute search as an oracle was too broad.
`test/neighbors.jl` independently organizes dense computation but uses the same
`hypot` and exact floating comparison: it is an implementation cross-check,
not an independent mathematical boundary oracle.

## Reproduction and procedure log

Run from the repository root:

```sh
julia --threads=2 --project=. --startup-file=no validation/consistency/panel/hpc/probe.jl
```

The procedure was: read the original simulation report, source and tests; consult
the primary LOF/local-scaling papers and Julia/NearestNeighbors documentation;
build a rational graph oracle; reproduce rank and grid cases; use seed 20261018
to find a small integer feature-permutation witness; freeze preprocessing to
isolate causality; record topology and raw/normalized changes separately.
A nested-function scope-capture mistake in the new diagnostic was corrected
before the final run. Earlier preliminary frozen-coordinate observations are
superseded by the saved final artifacts and assertions. Existing audit files
were never overwritten.
