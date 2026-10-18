# Numerical repair: neighborhoods and full-feature graphs

Repair investigation: 2026-10-18, Julia 1.13.2, two threads. This note records
implementation decisions and validation for the HPC member of the panel.
Human-perception validity remains outside scope.

## Cause and repair

The old implementation used rounded Euclidean norms both as numerical outputs
and as predicates deciding which neighbors belong to a kth-distance boundary.
This conflated two different tasks. Sequential `hypot` depended on feature order;
rounded norm equality also included genuinely unequal squared distances. A
rounded radius cannot recover the exact boundary after information has been lost.

The repaired search retains each row's boundary vertex and compares squared
Euclidean distances. Every input Float64 is interpreted as its exact binary
rational value. A floating filter first encloses the exact squared distance by
rounding each subtraction, square and addition outward. Disjoint enclosures
certify ordering. Overlapping enclosures invoke exact `Rational{BigInt}`
arithmetic; this is an exact decision, without an epsilon or approximate tie
bucket. Subnormal values and overflow widen the filter rather than invalidate
its decision. Scratch storage is linear in the number of observations; there is
no retained dense distance matrix.

The returned norm is separately rounded to Float64. A scaled floating estimate
supplies a seed, then comparisons against **exact squared rounding midpoints**
select the correct Float64 square root, including ties-to-even, subnormals and
overflow. kNN restores the original units before that final rounding. No mutable
MPFR precision or rounding-mode state is used, so threaded pairs cannot interfere.

Unweighted rank geometry now reaches neighbor scoring as doubled integer
midranks. The common translation and scale cancel from LOF and graph affinities;
kNN converts distances back by dividing by `2n`. This preserves mathematical rank
lattice ties that were lost in rounded unit-square coordinates. The matrix has
feature-by-entity orientation. Raw and z-score geometry use exact predicates on
the transformed binary coordinates, not an assumed ideal integer grid. Parent
repairs handle the stability of those transformations themselves.

Graph affinity uses the exact canonical rational quantity
`distance_squared^2 / (radius_squared_i * radius_squared_j)`. Its square root,
followed by `exp(-value)`, evaluates the original self-tuning affinity. Equal
geometric expressions therefore produce identical floating affinities even when
separate rounded distance/radius calculations would disagree. Affinities,
degrees and stationary mass remain Float64 numerical outputs; this is not a claim
that transcendental values are represented exactly. Underflowed affinities retain
the existing zero-edge behavior. Graph metadata records entity IDs, rank/binary
geometry, requested backend, effective backend and exact predicate policy.

## Algorithm screen and tradeoffs

[Shewchuk's adaptive-predicate work](https://www.cs.cmu.edu/~quake/robust.html)
motivates separating a cheap certified filter from an exact fallback. His
orientation/incircle kernels are not reused as if they were squared-distance
predicates; this implementation derives enclosures for its own operations.
[Julia's arbitrary precision documentation](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic)
supports the integer arithmetic choice; the implementation avoids thread-unsafe
changes to global BigFloat precision.

The [original LOF definitions](https://cs.fit.edu/~pkc/ml/related/breunig00sigmod.pdf)
require every kth-distance tie. Approximate neighbor retrieval and tolerance
inclusion were rejected because they change this estimand. The
[Zelnik-Manor–Perona affinity](https://papers.nips.cc/paper/2004/file/40173ea48d9567f1f393b20c855bb40b-Paper.pdf)
is retained; sparse union topology and degree typicality remain package choices.

Neither NearestNeighbors.jl's rounded KD-tree pruning nor the earlier radius
inflation heuristic has a proof that it retains every candidate required by the
exact binary predicate. Consequently **both requested backends conservatively
use the filtered exact scan**. `backend` and `effective_backend` report `:brute`;
`requested_backend` retains the caller's choice, and
`distance_policy=:filtered_exact` states the algorithm. This costs quadratic
candidate work. A future certified KD-tree accelerator must match the independent
rational oracle before replacing the fallback. Adding a larger arbitrary guard
would not establish that guarantee.

## Independent regression evidence

`test/numerical_neighbors.jl` uses a separate dense rational oracle for small
fixtures and 256-bit BigFloat affinity evaluation. It covers:

- The eight-entity integer fixture where entity 6 must include both squared
  distance 26 boundary neighbors; graph edge `(3,6)` is restored.
- Five-entity rank-lattice membership and pair-score reflection/feature symmetry.
- A near-but-unequal fixture with squared distances `1` and `1 + 2^-54`, whose
  rounded `hypot` values are equal; only the truly nearer point is included.
- A frozen non-dyadic 64-point grid, checked against its literal Float64 geometry.
- Coordinates spanning `1e-300` to `1e300`, both requested backends, row/feature
  permutations and reflection, plus graph affinity and metadata assertions.
- Exact square-root rounding at normal, subnormal, overflow and midpoint cases.

The existing neighbor suite exposed an additional **test-oracle mistake**:
arbitrary decimal rescaling can round coordinates and genuinely break an exact
binary distance tie. The previous test treated every decimal scaling as an exact
Euclidean similarity. Exact powers of two are the appropriate invariance control;
rounded decimal inputs require their own rational geometry oracle. All five old
decimal-scale fixtures are retained in the new suite with independent rational
membership and 256-bit LOF checks for their represented coordinates.

The focused suites pass **150 existing neighbor checks, 219 existing graph
checks, and 162 new numerical checks** on Julia 1.13.2 with two threads.

Run the focused checks from the repository root:

```sh
julia --threads=2 --startup-file=no --project=. -e \
  'using Platypus, Test, DataFrames; include("test/neighbors.jl"); include("test/graph.jl"); include("test/numerical_neighbors.jl")'
```

## Warmed cost measurements

One warmed graph fit per configuration, seeded normal inputs (`MersenneTwister`
20261018), `transform=:none`, `k=11`, Julia 1.13.2, two threads. Construction and
compilation were warmed before timing, and garbage collection ran before each
measurement. Figures are diagnostic measurements on this workspace, not stable
performance guarantees. Allocation is cumulative allocated bytes, **not peak
resident memory**.

| Entities | Features | Seconds | Allocated MB | Undirected edges |
| ---: | ---: | ---: | ---: | ---: |
| 128 | 2 | 0.0084 | 9.85 | 870 |
| 512 | 2 | 0.0298 | 40.80 | 3,358 |
| 128 | 12 | 0.0159 | 27.02 | 1,027 |
| 512 | 12 | 0.1362 | 124.29 | 4,312 |

Exact ties can still yield dense graphs. `max_edges` limits retained edges;
quadratic candidate work and temporary exact arithmetic allocations remain.
The benchmark does not measure worst-case tie multiplicity, peak memory or a
speedup over the previous uncertified KD-tree implementation.

## Procedure record

1. Read the panel HPC review, saved integer fixture, rational probe, implementation
   and existing tests; distinguished ideal rank geometry from frozen Float64 data.
2. Revisited the primary robust-predicate, LOF and self-tuning references above.
3. Replaced rounded-norm membership with an outward filter plus exact fallback;
   introduced a separate exact-rounding algorithm for returned distance values.
4. Integrated parent-provided integer rank coordinates, canonicalized graph
   affinity arithmetic, and made conservative backend fallback explicit.
5. Added independent mathematical regressions, ran old and new focused tests,
   and identified the invalid decimal-scaling invariant in the old test oracle.
6. Timed warmed corrected fits after geometry checks. Repository-wide tests and
   the broad before/after consistency audit are recorded in the panel repair report.
