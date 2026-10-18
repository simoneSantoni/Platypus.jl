# Exact neighborhood computation: HPC and Julia review

Date: 2026-10-18. Scope: numerical correctness, deterministic definitions,
memory use, and computational cost. Human perceptions are outside this review.

## Procedure and evidence

1. Inspected `src/estimation.jl`, `Project.toml`, and the existing test suite.
   The previous implementation allocated a dense distance matrix and sorted a
   complete row for every observation. LOF selected exactly k neighbors, even
   when additional observations had the same boundary distance; fixed `1e-12`
   additions also introduced a length scale into otherwise scale-free ratios.
2. Searched the primary LOF publication, official NearestNeighbors.jl repository
   and versioned source, and primary descriptions of LoOP, INFLO, flexible KDE
   outlier detection, and HNSW. Search phrases included “NearestNeighbors.jl
   KDTree knn inrange documentation”, “Breunig Kriegel Ng Sander LOF ties”,
   “Kriegel Schubert Zimek LoOP”, and “Jin Tung Han Wang symmetric neighborhood”.
3. Compared candidate search structures and estimator definitions separately.
   Replacing a search structure should preserve scores. Replacing LOF with
   another estimator changes the scientific quantity and needs its own API.
4. Implemented the chosen search backend in `src/neighbors.jl`, an independent
   dense oracle and invariance fixtures in `test/neighbors.jl`, and a reproducible
   benchmark in `benchmark/neighbors.jl`. The oracle is deliberately small and
   straightforward; production code never constructs its distance matrix.
5. Ran the standalone kernels before integration. Validation results and timings
   appear below; they describe this execution environment rather than a general
   hardware guarantee.

## Algorithm screen and decisions

| Candidate | Relevant properties | Decision for this phase |
| --- | --- | --- |
| Dense all-pairs distances | Simple reference implementation; quadratic memory | Keep only in small independent tests |
| Streaming brute-force search | Exact, linear working memory; quadratic distance work | Retain `backend=:brute` internally for verification and extreme numeric range |
| Exact KD-tree search | Suitable for the package's two-dimensional Euclidean projections | Adopt with a maintained Julia library |
| Ball tree | Supports broader metrics; unnecessary generality for these projections | Do not add a second production backend without workload evidence |
| HNSW | Approximate graph search, useful for high-dimensional retrieval | Defer: missing neighbors or ties change the target score |
| Tie-inclusive original LOF | Relative local reachability, already a promised estimator | Correct and retain |
| LoOP | Changes local-distance statistics and maps to a probability-like scale | Defer as a distinct estimator; do not claim rank scores are calibrated probabilities |
| INFLO | Adds reverse-neighbor information to local comparisons | Defer as a distinct estimator; it is not a drop-in LOF optimization |
| Flexible KDE outlier factors | Separates density estimation from density comparison | Retain this design distinction; density improvements are documented separately |

NearestNeighbors.jl documents KDTree, BallTree, and brute search, with reusable
`knn!` and `inrange!` buffers. The implementation uses these longstanding calls,
not recently added all-points convenience APIs. The installed version tested was
0.4.29 on Julia 1.10.12; the declared 0.4 compatibility permits Julia's resolver
to choose a compatible earlier release for older supported Julia versions.
[Official repository](https://github.com/KristofferC/NearestNeighbors.jl),
[versioned kNN API](https://raw.githubusercontent.com/KristofferC/NearestNeighbors.jl/v0.4.22/src/knn.jl).

The original LOF definitions include every observation at the kth distance.
Reachability uses the neighbor's kth distance, not the focal observation's.
The paper explicitly excludes duplicate configurations with infinite local
reachability density. Our duplicate extension is therefore a package decision,
not a claim that the publication defines this case: a zero-reachability cluster
receives LOF 1; an observation with positive reachability and a zero-reachability
neighbor receives LOF infinity, hence inverse-LOF typicality zero.
[Breunig, Kriegel, Ng and Sander, 2000, definitions 3–7](https://cs.fit.edu/~pkc/ml/related/breunig00sigmod.pdf).

LoOP and INFLO provide substantive alternatives developed by established
researchers. LoOP adds a probabilistic interpretation under its modeling
construction; INFLO uses both direct and reverse neighbors. Neither establishes
that it uniformly improves aggregated pairwise atypicality. Adopting them merely
because they are more elaborate would change the estimand without resolving the
identified correctness defects.
[Kriegel, Kröger, Schubert and Zimek, 2009](https://doi.org/10.1145/1645953.1646195),
[Jin, Tung, Han and Wang, 2006](https://hanj.cs.illinois.edu/pdf/pakdd06_outlier.pdf).

Schubert, Zimek and Kriegel's flexible KDE framework emphasizes separating
density estimation from the outlier comparison rule. That distinction informs
the implementation's separate kernels and documented score scales.
[Schubert, Zimek and Kriegel, 2014](https://epubs.siam.org/doi/10.1137/1.9781611973440.63).

HNSW is explicitly approximate. Its retrieval speed does not imply exact
neighborhood recovery, especially of all ties. That tradeoff is unattractive for
two-dimensional projections during this correctness-focused phase.
[Malkov and Yashunin](https://arxiv.org/abs/1603.09320).

## Numerical and memory contract

- Kernels reject invalid k, nonfinite coordinates, incompatible lengths, and
  unrepresentable distances. They accept vector views.
- A shared power-of-two coordinate rescaling protects squared-distance tree
  calculations without introducing an arbitrary epsilon into the estimator.
  Extremely heterogeneous coordinate exponents trigger robust `hypot` scans.
- A few floating-point successor steps enlarge only the candidate search radius;
  final tie membership uses `hypot(distance) <= kth_distance` exactly. Candidate
  distances are recomputed before selecting the kth radius.
- Sorted scalar reductions make results independent of neighbor traversal order.
  Zero-radius duplicate neighborhoods use their analytical value directly.
- kNN exposes the kth distance, including zero. Integration can use its negative
  as a finite typicality-oriented score; reciprocating zero is unnecessary.
- LOF computes ratios of mean reachability distances directly. This avoids
  creating infinite local densities before taking their ratios.
- Production working storage is linear in n plus a reusable neighborhood buffer.
  A buffer can contain n-1 observations when ties require it; neighborhoods are
  recomputed rather than all retained. Adverse geometry or extensive ties can
  require quadratic distance work, with additional sorting for deterministic
  reductions (up to O(N² log N) work). No universal subquadratic claim is made.
- Kernels use local working buffers and no shared mutable state. The integrated
  pairwise API can distribute independent feature pairs over Julia threads with
  fixed result columns. Each concurrent pair needs its own working storage;
  kernel-level threading and distributed execution are deferred.

## Validation and benchmark results

Standalone kernel tests: **150/150 passed**, Julia 1.10.12 and NearestNeighbors
0.4.29. Fixtures cover an independent dense LOF oracle, exact agreement of KD
and brute backends, row permutations, ties, duplicates, k=1 and k=n-1, scales
from 1e-200 to 1e200, extreme dynamic range, and invalid inputs.

Warm-up followed by one timed measurement per case, one Julia thread, CPU
reported as `goldmont`, uniform seeded cloud, k=20:

| n | Backend | Seconds | Cumulative allocation | Dense distance matrix alone |
| ---: | --- | ---: | ---: | ---: |
| 250 | KD tree | 0.008224 | 715,856 bytes | 500,000 bytes |
| 250 | Brute | 0.001549 | 543,664 bytes | 500,000 bytes |
| 1,000 | KD tree | 0.036459 | 2,867,808 bytes | 8,000,000 bytes |
| 1,000 | Brute | 0.017378 | 8,192,336 bytes | 8,000,000 bytes |
| 10,000 | KD tree | 0.332465 | 28,791,056 bytes | 800,000,000 bytes |

These are single-pair kernel timings, not full-package throughput or peak RAM.
Brute search wins on the small fixtures; the adopted tree does not promise a
speedup at every n. Cumulative allocation includes temporary sorting storage
and differs from peak live working memory. The quadratic-memory dense matrix
is absent from both production backends.

A subsequent integrated-package run during concurrent panel activity measured
0.826075 seconds for the 10,000-point KD case, with the same 28,791,056 allocated
bytes. This timing spread reinforces that the single-run values are workload
checks, not controlled comparative speed claims. The test counts remained
150/150 after the graph estimator generalized the shared distance helper.
