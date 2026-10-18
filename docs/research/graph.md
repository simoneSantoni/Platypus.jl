# Full-feature graph estimator: research and implementation record

Date: 2026-10-18. This extends the HPC review with a mathematically explicit
implementation of the previously advertised graph entrypoint. The objective
is internal consistency and an alternative representation capable of retaining
joint feature information. Human perception is outside scope.

## Screening and selection

Online searches covered “self-tuning spectral clustering Zelnik-Manor Perona”
and “Lovasz random walks graphs stationary distribution degree”. The primary
self-tuning paper defines an affinity with observation-specific length scales:
`exp(-distance(i,j)^2 / (sigma[i]*sigma[j]))`, where scales can be selected from
neighbor distances. This motivates local scaling in heterogeneous collections;
we adopt the affinity, not the paper's clustering algorithm.
[Zelnik-Manor and Perona, 2004](https://papers.neurips.cc/paper/2619-self-tuning-spectral-clustering.pdf).

For a symmetric weighted graph, degree divided by total degree is stationary
for the ordinary random walk. Disconnection makes stationarity nonunique;
bipartiteness can prevent convergence from arbitrary starts without invalidating
the stationary vector. These distinctions allow direct calculation without
eigenvectors or a convergence-dependent iterative approximation.
[Lovász, 1993](https://www.cs.cmu.edu/afs/cs/academic/class/15859n-f16/RelatedWork/random-walks-on-graphs.pdf).

| Option | Decision and reason |
| --- | --- |
| Fixed global Gaussian scale | Not selected; a single scale adds an unvalidated tuning choice across heterogeneous feature clouds |
| Local kth-distance scales | Adopted; local neighborhoods and scale derive from the same explicit k |
| Directed kNN graph | Not selected; stationary behavior depends on directed recurrent classes |
| Union-symmetrized kNN graph | Adopted; preserves either endpoint's inclusion and gives a transparent reversible walk |
| Mutual kNN graph | Deferred; rejects one-sided neighbors and increases disconnection |
| PageRank teleportation | Considered, then omitted; uniqueness is purchased by a new mixing parameter and altered estimand |
| Spectral eigensolver | Unnecessary for degree-proportional stationarity; adds numerical and convergence concerns |

These decisions are package design choices. The cited papers do not establish
that graph connectivity is universally preferable to density or pairwise
unusualness. Self-tuning kernels can equalize regions with different density;
their graph score must therefore be described as connectivity.

## Implemented definition

1. Validate the full `EntitySet`, transform every selected feature, and retain
   all dimensions in the distance calculation. Default transform is empirical
   marginal mid-rank; `:none` and `:zscore` follow the same preprocessing contract
   as pairwise scoring.
2. Find exact nonself kth distances, including every boundary tie. Reuse the
   bounded-memory KD-tree/brute kernels documented in [HPC notes](hpc.md).
3. Form the union of directed neighborhoods, assign symmetric self-tuning
   affinities, and retain a sparse adjacency matrix. A zero scale gives affinity
   one to coincident points and zero to noncoincident points, an explicit
   zero-scale extension rather than an arbitrary numerical epsilon.
4. Calculate weighted degrees with sorted reductions. Normalize degrees to
   obtain stationary mass; use empirical mid-percentiles for graph typicality
   and one minus these percentiles for atypicality.
5. Record feature names, transformation, effective and requested k, actual
   backend, affinity family, symmetrization, component count, stationary
   selection, isolation convention, edge limit, reference weights, and caller
   provenance. Copy input identifiers and labels into the fitted result.

If edges exist, component masses are proportional to their total weighted
degree, and isolated vertices have zero mass. An isolated vertex's conceptual
transition is a self-loop; adjacency and degree remain zero, so these fields
stay consistent. If the entire graph is edgeless, the transition is identity
and the selected stationary distribution is uniform. The result explicitly
reports whether stationarity is unique.

Nonuniform reference weights are rejected: a weighted transition/reference
model would require a separate derivation. The result is `GraphTypicalityResult`,
with no invented pair scores, pairwise conventionality, novelty, or dispersion.

## Memory and validation

Exact boundary ties can turn a kNN graph into a complete graph. The default
`max_edges=1_000_000` bounds stored undirected edges; exceeding it raises an
actionable error. It does not silently drop ties. Sparse storage is O(nk) only
when ties are limited. High-dimensional exact search may approach brute-force
time, so this is not a claim to eliminate the curse of dimensionality.

Validation includes an independent dense affinity oracle, direct stationarity
checks using a constructed transition matrix, KD/brute equivalence, exact row
permutation invariance, symmetric square and constant clouds, disconnected
components, duplicate clusters plus an isolated vertex, edge-limit rejection,
uniform weight rescaling, monotone-feature rank invariance, and detached result
metadata. A three-feature parity construction confirms that this full-feature
representation can distinguish configurations that pairwise LOF cannot in that
fixture; this is a representation check, not a general detection guarantee.

Integrated execution: **219/219 graph checks passed** on Julia 1.10.12 with
NearestNeighbors 0.4.29. After generalizing the common distance helper for
full-feature graphs, the existing **150/150 neighbor checks also passed**.
Command: `julia --project -e 'using Platypus, Test, DataFrames;
include("test/neighbors.jl"); include("test/graph.jl")'`.

A warm single-thread execution on the reported `goldmont` CPU, 1,000 observations,
eight uniform seeded features and k=20 produced 12,745 undirected edges in one
component. It took 0.078666 seconds, allocated 6,447,408 bytes cumulatively, and
returned a 472,648-byte result. Reproduce this workload with
`julia --project benchmark/graph.jl`. This is a workload check, not a controlled
cross-hardware or cross-algorithm speed claim.
