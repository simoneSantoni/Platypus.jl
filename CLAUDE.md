# Repository guidance

## Package and scope

Platypus.jl measures descriptive entity atypicality relative to a supplied
collection and selected feature schema. The working version is 0.2.0 (unreleased).
Pairwise methods score all feature pairs and aggregate normalized typicality;
a separate full-feature graph estimator measures connectivity. The median and
lower-tail summaries are inspired by Uzzi et al. (2013); IQR is a package
extension. This is not a reproduction of their randomized co-citation model.
Human-perception validation is outside the current correctness phase.

## Development commands

Use Julia **1.13.2** by default, recorded in `.julia-version` and pinned in CI
and documentation builds. Configure a checkout once with
`juliaup add 1.13.2` and `juliaup override set 1.13.2`, or use
`julia +1.13.2` explicitly. The package retains Julia 1.9 compatibility and a
minimum-version CI job. Historical validation records keep their actual runtime
versions; do not relabel older results as runs on the new default.
Use version-specific local manifests when switching runtime series. In this
checkout, `Manifest-v1.10.toml` preserves the previous environment and
`Manifest-v1.13.toml` records the new one. These generated files are ignored by
git. Reusing the old generic manifest under Julia 1.13 caused incompatible
standard-library pins and precompilation failures.

```bash
julia --project -e 'using Pkg; Pkg.test()'
julia --threads=2 --project -e 'include("test/runtests.jl")'
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'
julia --project validation/experiments.jl
julia --project benchmark/neighbors.jl
julia --project benchmark/graph.jl
```

Tests are organized into focused files included by `test/runtests.jl`.
Documentation builds are strict, execute examples, and check exported docstrings.
Local builds do not deploy; deployment requires `PLATYPUS_DEPLOY_DOCS=true`.
New manual pages must be registered in `docs/make.jl`.

## Architecture and contracts

- `src/types.jl`, `entities.jl`: finite validated inputs, unique IDs and feature
  names, one canonical `EntitySet.matrix`. Accessors and entity iteration return
  detached snapshots. Direct mutation of public storage remains possible and is
  validated at downstream entry points; do not promise immutable containers.
- `src/pairs.jl`: rank (default), z-score or raw transformations; shared feature
  coordinate vectors reduce pair storage. `PairSet` retains source data and
  transformation provenance. Exact rational weighted ranks and centered moments
  preserve represented reference ratios and zero-weight exclusion. Source-backed
  pair collections must contain every unordered pair exactly once, with consistent
  source/transform metadata; manual legacy coordinates retain unknown provenance.
- `src/eda.jl`: descriptive summaries, correlations and nearest entities.
  `pair_redundancy` reports absolute correlation, not valid aggregation weights.
- `src/neighbors.jl`: outward floating filters plus exact rational fallback for
  squared-distance membership, tie-inclusive LOF, explicit duplicate handling.
  BOTH requested backends currently execute conservative candidate scans:
  requested `:kdtree` does not enable KD-tree pruning. Candidate work is quadratic;
  scratch storage avoids a dense distance matrix. Rank geometry uses doubled
  integer midranks; raw/z-score predicates use transformed binary coordinates.
- `src/numerics.jl`: exact rational square-root rounding decisions without mutable
  global BigFloat precision or rounding contexts.
- `src/depth.jl`: exact rational 2×2 covariance, algebraic rank and strict affine
  support. No spectral truncation or tolerance tube. Numeric `depth_rtol` is
  rejected; omit it. Source-backed raw/z-score depth uses original source
  coordinates; rank depth uses exact reference ranks.
- `src/density.jl`: reflected Epanechnikov rank-space density and Gaussian
  real-space density; kernel-contribution leave-one-out holds reference marginal
  transformations and bandwidth fixed. Gaussian normal-score and empirical beta
  copulas have separate assumptions and are not interchangeable with KDE.
  Constant-support KDE and coordinate/reflection symmetries have explicit
  identity-preserving evaluations. Rational centered moments preserve cancellation
  identities; normal quantiles, kernels, exponentials and final arrays remain
  floating values. Do not claim every transcendental tie is certified.
- `src/estimation.jl`: `fit_typicality` accepts the closed set
  `TYPICALITY_METHODS = (:knn, :lof, :kde, :depth)`; `fit_copula` is separate.
  Always aggregate reference-weighted mid-percentiles. Larger pair scores mean
  more typical; `novelty` is a lower quantile, so smaller means more unusual.
  Conventionality is the median, dispersion the IQR, atypicality one minus the
  median. `normalize=false` is rejected; use `raw_scores` or
  `scoretable(result; raw=true)` for method-specific scores.
  The mid-CDF consumes exposed `normalization_scores`. These equal raw scores
  except real-coordinate KDE, which supplies dense ordinal labels of scale-free
  kernel masses. Raw physical log-density ties can intentionally have different
  percentiles. Labels have no density units or cross-pair/fit magnitude meaning.
  Export them with `scoretable(result; normalization=true)`; this option and
  `raw=true` are mutually exclusive. Record `normalization_basis` as
  `:raw_score` or `:scale_free_kernel_order`.
- `src/graph.jl`: sparse, exact union-kNN self-tuning affinity graph. Explicit
  stationary-distribution selection for disconnected graphs, edge-count limit,
  separate `GraphTypicalityResult`; no claim that connectivity is density.
- `src/outcome.jl`: descriptive OLS with pivoted QR and rank checks. No inference,
  standard errors, temporal validation or causal interpretation is implied.
- `src/visualization.jl`: all plot functions construct `TypicalityPlot` recipes.
  Rendering is optional: `using Plots; plot(plot_pair_cloud(ps, 1))`.

## Scientific and engineering conventions

- The supplied collection is the endogenous reference. Do not silently substitute
  an external reference or claim a frozen historical query interface exists.
- Labels are descriptive and must never alter estimation.
- Keep `:rank` the default, but describe its removal of marginal magnitude and
  its limitations. Pairwise representations do not identify every joint pattern.
- `ReferenceWeights` means graded reference membership, not survey, inverse
  probability, frequency-replication or outcome regression weights. Support
  nonuniform weights only where the full estimator contract supports them;
  reject unsupported combinations rather than ignoring weights.
- Percentiles are descriptive, not null probabilities. IQR across dependent
  pairs is not uncertainty. Cross-cohort scores require explicit interpretation.
- Preserve raw scores and detached metadata, including ordered entity IDs,
  normalized and original represented weights, feature schema, package version,
  geometry, requested/effective backend, depth rank/support and provenance.
  Numerical policy is part of the recorded analysis, not an implicit detail.
- Prioritize independent rational/polynomial oracles and regression controls over
  speed. Agreement between implementations that share a predicate is insufficient.
  Preserve true identities AND genuinely close-but-unequal controls; never widen
  tolerances or merge nearby scores merely to pass an audit.
- Distinguish exact represented-input geometry from an ideal transform lattice.
  Arbitrary decimal rescaling may alter represented data. Strong invariance
  assertions require that the input transformation preserve the intended geometry.
- Candidate accelerators need certified pruning or conservative fallback.
  Deterministic reduction order is not a proof of exact arithmetic or symmetry.
  Document singular cases, bandwidth assumptions, and finite-output limits.
- Record research sources, comparisons, adoption decisions and executed checks
  under `docs/research/`; keep failed simulations visible. The deterministic
  validation suite is evidence about its scenarios, not universal superiority.
- Do not restore exported stubs or misleading unsupported claims. API changes in
  this development cycle are documented in `CHANGELOG.md`.
- The final 84-dataset audit passed 14,637 checks with zero violations and zero
  unexpected comparison errors; baseline unexpected errors were zero. Package
  assertions passed 3,701/3,701, independent mathematical comparator checks
  976/976, and sociology checks 33/33. The real-KDE normalization ordering is a
  deliberate public-contract change. Preserve that distinction when reporting
  successful validation; see `docs/research/numerical_repair.md` for evidence.
  Never rank all-underflow mixtures categorically below positive represented
  masses: unequal weights can reverse that ordering.
  If any direct kernel underflows, shift unweighted kernel exponents for the
  entire row, keeping original exact weight ratios outside exp/log. A rational
  square of the represented half-exponential extends the range when needed.
  Compare exp(offset)*rationalmass expressions using exact equal-offset decisions
  or adaptive local directed MPFR bounds; no global precision mutation. This
  certifies represented expressions, not exact original Gaussian mixtures;
  shifted terms may still underflow. Keep first_full, interrupted runs, original
  failures, final results, focused tests and integrated validation distinct.
- Audit the mid-CDF against exposed `normalization_scores`, not raw logs for
  real-coordinate KDE. State this changed oracle premise explicitly; do not
  claim every test contract stayed unchanged. Unit-invariance comparisons and
  tolerances are unchanged. Direct Gaussian terms remain Float64; exact sums do
  not remove exponential rounding or underflow limitations.
