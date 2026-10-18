# Changelog

## 0.2.0 — development

This entry describes the current source changes. It does not announce a
published release tag.

### Measurement and numerical behavior

- Normalize each pair's exposed `normalization_scores` to reference-weighted
  mid-percentiles. If those values are all equal, every percentile is 0.5.
  Summaries always use normalized scores;
  `normalize=false` raises an error. Inspect `raw_scores` or
  `scoretable(result; raw=true)` instead.
- Add `normalization_scores` and `scoretable(result; normalization=true)`.
  These equal raw scores except for Gaussian KDE on real coordinates: it uses
  Float64 ordinal labels of dimensionless leave-one-out kernel masses.
  Physical log densities remain in `raw_scores`. Rounded raw ties can have
  distinct normalization labels and percentiles. The labels encode order/ties
  only; they are not density values or cross-pair/fit magnitudes.
  Record `normalization_basis=:scale_free_kernel_order` or `:raw_score`.
  This explicitly changes the original mid-CDF-of-raw-Float64 contract.
- Define `novelty` explicitly as lower-tail typicality: smaller means more
  unusual. It is not the minimum, historical newness, or a significance level.
- Apply supported reference weights to transformations, fitting, and
  normalization. Zero-weight rows remain queries without contributing reference
  mass. Nonuniform weights now raise for LOF/kNN and graph scoring.
- Use filtered exact squared-distance predicates with rational fallback for
  tie-inclusive LOF/kNN and graph neighborhoods. Preserve unweighted rank-lattice
  geometry before coordinate rounding. Define duplicate handling without a
  fixed distance epsilon. kNN raw scores are negative kth distances rather than
  reciprocal distances; larger remains more typical.
- Replace rank-space Gaussian KDE with reflected Epanechnikov product kernels.
  Real-coordinate KDE retains Gaussian kernels. Return raw log densities;
  `-Inf` represents zero density. Expose bandwidths and effective reference size.
- Replace the earlier spectral depth implementation with exact rational 2×2
  moments, algebraic covariance rank, and strict affine support. Off-support
  queries receive zero. Source-backed raw/z-score depth uses original coordinates;
  rank depth uses exact reference ranks. Numeric `depth_rtol` arguments now raise.
- Preserve original represented weight ratios through exact ranks, moments, and
  mid-CDF accumulation. Preserve constant-support KDE and coordinate/reflection
  identities using explicit algebraic or canonical evaluations. Retain genuine
  close-score distinctions rather than applying a blanket tie tolerance.
  Transcendental evaluations and final result arrays remain floating point.
- Implement Gaussian normal-score copulas and empirical beta copulas. Gaussian
  singularity and endpoint queries require explicit modeling choices; beta
  copulas require untied margins and uniform positive weights.
- Implement full-feature sparse self-tuning graph scoring, including explicit
  disconnected-component stationary-mass conventions and a `max_edges` limit.
- Implement descriptive OLS with pivoted QR, optional interaction, and
  rank-deficiency rejection. No causal interpretation or default inferential
  standard errors are supplied.

### Data, performance, and visualization

- Reject nonfinite features, invalid weights, duplicate IDs, duplicate feature
  names, and invalid neighborhood parameters.
- Make `EntitySet.matrix` authoritative, return accessor snapshots, and retain
  detached raw-data snapshots for weighted pair fitting.
- Replace dense distance matrices with filtered exact candidate scans. Both
  `:kdtree` and `:brute` requests currently use `effective_backend=:brute`:
  uncertified KD-tree pruning is disabled. Candidate work is quadratic; exact
  arithmetic adds allocation and integer-bit costs. Support independent-pair
  threading without promising a KD-tree speedup.
- Share transformed feature vectors among pairs. Score matrices still require
  storage proportional to entities × feature pairs, and worst-case neighbor
  search remains expensive.
- Add `analysis_metadata`, raw score exports, and detached configuration.
  Bind reference weights to ordered entity IDs, retain original represented
  weights for pairwise fits, distinguish requested/effective backends, and record
  geometry, normalization basis, and depth rank/support policies. Pairwise results
  retain three Float64 score matrices (approximately `24*N*P` bytes), before
  metadata and temporary exact-arithmetic storage.
- Validate that each `PairSet` contains every unordered feature pair exactly
  once, with consistent source/transform provenance; reject duplicate or missing
  pairs rather than silently changing the aggregation.
- Make `pair_redundancy` a diagnostic-only correlation table: remove the
  unjustified `weight` column and retain undefined correlations as `NaN`.
- Implement all nine plot helpers as detached `TypicalityPlot` RecipesBase
  objects. Rendering requires a separately installed plotting frontend such as
  Plots: `plot(plot_novelty_conventionality(result))`.

### Documentation and validation

- Set Julia 1.13.2 as the development, primary CI and documentation-build default;
  record the runtime in `.julia-version` and document checkout-local Juliaup setup.
  Retain the Julia 1.9 compatibility floor and minimum-version CI coverage.

- Replace unseeded tutorial claims with deterministic, executed examples and
  assertions. Explain score direction, support assumptions, memory costs,
  data selection, and reproducibility.
- Separate descriptive measurement from null-model tests, historical novelty,
  causal effects, and unperformed human-perception validation.
- Add algorithm research records, numerical/property tests, synthetic validation,
  and reproducible benchmarks.
- Remove development-as-stable documentation routing, use a portable docs
  environment, and make local documentation builds separate from deployment.

### Numerical repair validation status

The final 84-dataset numerical audit passed **14,637 checks with zero violations
and zero unexpected comparison errors**; baseline unexpected errors were also
zero. The package suite passed 3,701/3,701 assertions, independent mathematical
comparator checks 976/976, and sociology checks 33/33. The
[numerical repair report](docs/research/numerical_repair.md) records the final
source fingerprint, evidence, and preserved failed iterations.

Real-coordinate KDE represents each row as an exponential offset times rational
mass. Any underflowed direct kernel triggers an unweighted-exponent shift, with
original reference-weight ratios retained exactly outside the exponential; a
rational half-exponential square extends the represented range. Cross-offset
comparisons use adaptive local directed bounds without changing global precision
state. Independent controls cover mixed underflow with unequal weights and
genuine close-density contrasts across unit changes.

The independent mid-CDF check now consumes exposed `normalization_scores`;
unit-invariance checks and tolerances are unchanged. This is an explicit public
contract change. The comparator certifies ordering of represented expressions;
it does not eliminate further shifted-term underflow or certify exact Gaussian
densities. Passing these scenarios is not a universal all-input correctness claim.

### Migration from earlier 0.2 development snapshots

Omit `depth_rtol`; numerical covariance-rank truncation is no longer available.
Thin full-rank data is treated as algebraically full rank. Refit analyses that
used the earlier spectral policy. Both requested neighbor backends now scan
candidates, so remeasure runtime and memory before scaling up. Preserve the
new ID-bound weight mapping and requested/effective backend fields with exports.

Code that reconstructs percentiles must now use `normalization_scores` and the
recorded original `source_reference_weights`, not `raw_scores` unconditionally. Include the new
ordering matrix when directly constructing or serializing `TypicalityResult`;
prefer estimator-returned results. Export raw, ordering, and normalized columns
with their distinct meanings. Raw and normalization export switches cannot both
be true.

### Migration from 0.1

Recompute and revalidate old scores: normalization, LOF ties, KDE kernels, depth
degeneracy handling, and weighted references have changed. Use `raw_scores`
rather than `normalize=false`; remove uses of `pair_redundancy(...).weight`.
Treat plotting helpers as recipe constructors. Rebuild pairs after data edits
and save `analysis_metadata` with exports. Previously reserved estimators are
implemented, with explicitly documented assumptions and distinct result types.

## 0.1.0

Initial pairwise data model, exploratory helpers, and four typicality scoring
rules. Additional estimators and all visualization functions were exported
placeholders.
