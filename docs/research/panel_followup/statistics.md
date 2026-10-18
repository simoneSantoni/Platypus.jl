# Statistics panel follow-up: density, copulas, and normalization

Reviewed 2026-10-18 on Julia 1.13.2. This is a review and diagnostic investigation;
no estimator or test-suite source was changed. Inputs are the preserved 84-dataset
audit and fresh small mathematical fixtures. Human-perception validation is out
of scope.

## Conclusion

The KDE and copula flags expose real violations of the package's mathematical
contracts. They are not evidence that the underlying KDE, Gaussian-copula, or
empirical-beta-copula definitions are unsuitable. Three distinct numerical
mechanisms interact with exact floating-point ranking:

1. Algebraic identities are not preserved by the computation, including constant
   KDE density and zero correlation on an independent Cartesian product.
2. Arithmetic expressions depend on which coordinate is evaluated first.
3. Rank transformation and density evaluation do not preserve reflection ties.

The normalizer correctly implements the mid-CDF of the *computed* values, but
there is no numerical contract connecting that CDF to the intended mathematical
scores. Consequently, tiny raw errors produce large categorical order changes.
This is a high-priority correctness problem for the normalized public output.

## S1. High: weighted constant KDE violates an exact identity

Location: `src/density.jl:113` through `:140`, particularly the final subtraction
of two independently rounded logarithmic expressions at line 140.

Three identical points `(1,1)` with weights `[1,2,3]` have both reference
coordinates constant. The documented reduced-support model drops both
coordinates. For every leave-one-out reference, density must be exactly 1 and
log density exactly 0, independent of the positive weights.

Observed log densities are `2.7755575615628914e-17`,
`1.6653345369377348e-16`, and `-1.1102230246251565e-16`; normalized scores are
`[0.5833333333333334, 0.8333333333333333, 0.25]` instead of `[0.5,0.5,0.5]`.
The estimator evaluates `logsumexp(log(p[j])) - log(sum(p[j]))`; that identity
is not exact in floating arithmetic. This explains why simply rescaling weights
or changing normalization details can change the order on a constant cloud.

Repair direction: preserve the zero-dimensional-support identity explicitly
after the existing off-support-query checks. Return log density 0 for on-support
queries, and `-Inf` for queries violating a constant coordinate. Retain the
positive remaining-reference-mass validation. A direct shortcut here is an
identity-preserving evaluation, not a statistical change or an epsilon rule.

Acceptance: arbitrary positive, unequal, rescaled weights; zero-weight duplicate
queries; off-support queries; every transform; both log/raw internal output
options. Every on-support score must be 0.5. Check at least one nonconstant
coordinate control so the shortcut cannot accidentally erase real densities.

## S2. High: feature-order arithmetic explains 91 existing flags

Locations: `src/density.jl:133`, `:137`, `:163`, and `:223`.

For a product KDE, exchanging coordinates must preserve density when the
bandwidths are exchanged with them. The audit uses automatic bandwidths, so this
condition holds. A Gaussian copula is also symmetric under exchanging its two
coordinates, with unchanged correlation. The current expressions evaluate:

```julia
log(p[j]) + lu + lv
p[i] * (zu[i]-mu) * (zv[i]-mv)
2rho * z1 * z2
```

The first two operations in each expression round before the third operand is
used. Exchanging `u` and `v` changes the intermediate result. Sorting terms before
their final sum preserves row order, but cannot make already different terms
equal. The same problem affects both estimated correlation and its density.

An isolated diagnostic module copied the current density helpers and changed
only these groupings to:

```julia
log(p[j]) + (lu + lv)
p[i] * ((zu[i]-mu) * (zv[i]-mv))
2rho * (z1 * z2)
```

It replayed all **91** existing KDE/Gaussian feature-order flagged
configurations: 13 KDE, 37 unregularized Gaussian, and 41 regularized Gaussian.
The original helpers reproduced all 91 normalized-score violations. The
isolated prototype had **zero violations and exactly zero raw differences**
between coordinate orders on every configuration. This is direct causal
evidence, not a tolerance change. The prototype does not address reflection,
constant support, zero-weight preprocessing, or moment accuracy.

Repair direction: adopt coordinate-symmetric evaluation and test both raw and
normalized equality. Also audit similarly commutative expressions outside these
helpers. The diagnostic source is dynamically copied into an isolated module;
it is not a maintained replacement estimator.

Acceptance: all 91 saved configurations plus nonuniform references, explicit
unequal bandwidths exchanged with coordinates, zero-density KDE results, and
Gaussian correlations near their supported limits.

## S3. High: independent product grids acquire spurious Gaussian ranking

Locations: `src/pairs.jl:96`, `:105` through `:124`;
`src/estimation.jl:174` through `:176`; `src/density.jl:218` through `:230`.

Use the 25 observations in `{1,2,3,4,5} × {1,2,3,4,5}`, equally weighted.
This empirical distribution factors exactly. Any deterministic marginal
transformation preserves the Cartesian product, hence its normal-score
covariance is exactly zero. Gaussian copula density is therefore 1 everywhere,
including with the explicit shrinkage 0.05 used by the probe.

The plain fit returns correlation 0 and all scores 0.5. Append two copies of
existing observations with weight zero. They must not change the reference.
The augmented fit returns correlation **3.4380156845217456e-17** and changes
normalized scores by **0.4600000000000002**.

The probe converts the *actual rounded normal scores* of the augmented fit to
`Rational{BigInt}` and recomputes their population covariance independently.
It is **exactly `0//1`**. Thus, even after allowing the marginal approximation,
the estimated nonzero dependence is a reduction-arithmetic error. Adding zero
rows also switches uniform preprocessing to the weighted rank path, whose
division and cumulative-mass operations differ from `(rank-0.5)/n`.
Coordinate-symmetric grouping alone does not establish independence accuracy.

Repair direction: reconcile uniform positive-support ranks with the unweighted
rank formula, use cancellation-aware moment evaluation, and give unresolved
near-zero covariance/score differences a documented numerical policy. A blanket
`abs(rho)<epsilon` shortcut would also erase genuine weak dependence; an error
bound, identity recognition, or higher-precision fallback needs explicit design.
Likewise, increasing precision alone reduces errors without necessarily
preventing exact-ranking amplification of the remaining errors.

Acceptance: independently factorized finite supports with uniform and product
weights, with and without zero queries; mathematically nonzero weak-correlation
controls; feature exchange and row permutation; direct checks on correlation,
log density, and final percentile. Retain the 125-point audit grid, whose existing
zero-query flag is 0.94, in addition to this smaller fixture.

## S4. High: beta-copula reflection breaks provably equal densities

Locations: `src/pairs.jl:63`; `src/density.jl:197` through `:201`;
`src/estimation.jl:90`.

For `(1,1)`, `(2,2)`, `(3,3)`, the rank coordinates are exactly intended to be
`1/6`, `1/2`, and `5/6`. Differentiating the empirical-beta mixture gives the
independent rational-polynomial oracle

```text
b(n,r,u) = n * binomial(n-1,r-1) * u^(r-1) * (1-u)^(n-r)
c(u,v)   = sum(b(n,r_j,u) * b(n,s_j,v), j=1:n) / n
```

Its exact densities are **[121/72, 9/8, 121/72]**, so expected mid-percentiles
are **[2/3, 1/6, 2/3]**. The package returns approximately
**[0.5, 1/6, 5/6]**. Reflecting the original coordinates exchanges the two
endpoint percentiles, a change of 1/3. The log densities differ by only
`2.220446049250313e-16` between orientations.

The mixture definition and untied-rank restriction agree with the primary
empirical-beta construction. Reflection follows from
`b(n,r,u) = b(n,n+1-r,1-u)` and reversing both rank vectors; it does not require
the original dependence pattern itself to be symmetric. The three-point
fixture additionally has equal endpoint densities within each fit. This proof
uses exact rational arithmetic and does not depend on `Distributions.logpdf`.
[Segers, Sibuya and Tsukahara, equations 2.1–2.2](https://arxiv.org/html/1607.04430).

Repair direction: retain integer/tie-mass rank provenance and evaluate reflected
kernel configurations canonically where possible; establish a common policy
for mathematical versus numerically unresolved ties. A polynomial oracle is
practical for tiny regression fixtures; large production fits still need stable
log-space evaluation. An arbitrary global decimal-rounding step is not a
principled fix.

Acceptance: the exact three-point oracle, reflected coordinates, exchanged
features, row permutations, and small untied permutations checked by exact
polynomials. Include genuinely unequal but close densities to detect overmerging.

## S5. Shared high-priority design gap: ranking has no error contract

Location: `src/estimation.jl:80` through `:101`.

The reference mid-CDF is mathematically correct conditional on supplied floating
scores. Therefore it should not be called an incorrectly implemented CDF. Its
exact-equality grouping nevertheless transforms tiny estimator errors at a tie
into differences proportional to the mass of that tie. Conditional CDF tests
cannot validate raw-score accuracy or the resulting scientific ordering.

Possible designs have different contracts: preserve estimator identities and
symmetries; recompute uncertain comparisons adaptively; expose interval/rank
uncertainty; or introduce documented, deterministic tie bins tied to numerical
error bounds. Approximate pairwise equality is not transitive, so a vague
`isapprox` grouping rule can make scores depend on sort order or chains of close
values. Repairs need both exact-oracle fixtures and close-but-unequal controls.

## Statistical contract review and limits

- KDE is kernel-contribution leave-one-out with full-reference ranks and
  bandwidths held fixed. Deleting an entity and rebuilding all transforms is a
  different estimator, so equivalence to a full refit is not an acceptance
  requirement. Splitting a reference weight across duplicate rows likewise
  changes which kernel contribution is excluded and changes effective sample
  size; replication invariance should not be assumed.
- Weights represent reference membership, not survey-design or frequency
  replication. Uniform multiplication should preserve the model; finite-precision
  scaling changes ratios slightly. Constant-cloud and exact-product fixtures
  avoid mistaking that caveat for an explanation of their violations.
- Gaussian density, normal-score correlation fitting, and explicit shrinkage
  are internally compatible, provided they are described as the documented
  descriptive estimator. It is not Gaussian-copula maximum likelihood. As a
  formula cross-check, the maintained R copula implementation evaluates normal
  copula density from a multivariate normal density divided by marginal normal
  densities. The two-dimensional analytic form in Platypus is algebraically
  consistent with that identity. [R copula implementation](https://github.com/cran/copula/blob/master/R/normalCopula.R).
- Rejecting constant Gaussian margins, unregularized singular dependence, beta
  ties, and nonuniform beta weights is consistent with the current contracts.
  These expected rejections should not be pooled with numerical failures.
- KDE's active-dimension convention and bivariate automatic bandwidth rule are
  declared modeling choices. A one-dimensional optimal-bandwidth discussion is
  separate from these invariance bugs. Likewise, Gaussian/beta in-sample scoring
  and KDE leave-one-out scoring need not agree.
- New smooth, data-adaptive empirical-copula estimators are legitimate later
  candidates, but the screened primary paper evaluates its estimators in its
  own Monte Carlo settings; it does not establish correctness of Platypus
  density-derived ranks or remove floating-point tie problems. Do not swap the
  estimator to hide these defects. [Kojadinovic and Yi](https://arxiv.org/abs/2106.10726).

## Procedure and artifacts

1. Read the existing audit, all 630 check flags, density helpers, weighted-rank
   preprocessing, weight validation, and normalization. Classified flags by
   method/transformation before constructing probes.
2. Revisited the original empirical-beta definition, a maintained Gaussian
   copula implementation, and the data-adaptive copula research above. Derived
   the small rational oracles independently.
3. Reproduced the constant-cloud identity and constructed a 25-point product
   grid. Proved its normal-score covariance is exactly zero even when its
   rounded normal scores are treated as exact inputs.
4. Built an isolated arithmetic-grouping prototype and replayed all 91 existing
   KDE/Gaussian feature-order failures. Recorded original/prototype raw and
   normalized differences, without modifying package source or original audit
   output.
5. Ran the probe script successfully on Julia 1.13.2. Its rational assertions
   validate the oracle premises; output CSVs retain current package violations.

Reproduce from the repository root:

```sh
julia --startup-file=no --project=. validation/consistency/panel/statistics/probes.jl
```

Artifacts in `validation/consistency/panel/statistics/` are `probes.jl`,
`environment.txt` (runtime/source hashes), `oracles.csv`, `beta_exact.csv`,
`gaussian_grid.csv`, `feature_exchange_attribution.csv`, and `summary.txt`.
The prototype result establishes this particular feature-order cause; it is not
a claim that all 630 audit flags are fixed or individually explained.
