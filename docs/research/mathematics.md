# Mathematical audit and independent oracles

Date: 2026-10-18. Scope: definitions, invariances, singular cases, numerical
arithmetic, and the information retained by pairwise projections. Human
perception is outside this work. The checks establish agreement with explicit
mathematical definitions; they do not establish social construct validity.

## Research procedure and decisions

The initial audit reproduced two failures: multiplying all membership weights
by a positive constant changed KDE rankings, and changing one coordinate's
units reversed Mahalanobis-depth rankings under `transform=:none`. The causes
were dimensionful additive epsilons. The revised implementation normalizes
membership weights before arithmetic, handles empty support explicitly, and
solves depth in dimensionless coordinates with a declared singularity tolerance.

Primary-source searches used “Segers Sibuya Tsukahara empirical beta copula”,
“Zuo Serfling general notions statistical depth affine invariance”, and
“Breunig LOF k-distance neighborhood”. The original papers were inspected,
including the beta copula's equations 2.1–2.2 and its discussion of ties. This
was a targeted mathematical audit, not a systematic review or a claim that the
selected methods dominate alternatives.

| Candidate or mathematical property | Decision and reason |
|---|---|
| Original tie-inclusive LOF | Retain the original distance-neighborhood definition. Check row permutations against an independent small dense oracle. |
| Mahalanobis depth | Retain as global centrality. Solve covariance geometry without a fixed length-scale ridge; declare how singular support is treated. |
| Gaussian copula | Adopt as an explicit parametric comparator. Independence has constant density one, which supplies an exact oracle. |
| Empirical beta copula | Adopt only with equal positive weights and untied marginal ranks. These restrictions preserve its finite-sample copula identity. |
| Reflected compact KDE | Adopt as a density on bounded rank coordinates, without claiming its estimated margins are exactly uniform. |
| Higher-order-only atypicality | Document as outside the information available to a pairwise method, rather than inventing a correction inside aggregation. |

Zuo and Serfling characterize depth using affine invariance, central maximality,
center-outward monotonicity, and decay at infinity. That is a centrality target,
not a local-density definition. For example, a uniform independent square has
constant copula density while its Mahalanobis depth decreases away from its
center. Neither behavior makes the other method incorrect; the estimands
differ. [Original depth paper](https://www.stt.msu.edu/~zuo/papers_html/annals001.pdf).

The empirical beta copula is an equal mixture indexed by the original marginal
ranks. Under untied ranks, every rank occurs exactly once in each margin, so
the mixture has uniform margins. Arbitrary mixture weights or repeated ranks
break that identity. The package rejects those inputs and evaluates the full
in-sample mixture; a leave-one-component-out mixture would require a different
claim. [Segers, Sibuya and Tsukahara, section 2.1](https://arxiv.org/html/1607.04430).

The original LOF neighborhood includes every observation within the kth
nonself distance. Keeping exactly k observations after sorting is not the same
definition when distances tie. The package now includes all boundary ties and
declares its extension for duplicate clusters with zero reachability distance.
[Breunig et al., section 4](https://webdocs.cs.ualberta.ca/~zaiane/pub/check/breunig.pdf).

## Algebraic oracles

The following are direct derivations used to check the implementation.

**Membership weights.** Normalize finite nonnegative weights to probabilities
`p`. For a common positive scale `a`, `aw/sum(aw)=w/sum(w)`. Consequently,
transforms, fitted parameters, raw scores, and final normalization must be
unchanged apart from floating-point rounding. Zero-weight rows are queries;
appending them must not alter existing scores. A leave-one-out denominator is
the remaining positive mass, not an epsilon-regularized approximation to it.

**Mid-percentiles.** For raw scores `s`, compute

```text
T_i = sum_j p_j * (indicator(s_j<s_i) + indicator(s_j==s_i)/2).
```

This definition gives exactly one half for a tied reference distribution.
The test oracle uses a direct quadratic comparison loop, independently of the
production sorting implementation. Median, lower quantile, and IQR are then
summaries of `T`, while raw scores remain separately available. Floating-point
noise can split mathematically tied scores; ordering and reductions must be
tested explicitly because rank normalization can magnify that noise.

**Nonsingular depth.** For weighted mean `mu` and population covariance `C`,
the reference calculation is `1/(1+(x-mu)'*(C\(x-mu)))`. For any invertible
matrix `A`, covariance transforms to `A*C*A'`; substitution proves invariance
under a nonsingular affine coordinate change. Tests compare the production
spectral calculation with a direct covariance solve, coordinate rescalings,
reflections, translations, and a shear. At numerical rank boundaries, the
declared spectral tolerance determines which support is resolved; this is not
a promise of exact affine invariance under arbitrarily ill-conditioned maps.

**Singular depth.** A reference supported at `(-1,-2),(0,0),(1,2)` has intrinsic
coordinate variance `2/3`. Its depths are `.4,1,.4`; query `(2,4)` has depth
`1/7`. Query `(0,1)` is outside the affine span and receives zero by the
package's explicit extension. A point mass has depth one at its support and
zero elsewhere. These conventions avoid treating off-support observations as
central merely because a pseudoinverse ignores null-space coordinates.

**Copula checks.** Normal-score ranks `[1,2,3,4]` and `[2,4,1,3]` have cancelling
cross products and zero correlation. The Gaussian copula therefore has log
density zero at every point. For two concordant empirical beta ranks,
`c(u,v)=2*((1-u)*(1-v)+u*v)`; at the two midrank observations its density is
exactly `1.25`. This polynomial oracle does not call a beta-distribution density
implementation. Singular Gaussian dependence requires explicit shrinkage,
whereas ties and nonuniform weights are rejected for the beta family.

**Bounded density.** The reflected Epanechnikov kernel conserves mass on the
unit interval when the bandwidth is at most one. Reflection is symmetric in
evaluation and center coordinates, so averaging it over a continuous uniform
reference also gives density one. A dense Cartesian midpoint grid approximates
this integral. Its finite grid and omitted self contribution introduce real
errors; the test checks raw density near one with an explicit tolerance instead
of demanding equal finite-sample ranks.

## Information-loss counterexample

Take 100 copies of each binary triple `000,011,101,110`, and one copy of each
`001,010,100,111`. Every bivariate cell contains 101 observations. With k=20,
kNN and LOF cannot distinguish the four rare triple combinations from common
ones: every projected point has 100 coincident neighbors. The test asserts this
blind spot to prevent documentation from implying general joint detection.
This finite binary construction illustrates information loss; it is not an
example satisfying the untied continuous assumptions of the beta copula.

The limitation also holds for continuous distributions. On the unit cube,
`f(x,y,z)=1+theta*(2x-1)*(2y-1)*(2z-1)`, with `abs(theta)<=1`, is a nonnegative
density. Integrating out any coordinate removes the product term, leaving
every bivariate marginal equal to one. Thus different values of theta have
identical bivariate distributions despite different three-way dependence.
This is a direct integration argument, not an empirical performance claim.

## Execution record

The public checks are in `test/mathematical_properties.jl`. The package test
entrypoint includes this file; numerical density and neighbor kernels also have
separate independent suites. The public mathematical suite passed **137/137**
checks on Julia 1.10.12, using `--project=. --compiled-modules=no`, on 2026-10-18.
The initial run exposed normalization roundoff in the all-tied 404-row case;
dividing the weighted mid-CDF by its actual accumulated total restored exactly
one half. An additional adversarial weight vector `[1; fill(1e-16,9)]` exposed
row-order-dependent floating reductions. Sorting reductions and directly
detecting constant reference support fixed the reproduced cases. Those
regressions are included in the passing count. The integration validation
record provides the full package totals. No inferential calibration or
agreement with human judgment is implied by passing these tests.
