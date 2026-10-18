# Numerical repair: statistical density estimators

The repair preserves the existing KDE and copula definitions. It changes their
numerical evaluation; it does not introduce an approximate-tie tolerance or
change the estimator to hide a failing invariance. Human perception remains
outside the validation scope.

## Procedure and decisions

1. Read the statistics and mathematics follow-up reports, their independent
   oracles, and the preserved simulation audit. Distinguished exact identities
   from genuinely close but different values.
2. Rechecked the empirical-beta definition against [Segers, Sibuya and
   Tsukahara, equations 2.1–2.2](https://arxiv.org/html/1607.04430), and Julia's
   [arbitrary-precision arithmetic documentation](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/).
   The beta density is the derivative of the Bernstein mixture. Its reflection
   identity allows an identical canonical computation in both orientations.
3. Implemented the mechanisms below in `src/density.jl`. Added independent
   polynomial oracles and exact identity tests in `test/numerical_density.jl`.
4. Ran the focused tests on Julia 1.13.2. All 177 focused assertions passed,
   including genuine-near-tie controls and 70 small exact-polynomial mixture
   comparisons. A wider separate sweep of 350 such mixtures also had zero
   ordering or tie disagreements. The integrated audit and complete suite
   are recorded in the overall numerical repair report.

## Root causes and implementation

**Constant support.** After validating positive leave-one-out reference mass
and checking off-support queries, two constant coordinates return log density
exactly zero. This implements the documented zero-dimensional density identity;
it avoids subtracting independently rounded logarithms of the same mass.

**Coordinate exchange.** Product log kernels group their two coordinate
contributions before adding log weight. Gaussian cross-products are grouped
symmetrically. Sorting the outer reduction alone cannot repair an asymmetric
inner expression.

**Moments.** Bandwidth variance is calculated as exact rational centered
moments of the supplied represented coordinates and original reference weights.
Gaussian covariance is calculated as exact rational centered moments of the
represented normal scores. This preserves exactly zero covariance on a finite
Cartesian product, including product weights, while retaining arbitrarily weak
nonzero represented dependence. Bandwidth square roots and squared correlation ratios are normalized by an
exact binary exponent before conversion to Float64 and square root evaluation.
This avoids intermediate overflow and underflow without depending on or changing
the caller’s BigFloat precision or rounding context. The shared numerical
helper then verifies the square-root rounding interval by exact rational
comparisons. Rational-to-Float64 conversions likewise use explicit integer
rounding, because Julia’s generic conversion may otherwise use ambient
BigFloat arithmetic.
Centering occurs before conversion, preserving represented small differences
at large offsets and avoiding overflow from squaring huge Float64 coordinates.

**Reflection.** The package passes exact weighted rank coordinates to the
estimators. Reflected Epanechnikov image distances are computed from those
coordinates, sorted canonically, and only then converted for kernel evaluation.
Complementary Gaussian ranks use one lower-half quantile approximation with an
explicit sign. Beta kernels use the lower half interval and reflect their
integer shape rank using `r -> n+1-r`. These are mathematical symmetries, not
numerical rounding rules. For real-space KDE, exact centered variance and the
sign-invariant squared coordinate difference preserve reflection.

**Uniform and zero-weight references.** Moment accumulation uses original
weights normalized by an exact rational sum; zero-weight rows cannot change
those moments. Matching exact weighted ranks in preprocessing ensures the
reference geometry does not change when zero-weight queries are appended.

## Acceptance controls and limits

The tests independently derive the three-point beta densities
`[121/72, 9/8, 121/72]`, exercise reflection in one or both margins and feature
exchange, and retain a perturbation of size `2^-40` that creates genuinely
unequal beta densities. Gaussian product-grid covariance is exactly zero for
uniform and product weights, with and without zero-weight queries. Increasing
one weight by `2^-40` must retain nonzero correlation. Constant weighted KDE,
off-support queries, nonconstant controls, large offsets, row permutations,
feature exchange, and appended zero queries are also covered.

Exact rational reductions incur allocation and integer-arithmetic costs. The
implementation uses them where cancellation or combinatorial rank provenance
matters, and retains Float64 transcendental kernel evaluation. This does not
claim an exact real-arithmetic evaluation of arbitrary transcendental scores.
Likewise, floating multiplication of weights by an arbitrary nonbinary factor
can change their represented relative ratios; that is distinguishable from an
identity failure. A universal absolute epsilon would conceal real distinctions
and is deliberately absent.

## Final replay: physical log units and distinguishable ordering

The first full repaired replay left two flags: `zero_inflated_2026`, real-space
KDE, with common binary unit changes `2^-200` and `2^200`. The 512-bit oracle
found that rows 41 and 49 have **genuinely unequal** dimensionless kernel masses,
differing by approximately `2.29158e-18`. Their log-density difference is about
`3.92195e-16`. A unit change adds a common constant to log densities; at the
resulting magnitude, one Float64 ulp exceeds that genuine difference. Both
stored raw log densities therefore round to the same number. Factoring the
bandwidth normalizer outside log-sum-exp does not solve the loss: the
scale-free Float64 logarithms also tie.

Consequently, keeping the old rule “normalize precisely the stored Float64 raw
log densities” cannot simultaneously preserve this real distinction and be
invariant to units. The repair now exposes a separate normalization input.
Real-space KDE accumulates dimensionless Gaussian kernel masses directly,
using exact rational sums of the represented exponential terms and original
reference weights. It divides by exact remaining reference mass and compares
those rational results before any conversion can merge distinct sums. Dense
Float64 ordinal labels encode the resulting order and equality; their
magnitudes have **no density units or quantitative interpretation**. The public
`normalization_scores` field exposes these labels and metadata identifies
`:scale_free_kernel_order`. `raw_scores` continues to hold log density in the
original feature units. Rank-space KDE retains `:raw_score` normalization.

A further independent CSS probe exposed an invalid initial fallback: an
underflowed unweighted kernel is not necessarily smaller than a represented
kernel multiplied by a tiny reference weight. The first mixed-underflow
ordering also lost a real conditional-weight contrast of about `1e-100` by
logging already rounded weights. Both failures are preserved in the iteration
records and corrected as follows.

Every row is represented for comparison as `exp(offset) * mass`, where `mass`
is a nonnegative rational. Ordinary rows keep their direct represented kernel
sum with offset zero. If any kernel underflows, the implementation shifts by
the maximum **unweighted** Gaussian exponent, then accumulates the shifted
kernel terms with the original exact reference weights and divides by exact
remaining weight. Weights remain outside logarithms, preserving tiny contrasts
near conditional weight one and the exact identity of a single-reference
leave-one-out fit. When a shifted exponential still underflows, its rational
value is formed by squaring the represented half-exponential. This extends the
range without a floating product or enormous exponent-sized integers.

Pairs with equal offsets use exact rational comparison. Otherwise, comparison
reduces to the sign of `log(mass_a/mass_b) - (offset_b-offset_a)`. Directed local
MPFR arithmetic encloses both the rational inputs and logarithm; precision
starts at 128 bits and doubles until the intervals separate. Every destination
has explicit precision and rounding. The implementation neither reads nor
changes the caller’s BigFloat arithmetic context. Exact equality cases are
handled first; for the remaining positive-rational ratio and nonzero rational
offset difference, equality is impossible, so refinement eventually determines
the order. This produces one strict weak order across direct and rescaled rows,
with no artificial bucket boundary or epsilon tie policy. The procedure follows
[MPFR’s directed rounding semantics](https://www.mpfr.org/mpfr-current/mpfr.html#Rounding-Modes).

The comparator certifies the order of the **represented** expressions, not the
exact original transcendental mixture. Very small shifted contributions can
still fall below the extended represented range. If every Gaussian exponent
is already `-Inf`, the represented mass is zero. These finite-precision limits
are explicit and do not justify changing an observable, resolvable ordering.

The final focused suite includes both complete four-row weighted CSS fixtures
against a 2048-bit mixture oracle, at all three unit scales; a within-row
underflowed heavy-kernel control; preservation of the tiny conditional-weight
contrast; a 400-bit close comparison that requires adaptive refinement; and
antisymmetry, transitivity, threaded execution and low ambient precision checks.
An independent mathematics acceptance probe passed 976 assertions on the same
implementation. All 177 focused density assertions pass on Julia 1.13.2 with
four threads.

The compact four-point fixture `[-2 0; 0 0; 1e-14 0; 1 0]` reproduces the issue
without depending on saved audit data. An independent 512-bit oracle confirms
the unequal scores. At both binary unit changes, the public raw logs collide
while distinct ordering labels and correct percentiles remain unchanged.
Additional tests exercise all-underflow zero-weight queries. The original
96-row counterexample also retains identical percentiles under both unit
changes after this repair.
