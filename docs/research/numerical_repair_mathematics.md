# Mathematical repair: exact pairwise Mahalanobis depth

The previous depth implementation lost information while centering large-offset
data, split exact mathematical score ties, and used incompatible tests of
singular support. A floating eigensolve could improve conditioning but could not
certify the equalities needed by exact mid-percentile normalization. The repair
uses the special structure of the problem: every fitted covariance is 2 by 2.

## Algorithm and policy decision

Represent supplied finite binary coordinates and reference weights as exact
`Rational{BigInt}` values. Exact weighted moments produce

```
C = [a b; b c], determinant = a*c - b*b
q = (c*x² - 2*b*x*y + a*y²) / determinant
depth = determinant / (determinant + c*x² - 2*b*x*y + a*y²)
```

Here `x,y` are exactly centered query coordinates. Convert the final rational
depth to `Float64` once. Equal mathematical depths then have identical stored
values. No tolerance merges nearby distinct scores. Differences smaller than
the final floating output's resolution can still round to the same stored value;
the public result remains a floating array.

If the exact determinant is zero, determine whether the reference support is a
line or point. For a line with `a > 0`, the exact test `a*y == b*x` decides support
and on-support depth is `a/(a+x²)`. If only `c > 0`, require `x == 0` and use the
analogous expression. For a point, both centered coordinates must be exactly
zero. Every exactly off-support query receives zero, including queries less
than a floating eigensolve's former support tolerance from the line.

The repaired default uses **exact algebraic covariance rank**, including very
thin but full-rank data. This deliberately removes the previous approximate
spectral truncation policy. Numeric `depth_rtol` arguments now fail with an
actionable message; omit that argument. Silently ignoring a supplied tuning
parameter would conceal a changed estimator. Metadata records the exact
effective rank, exact support policy, solver and absence of rank truncation.

Depth is affine invariant when the transformation preserves the represented
input exactly. The fit uses original source coordinates for `:none` and
`:zscore`, because an invertible coordinate standardization does not change
Mahalanobis depth. Rank fits use exact weighted rank coordinates. This also
avoids treating floating rounding of an intermediate transform as observed
structure. The solver cannot undo rounding introduced before data reaches the
package; arbitrary floating coordinate transformations are not exact identities.

## Algorithm screening and online sources

Reviewed the prior panel's floating covariance/eigendecomposition, shifted
moment and SVD recommendations. SVD/QR are established methods for rank-deficient
least squares, but their numerical rank decisions do not establish exact score
equality or strict affine support. See the primary
[LAPACK least-squares documentation](https://www.netlib.org/lapack/lug/node27.html).
The Julia manual documents rational types and their arithmetic; the repair uses
arbitrary-size integer numerators/denominators to avoid integer overflow. See
[Julia rational numbers](https://docs.julialang.org/en/v1/manual/complex-and-rational-numbers/#Rational-Numbers).
Both sources were opened during this repair. The closed 2-dimensional rational
solve is an implementation choice derived from covariance algebra, not a claim
that the sources prescribe this package's estimator.

The exact solve was selected because it establishes the required identities
directly and remains practical for this audit's small pairwise fits. A local
Julia 1.13.2 probe, after warming compilation and averaging ten calls, measured
approximately 0.46 ms for 32 rows, 3.6 ms for 125 rows and 20 ms for 1,000 rows
with random floating coordinates and weights. These are single-machine timing
observations, not performance guarantees or a comparative benchmark. Exact
arithmetic has larger allocation and integer-bit complexity than a floating
solve; very large data or high-dimensional pair sets need separate measurement.

## Independent acceptance evidence

`test/numerical_depth.jl` uses a separate rational matrix covariance and generic
linear solve as an oracle, plus fixtures with analytic answers:

- Three noncollinear equally weighted observations: all depths are exactly
  one third and all normalized scores are one half.
- A perturbed quadrilateral: two genuine depth groups about `1e-13` apart are
  retained, with normalized scores `[0.25,0.75,0.75,0.25]`.
- Exactly represented large offsets, common binary unit changes, row and
  feature permutations and reflections preserve the depth values.
- Exactly off-line zero-weight queries remain zero after an invertible shear.
- Opposite-signed extreme floating values do not overflow intermediate moments.
- Tiny but full-rank triangles remain rank two; numeric truncation arguments
  are explicitly rejected.
- Adding zero-weight copies leaves the existing reference scores unchanged.

The full integrated test and simulation results are recorded by the repair's
parent validation procedure. Original investigation outputs are retained as
historical evidence rather than overwritten with repaired results.
