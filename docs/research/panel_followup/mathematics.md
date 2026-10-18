# Mathematics panel: numerical identities, support, and normalization

Review date: 2026-10-18. Julia 1.13.2. This is a review of the current working
tree, not an estimator repair. No source, package tests, or original audit
outputs were changed. The independent probe and its observations are in
[`validation/consistency/panel/mathematics/`](../../../validation/consistency/panel/mathematics/).

## Assessment

The original depth counterexamples are valid. A second, previously untested
depth defect concerns support classification: the same off-support query has
depth almost one or exactly zero depending on whether an exact invertible
coordinate change aligns a singular reference line with an axis. A blanket
approximate-score tie rule is not a correct general repair: an exact rational
control shows that it would erase genuine distinctions which the current code
correctly preserves.

Severity here describes consequences for the package's mathematical contract,
not prevalence in an application. The saved 630 flags are not 630 independent
bugs. Of these, 525 have `raw_max_relative <= 1e-10`, and 105 exceed it. Neither
count establishes how many mathematical ties exist.

| Finding | Severity | Evidence status |
|---|---|---|
| Exact equal depths become unequal percentiles | High: arbitrary distinctions in exported scores | Independently confirmed by rational covariance oracle |
| Off-support depth changes under an exact invertible shear | High: raw depth changes from 1 to 0 | New independently confirmed counterexample |
| Large common offsets degrade depth despite exactly retained data differences | High for raw-score correctness | Reconfirmed; centering probe isolates an avoidable loss |
| A universal approximate-score tie rule erases genuine close scores | High risk in proposed remediation | Independent exact control currently passes |
| Effective numerical rank and support tolerance are insufficiently described | Medium: interpretation and reproducibility gap | Confirmed source/documentation mismatch |

## 1. Exact ties and the discontinuity of empirical percentiles

For a noncollinear triangle, `C = Z'Z/3` implies each squared Mahalanobis
distance is 2. The three depths are therefore exactly 1/3 and their reference
mid-percentiles exactly 1/2. The independent probe computes the covariance,
linear solve, and depths in `Rational{BigInt}`; it does not use package internals
for the expected answer. This identity holds for any noncollinear three-point
geometry, including the rank-transformed triangle.

The package has maximum raw errors of `1.11e-16` with `:none` and `3.89e-16`
with `:rank`, but maximum percentile error **1/3**. `:zscore` also fails. The
failure has two parts:

- `src/estimation.jl:122–155` computes standardization, covariance,
  eigendecomposition and quadratic forms independently for each row. Their
  floating approximations need not preserve a mathematical equality.
- `src/estimation.jl:80–101`, particularly exact equality at line 90, correctly
  ranks the floating values it receives. This is not itself an incorrect
  empirical-CDF formula. It creates large distinctions when preceding raw
  errors split a true tie.

The comment at `src/estimation.jl:74–76` overstates the guarantee of sorting
summands: a deterministic reduction order does not generally preserve symmetry
through other arithmetic and eigensystem evaluation. It also is not an error
bound on a sum.

A weighted control strengthens the oracle: for a noncollinear triangle with
positive normalized weights `p_i`, each depth is exactly `p_i`. The fixture
with weights `[1,2,3]` satisfies that identity to floating precision and gets
the expected weighted percentiles. It should remain a passing control.

### Why global approximate ties are not a general fix

Consider four observations:

```text
(-1, -1), (1, -1), (1, 1), (-1, 1 + 2^-40)
```

The perturbation is exactly representable. Rational covariance calculations
give **two distinct exact depth groups**, separated by
`1.0105496686363572e-13`. Their exact mid-percentiles are
`[0.25, 0.75, 0.75, 0.25]`, which the current package returns correctly. A
global absolute tolerance of `1e-12` that merges these depths would return
`0.5` for all rows, introducing a new mathematical error. Similar examples
can be constructed at smaller scales.

Approximate equality also is not transitive: three scores at `0`, `0.75τ`,
and `1.5τ` have close adjacent pairs but distant endpoints. Sequential neighbor
merging can form arbitrarily wide chains; comparing to a group anchor instead
gives a different partition. Neither policy proves exact equality.

Recommended direction:

1. Preserve provable identities and symmetries in estimators, including constant
   support and analytically identifiable simplex cases, rather than smoothing
   all scores afterward.
2. Evaluate stable centered formulations and condition-aware error estimates.
   Adaptive precision can resolve some uncertain orderings, but extra precision
   alone is not a universal certificate that two computed values are exactly
   equal, particularly for transcendental estimators.
3. If uncertainty-aware score grouping is adopted, describe it as a changed
   normalization contract with explicit, deterministic grouping semantics.
   Overlapping numerical uncertainty intervals show unresolved ordering, not
   mathematical equality. Returning score intervals or unresolved-order
   diagnostics is another API choice, not a drop-in implementation detail.

Acceptance must include both the triangle's exact ties and the perturbed
quadrilateral's genuine distinctions. It must check raw values and normalized
outputs separately. A test that merely relaxes the audit tolerance, or applies
the same grouping code in its oracle, is insufficient.

## 2. New: inconsistent singular-support classification

Use references `(-1,-1)`, `(0,0)`, `(1,1)`, query `(0,δ)`, and weights
`[1,1,1,0]`, with `transform=:none`. For any nonzero `δ` the query is exactly
outside the reference affine span. Under the documented off-support convention
(`src/estimation.jl:212–214`) its depth is zero.

Apply the invertible integer shear `(x,y) -> (x,y-x)`. Its determinant is one;
for the dyadic probe values the transformed matrix is exact, verified by
rational arithmetic. The reference support becomes the x-axis, and the query
remains off support. There is no missing input information or ambiguous exact
rank in this example.

| Query δ | Original line: query depth | Sheared axis: query depth | Expected |
|---:|---:|---:|---:|
| `2^-30` | `1.0` | `0.0` | `0.0` |
| `2^-24` | `0.9999999999999987` | `0.0` | `0.0` |
| `2^-20` | `0.9999999999996589` | `0.0` | `0.0` |
| `2^-16` | `0.0` | `0.0` | `0.0` |

At `δ=2^-30`, query typicality changes from **5/6 to 0** after the shear.
The exact on-support reference depths are `[2/5,1,2/5]`.

The root cause is explicit in the code:

- A constant reference coordinate uses exact comparison to reject differing
  queries (`src/estimation.jl:116–118`).
- A singular eigendirection accepts residuals within
  `sqrt(rtol) * max(norm(z),1)` (`src/estimation.jl:150`). At default tolerance
  this introduces a tube of roughly `1e-6` in standardized coordinates and
  silently projects accepted queries onto the span.

This is a confirmed contract inconsistency, not evidence that every numerical
support tolerance is inappropriate. A tolerance can be a defensible model
choice, but it must be separated from eigenvalue rank truncation, applied
consistently to axis-aligned and oblique supports, and disclosed. A metric tube
also cannot promise unrestricted affine invariance without additional metric
semantics. Merely changing one comparison operator does not settle these
choices.

Acceptance criteria: exact on/off-support dyadic fixtures in both coordinate
representations; reflected/permuted versions; positive-reference and zero-weight
query cases; stated behavior at tolerance boundaries. If the strict support
contract remains, both queries above must return zero. If a numerical tube is
chosen, tests and documentation must explicitly describe that changed estimand.

## 3. Large-offset depth error is avoidable on the saved fixture

The five-row dyadic fixture from the original audit remains unchanged after
adding `1e14` and subtracting it again. Its exact rational depth is translation
invariant. The package's maximum raw error rises from `5.55e-17` to
**0.00564011** after translation.

`src/estimation.jl:122–126` divides raw coordinates by their large absolute
maximum before forming the mean and centered differences. This rounds away
useful precision in the small deviations. Subsequent sorted summation cannot
restore it. This is distinct from input rounding: the exact stored values
already retain all the original differences.

As a diagnostic only, the probe subtracts a represented reference row before
fitting. It verifies in rational arithmetic that this preserves the expected
depth and observes a maximum raw error of `1.11e-16`. This isolates an avoidable
centering loss; it is not a production implementation or universal guarantee.
An arbitrary first-row anchor can create row-order dependence, and subtracting
extreme opposite-signed inputs can overflow, so a production replacement needs
a deterministic safe origin and scaling strategy.

Shifted moments and corrected two-pass/pairwise algorithms are established
candidates, with different accuracy and parallel-reduction tradeoffs. The
relevant primary numerical analysis is Chan, Golub and LeVeque (1983),
[Algorithms for Computing the Sample Variance](https://www.math.pku.edu.cn/teachers/litj/notes/numer_anal/AmerStat_37_242_Chan_Variance.pdf).
The probe supports evaluating these directions; it does not benchmark or adopt
a complete replacement. The analogous preprocessing paths
`src/pairs.jl:71–79` and `129–150` also scale before centering and need coordinated
evaluation, especially when zero-weight rows switch the uniform/weighted path.

Acceptance criteria: centered/translated exact dyadic data against independent
weighted rational covariance; common binary unit changes; permutations and
reflections; uniform versus equal-positive-weight-plus-zero-query fitting;
opposite-signed extreme values without overflow; condition-aware raw-error
limits, with normalized ties assessed independently.

## 4. Numerical rank is a model/accuracy decision

`src/estimation.jl:137–150` forms a covariance of standardized coordinates and
truncates eigenvalues relative to the largest. Near this threshold, changes in
effective rank need not signify an implementation bug: truncation intentionally
changes the full-rank inverse problem. Eigenvalue cutoffs do not by themselves
define a query support tolerance. A covariance-based solve also squares the
condition number of the underlying centered weighted design matrix.

An SVD or rank-revealing QR of that design matrix is a candidate for comparison,
especially when a stable numerical rank and residual diagnostic are required.
LAPACK documents SVD/orthogonal-factorization approaches to
[rank-deficient least squares](https://www.netlib.org/lapack/lug/node27.html)
and treats thresholded singular values as regularization in its
[error-bound discussion](https://netlib.org/lapack/lug/node83.html). These sources
support algorithm screening; they do not prove that replacing a two-by-two
eigensolve alone resolves percentile instability. Centering and the score
contract remain necessary work.

Record effective rank, spectral scale/cutoff, and support-policy metadata, so
users can distinguish a numerical approximation from the full-rank estimator.
Do not demand unrestricted affine invariance from a deliberately thresholded
problem without specifying its conditioning domain.

## Audit interpretation and procedure

The audit's row/feature permutations, label changes, common reflections,
and exactly preserved rank recodings are appropriate mathematical consistency
checks. The integer-rank LOF oracle also correctly separates neighborhood
membership from subsequent floating evaluation. Agreement between two shared
floating implementations is weaker evidence than an independent oracle.

The raw threshold `abs(a-b)/max(1,abs(a),abs(b)) <= 1e-10` is a diagnostic rule,
not a universal numerical bound. Its small-raw flags should not all be called
exact ties: genuine close scores, conditioning, and floating changes in
normalized weights can also produce unstable orderings. In particular,
multiplication by `1e±150` can change represented relative weights slightly.
Identical weighted clouds still supply an exact KDE oracle for any positive
weights, but that argument must not be generalized to every weighted dataset.

Procedure performed:

1. Read the original report, full audit, reproducer, depth/normalization and
   preprocessing code, and source contracts.
2. Recomputed the original depth expectations with an independent rational
   covariance implementation; added the weighted simplex identity as a control.
3. Constructed a genuine-near-tie rational fixture to challenge approximate
   score merging rather than assume that smaller raw differences mean equality.
4. Inspected singular-support branches, constructed an exact shear with
   zero-weight off-support queries, and swept four dyadic distances.
5. Repeated the high-offset oracle and tested an exact origin shift as a
   diagnostic of the centering mechanism.
6. Ran the public-API probe successfully; all exact-arithmetic premises passed.
   Saved 16 case summaries and their per-entity inputs/observations. The probe
   records errors without pretending the unresolved package defects pass tests.
7. Recounted 630 original flags and the 525/105 raw-threshold split from CSV.
8. Consulted primary numerical-analysis sources linked above, plus Julia's
   [floating-point and arbitrary-precision manual](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/).

Reproduce from the repository root:

```sh
julia --project=. --startup-file=no validation/consistency/panel/mathematics/probe.jl
```

This panelist did not independently validate every estimator or assign a root
cause to each of the 630 flags. The highest-priority mathematical decisions are
the raw identity/ordering contract, a consistent singular-support policy, and
stable centered moments. Those should precede performance optimization of a
replacement implementation.
