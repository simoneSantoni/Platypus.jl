# Numerical repair: public contracts and communication

Date: 2026-10-18. Working-tree 0.2.0; numerical development runtime Julia 1.13.2.
This record covers documentation changes, not an independent numerical proof.
Human-perception validation remains outside scope.

## Status and evidence boundary

The earlier audit and renewed review are preserved as evidence about the
pre-repair implementation. Production repairs now address exact moment identities,
rank provenance, boundary decisions, reflection symmetry, support policy, and
metadata/collection integrity. That is different from the earlier isolated
91-case prototype: the present documentation describes integrated source changes.

The final frozen-source validation is complete: the 84-dataset audit passed
**14,637 checks, with zero violations and zero unexpected comparison errors**;
baseline unexpected errors were also zero. The package suite passed
**3,701/3,701** assertions, the independent mathematical comparator probe
**976/976**, and sociology checks **33/33**. The
[integration report](numerical_repair.md) records the evidence and failed
iterations. This communication record explains the resulting public contracts.

The successful audit checks weighted mid-CDF values against exposed
`normalization_scores`, rather than unconditionally against rounded physical
raw logs. Unit-invariance comparisons and tolerances are unchanged. This is an
explicit contract change, not evidence that every original test premise stayed
the same. Passing these scenarios is not a proof of universal correctness or
exact real-arithmetic values for arbitrary transcendental scores.

## Procedure

1. Read the [mathematics repair](numerical_repair_mathematics.md),
   [statistics repair](numerical_repair_statistics.md), and
   [HPC repair](numerical_repair_hpc.md), including their independent oracle
   evidence, policy decisions, costs, and primary-source records.
2. Inspect the integrated `src/depth.jl`, `numerics.jl`, `neighbors.jl`,
   `density.jl`, `pairs.jl`, `types.jl`, `estimation.jl`, and `graph.jl`
   contracts. Check metadata field names and the distinction between source
   coordinates, displayed Float64 pair coordinates, and scoring geometry.
3. Update README, changelog, manual introduction, estimation guide/API, and
   repository guidance. Parent integration owns historical supersession notices,
   the research indices, and final validation synthesis.
4. Search the assigned pages for stale spectral truncation, active KD-tree,
   implicit exactness, metadata, and completed-validation claims. Check edits
   with `git diff --check`. Numerical tests are owned by the implementing
   specialists; the root integration performs the strict documentation build.
5. Flag source-docstring omissions to the integration lead rather than editing
   source outside this task's ownership.

## Main public-contract changes

| Topic | Current documented contract |
| --- | --- |
| Neighbor membership | Outward floating bounds certify clear squared-distance comparisons; uncertain cases use exact rational fallback. Returned distances are separately rounded. |
| Backend and cost | Both requested `:kdtree` and `:brute` currently execute the conservative scan. Effective backend is `:brute`; candidate work is quadratic. No active KD-tree pruning speedup is claimed. |
| Geometry | Source-backed unweighted rank fits preserve doubled integer midranks. Raw/z-score neighbor predicates use represented transformed binary coordinates. Arbitrary decimal scaling can alter those inputs. |
| Depth | Exact rational 2×2 moments, algebraic rank, and strict affine support replace spectral truncation. Numeric `depth_rtol` raises; omit it. Thin full-rank inputs remain full rank. |
| Depth provenance | Source-backed raw/z-score depth uses original coordinates; rank depth uses exact weighted reference ranks. Affine identities require exactly preserved represented input. |
| Density/copula symmetries | Constant-support KDE and coordinate/reflection identities have explicit canonical evaluations. Rational moments preserve cancellation identities while quantiles, kernels and final arrays remain floating point. |
| Weight semantics | Original represented weight ratios reach the relevant transformations, moments, and mid-CDF masses. A new rounded input ratio is distinct from a failure to preserve the original reference. |
| Output precision | Exposed `normalization_scores` control percentile groups. Real-coordinate KDE uses ordinal labels of scale-free kernel masses; other methods use raw values. Raw physical-log ties may have distinct percentiles. No blanket epsilon merges nearby values. |
| Provenance | Ordered entity IDs bind weights to rows. Pairwise metadata retains original and normalized weights; neighbor fits distinguish requested/effective backends, and depth records rank/support policy. |
| Pair collection | Every unordered feature pair must occur exactly once; source and transform provenance must agree. |

The retained keyword `backend=:kdtree` could otherwise imply an acceleration that
is no longer active. The retained `depth_rtol` name could otherwise imply a
tuning knob that is silently ignored. Both changes therefore appear in the
README, method/API pages, and migration notes, alongside their behavior and cost.

The exact rational depth solution changes the earlier approximate rank policy.
It is not described merely as a faster or more stable implementation of the
same spectral truncation. The documentation explicitly states exact algebraic
rank, strict support, removed tolerance, and the need to refit earlier analyses.

## Real-KDE normalization contract and acceptance

The first full replay's two remaining flags were the `:none` Gaussian-KDE fit
on `zero_inflated_2026`, after common unit changes `2^200` and `2^-200`.
The density specialist's independent 512-bit calculation established a genuine
log-density difference of approximately `3.92195e-16` between rows 41 and 49.
The physical log adjustment has magnitude about 277, at which Float64 spacing
cannot retain that contrast. Calling those observations a mathematical tie, or
merging them with a tolerance, would be wrong. Demanding both exact preservation
of that distinction and a percentile obtained exclusively from the rounded
physical logs was an incompatible output contract.

The repair therefore exposes three separate matrices:

- `raw_scores`: physical or method-specific outputs, including Gaussian KDE log
  density with its bandwidth normalization.
- `normalization_scores`: the values defining order and equality for the
  reference mid-CDF. For real-coordinate KDE these are dense Float64 ordinal
  labels of dimensionless leave-one-out Gaussian kernel masses; for other methods,
  including rank-space KDE, they equal raw scores.
- `pair_scores`: weighted mid-percentiles of `normalization_scores`.

Gaussian kernel exponentials remain computed Float64 terms. Their weighted
dimensionless masses are accumulated exactly as represented, before a common
physical log normalization can erase a small difference. The ordinal labels
have no density units or meaningful spacing and cannot be compared across fits
or pairs. The initial fallback ranked all-underflow mixtures below positive
represented masses. The computational social science acceptance probe disproved
that rule: unequal weights allow a query with positive represented terms to have
less true weighted mass than a reference whose unweighted exponentials all
underflow. This is a correctness defect, not an acceptable alternative ordering.

`scoretable(result; normalization=true)` exposes the ordering; it cannot be
combined with `raw=true`. Metadata records `normalization_basis` as
`:scale_free_kernel_order` or `:raw_score`. Raw Float64 ties can intentionally
have different percentiles if ordering labels distinguish them. Public formulas
and explanations now define their mid-CDF on that exposed input.

This is a **deliberate contract change**, not a claim that all original test
premises were preserved. The independent mid-CDF audit now consumes
`normalization_scores`; unit-invariance comparisons and tolerances are unchanged.
The first full replay remains preserved in
[repair_iterations/first_full](../../validation/consistency/repair_iterations/first_full/).
The next full replay was stopped after the independent mixed-underflow fixture
blocked acceptance. The intermediate 3,671-assertion suite pass does not override
that counterexample; its logs are preserved by the integration lead. For
x=[0,100,38,139], y=0, bandwidth 1 and weights [1e-20,1,0,1e-10], query 3 received
percentile 1 instead of approximately 1e-20. The independent probe uses 2048-bit
mixture calculations and verifies a well-separated cross-regime ordering; it
does not require resolution of every exponentially tiny difference between
reference rows. See the
[mixed-underflow probe](../../validation/consistency/repair_panel/computational_social_science/mixed_underflow_probe.jl).

The density specialist implemented a replacement that compares represented
mixtures on a common basis, including previously omitted contributions within a
row. Its frozen representation is `exp(offset) * rational_mass`: ordinary rows
retain direct accumulation, while any underflowed direct kernel triggers a
shift of the whole row's **unweighted Gaussian exponents** by their maximum.
Original reference-weight ratios remain exact rational factors outside exp/log,
and the sum is divided by exact remaining reference mass. The interim weighted
log shift was rejected because it lost a resolvable conditional-weight contrast.
If a shifted exponential underflows, the implementation squares its represented
half-exponential as a rational, extending the represented range without a
floating product.

Equal offsets permit exact rational comparison; unequal offsets use adaptive,
directed MPFR bounds with explicit local precision and rounding, avoiding global
precision mutation. These comparisons certify the represented expressions, not
the exact original transcendental mixture. Even the extended shifted terms can
still round or underflow; rows with only negative-infinite Gaussian exponents
represent zero mass. The statistics specialist's final record documents 177
focused density assertions passing on Julia 1.13.2 with four threads. The
independent mathematics acceptance probe passed 976 assertions on Julia 1.10.12,
and the integration lead repeated those 976 successfully on Julia 1.13.2.
These focused results are complemented by the completed 3,701-assertion package
suite and the 14,637-check final audit, both passing. Earlier failures remain
historical evidence; they are not the final implementation's acceptance status.

Documentation changes cover README, changelog, repository guidance, introduction,
estimation guide/API and measurement-scope guide. Memory estimates now account
for three retained Float64 matrices, approximately `24*N*P` bytes. The method
guide adds an executed export example and the changelog describes migration for
users reconstructing percentiles or directly constructing result records.

## Claims deliberately avoided

- Passing the completed 84-dataset audit is evidence for those scenarios, not
  a guarantee for every input or a validation of the social construct.
- A floating filter is not approximate-neighbor retrieval: unresolved predicates
  fall back to exact arithmetic rather than widening a tie boundary.
- Exact membership does not imply exact representation of exponential affinities,
  LOF ratios, Gaussian quantiles, or every density equality.
- Exact depth preserves equal rational targets at final conversion, but two
  genuinely different depths below Float64 resolution can still round alike.
- Removing a dense distance matrix does not imply subquadratic runtime or small
  total allocations. Rational integer bit lengths and temporary storage matter.
- Metadata is an analysis record, not verification of user-declared datasets or
  temporal cutoffs. Preserve inputs, code, revision, and environment separately.

## Source and verification record

The mathematical and algorithmic decisions are grounded in the specialist
records linked above. Those records identify the primary Julia rational/
floating arithmetic documentation, LAPACK rank-policy background, empirical-beta
construction, original LOF definition, adaptive-predicate work, and self-tuning
affinity paper actually consulted during implementation. This communication pass
did not add a different algorithm or infer a new guarantee from those sources.

No source or test files were changed by this documentation task. Public examples
retain their deterministic inputs. All six estimation-guide blocks, including
the new ordering-export example, executed successfully on Julia 1.13.2 with
two threads after the contract update and again after the final kernel freeze.
The latter run used `/tmp/platypus-final-contract-docs.jl`, checked against all
six current `@example estimation` blocks (ignoring blank separators), and exited
successfully. The kernel review checked `src/density.jl` against the final
statistics record, including unweighted exponent shifts, original exact weights,
the rational half-exponential square, and local directed comparison bounds.
The integration lead identified the frozen source fingerprint as
`a60109d4a70036f7b7474d7a410e105b17f182f869b218332a45e5a0acb05f8c`.
`git diff --check` passed for the assigned files. The root's strict documentation
build passed against the final source, checking rendering, cross-references,
exported docstrings, and executable examples; the integration lead will rebuild
after these final status edits. Numerical acceptance includes the independent
controls and full before/after evidence described by the specialists, including
genuine-close-score and unequal-distance controls. Public status statements now
report the final results above and link the main integration report for history.
