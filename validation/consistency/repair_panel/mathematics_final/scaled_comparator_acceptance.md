# Independent acceptance of scaled kernel ordering

Accepted implementation: `src/density.jl` SHA256
`67f545b13c28e0716b6115d357ce92fa10e96aafba5a673801e101a3a9a29bca`.

Execution: **976/976 checks passed**, Julia 1.10.12, four threads, 2.3 seconds.
The probe includes the source helpers without importing the package test suite.
It does not modify source. See `scaled_comparator_probe.jl`,
`scaled_comparator_execution.log`, `scaled_comparator_environment.txt`, and
`scaled_comparator_mixed_oracle.csv` alongside this note.

## Mathematical comparator review

Each key denotes the positive real `exp(offset)*mass`, where the represented
offset is a finite rational and mass is an exact nonnegative rational. Zero
masses compare as zero independently of their offsets. Equal offsets reduce
comparison to exact rational mass comparison. Otherwise the comparison is
equivalent to comparing `log(mass_a/mass_b)` with `offset_b-offset_a`.

The implementation encloses rational inputs using separately rounded BigInt
numerators and denominators, crossing denominator bounds appropriately. For
negative rationals the denominator directions are reversed. Directed MPFR
division and logarithm then produce enclosing intervals. A result is returned
only after intervals separate. The exactly equal case is handled before the
loop: for a positive rational ratio and rational target, `log(ratio)=target`
requires ratio one and target zero. Otherwise adaptive refinement eventually
separates the represented values. Thus comparison defines one strict order
across ordinary and shifted keys, rather than treating underflow as a lower
bucket.

All MPFR destinations have explicit local precision and all relevant operations
have explicit rounding directions. The code avoids Julia's rational BigFloat
constructor, which internally changes global precision/rounding on Julia 1.10.
Exact-zero constructors do read the default rounding argument, but their value
is unaffected; the comparison does not depend on or mutate ambient settings.

Independent tests verify signed rational enclosures with magnitudes through
`2^±4096`; dyadic brackets around e at 400 fractional bits, forcing refinement
beyond the initial 128 bits; unchanged ambient precision and rounding under
64/256-bit settings and three rounding modes; agreement with a 2048-bit oracle
for mixed representations; antisymmetry, transitivity, equality, row permutation,
and concurrent comparisons.

## Kernel construction and regression evidence

Underflow-affected rows now shift **unweighted kernel exponents** and multiply
the resulting represented exponentials by exact original reference weights.
The remaining leave-one-out mass is divided exactly. This matters: logging
tiny conditional weights can turn a real `1-1e-100` contrast into exactly one,
even when the scaled comparator itself is correct.

The original limited 972-check run passed the comparator and separated-query
criteria. Expanding acceptance to every row of both weighted mixed-underflow
fixtures exposed two failures in the conditional-log-weight construction:
references 2 and 4 incorrectly tied. The final construction passes all 976
checks. Complete eight-row results and exact weighted-CDF comparisons are
recorded in `scaled_comparator_mixed_oracle.csv`. In the wide-weight case,
reference 2 has a log-density gap of approximately `-1e-100` relative to
reference 4; their final percentiles are correctly `.5` and `1`.

Both weighted fixtures retain their full ordering under common coordinate and
bandwidth rescaling by `2^±200`. The zero-inflated fixture also retains the
genuine entity-41-over-49 distinction and entity 41's percentile `8.5/96` at
all three unit scales.

## Scope of acceptance

The comparator proves ordering of its **represented scaled kernel masses**.
It does not make the underlying Float64 exponential evaluations exact, nor
promise to resolve every possible true Gaussian-density difference. The
half-exponential rational-square fallback extends represented range but does
not provide arbitrary-precision kernel evaluation. These limits are distinct
from the proven strict ordering of the represented keys.

Raw original-unit log densities can round to ties even when normalization keys
remain distinct. Documentation and audit oracles must use the explicitly
provided normalization scores to reproduce percentiles. This is a deliberate
normalization-contract change, not an assertion that the old raw-score CDF
contract remains unchanged.
