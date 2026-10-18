# Independent final KDE consistency diagnosis

Examined `first_full/datasets/zero_inflated_2026.csv`, pair `(f1,f2)`, 96
uniformly weighted entities. `probe.jl` uses Julia Base only; it does not import
the package estimator. Input decimal CSV values are first parsed as Float64,
then promoted exactly to 512-bit BigFloat. The oracle evaluates Gaussian kernel
sums at the represented bandwidth values and, separately, at high-precision
rule-of-thumb bandwidths. Both give the same strict ordering.

The culprit is entities **41 and 49**, not an actual density tie:

| Entity | f1 | f2 | Oracle log density |
|---|---:|---:|---:|
| 41 | 1.3304618062110518 | 2.235492214628469 | -3.6523572063482661256717223770905… |
| 49 | 1.3265525472336013 | 2.545792979854194 | -3.6523572063482665178669550292384… |

The strict difference `logdensity(41)-logdensity(49)` is
`3.9219523265214789216620617785e-16`. It remains positive using fully
high-precision bandwidth calculation (`3.9219523265214423917650226096e-16`).
There is therefore no basis for merging these observations with a numerical
tie tolerance. The correct midpercentile for entity 41 is `8.5/96`; merging the
two gives `8/96`, producing the observed `1/192` discrepancy.

Changing both units by `2^±200` adds a common log-density constant of magnitude
about 277. Float64 spacing there is too coarse to preserve the distinction.
The resulting raw displayed log densities can legitimately coincide after
rounding even when the underlying density ordering differs.

**Factoring out log bandwidth alone is insufficient.** The independent probe
computes the scale-free log mean kernel score and obtains the same Float64
value `-5.142518787131657` for entities 41 and 49 at all three unit scales.
This makes the invariant check pass by consistently losing a genuine
distinction; it does not agree with the mathematical oracle.

Computing the dimensionless kernel mean directly in the linear domain,

```text
S_i = sum_{j != i} exp(-0.5 * ((u_i-u_j)^2/hu^2 + (v_i-v_j)^2/hv^2)) / 95,
```

retains the distinction: `S_41=0.00584295400571043`,
`S_49=0.005842954005710428`. All 96 midpercentiles then agree with the
512-bit oracle at all three unit scales. This calculation forms the sum
directly; exponentiating the already-rounded log score cannot restore the
lost information.

For general inputs, a Float64 linear representation can underflow for distant
queries. A general repair therefore needs a representation or comparison
policy that also retains very small positive densities: higher-precision
normalization keys, or a scaled mantissa/exponent comparison with refinement
for unresolved comparisons. Approximate ties should not silently replace
mathematical ties. Preserve the original-unit raw logs for display and document
that any new normalization keys, rather than rounded raw logs, determine ranks.

`oracle.txt` records full 512-bit values and all three scale comparisons.
The Float64 formula in this independent probe is an algebraic reconstruction,
not a claim of bit-for-bit identity with every production evaluation order.
The high-precision ordering conclusion does not depend on that detail.
