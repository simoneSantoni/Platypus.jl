# Julia engineering: numerical and structural repair validation

This note records the engineering portion of the numerical-resolution work on
18 October 2026. The original review and its outputs remain unchanged in
`panel_followup/julia.md` and `validation/consistency/panel/julia/`.

## Pair coverage and source provenance

`_validate_pairset` now requires both the combinatorial pair count and unique
ordered index tuples with valid `i < j`. Together these conditions imply that
every unordered feature pair is present exactly once. Merely checking vector
length previously allowed one projection to replace another and silently
changed the aggregation weights. Validation accepts arbitrary ordering of the
complete pair set; this is not an imposed canonical ordering.

Known transforms require source data. Conversely, source-backed pairs require
a known transform. Coordinate/source agreement is still rechecked at fit time.
The historical six-argument constructor retains unknown provenance and uniform
typicality fitting, but cannot establish the raw data needed for weighted
preprocessing or copula fitting. Malformed eight-field objects fail with an
argument error at analysis entry rather than reaching an unrelated method error.

The focused regression exercises duplicate and missing pairs across kNN, LOF,
KDE, depth, Gaussian copula and beta copula. Reversing all pairs must preserve
raw and normalized scores after reversing their columns, and preserve every
aggregate exactly. This passing control prevents a repair that accidentally
rejects legitimate pair permutations.

## Original reference mass and export identity

Weight normalization sums the exact rational representations of the supplied
Float64 weights and rounds each final fraction once. This avoids reduction-order
changes, overflowing finite sums and changes caused by appending zero masses.
Unrepresentable positive normalized masses are rejected. The normalized Float64
vector is a convenient export representation; it is not reused as an exact
replacement for the original weights. Repeated floating normalization need not
be idempotent, and no such requirement is imposed.

The estimation integration retains original source weights for subsequent
geometry and score percentiles and records them as `source_reference_weights`.
Metadata also carries `entity_ids` in corresponding order, binding each mass to
an entity after the exported score table is reordered or filtered. The regression
uses nonuniform weights and deliberately unsorted identifiers, checks an explicit
identifier join, and verifies metadata snapshots remain detached. Outcome
measurement provenance inherits this identity binding; graph metadata separately
records its own identifiers. Graph fits currently use uniform reference mass.

Requested and effective neighbor backends are distinct metadata fields. The
current conservative exact-distance implementation records a brute-force
effective backend even when the caller requested a KD-tree. This repairs an
ambiguity in recorded settings; the earlier missing-pair defect was the separate
issue that directly corrupted aggregation.

## Precision and concurrency review

The interim exact-moment implementation exposed two additional numerical hazards.
Converting a squared standardized deviation to Float64 before taking its square
root can overflow although the root is representable: reference coordinates
`[0, 1]` and a zero-mass query at `1e200` are a direct example. Also,
`BigFloat(x; precision=256)` sets the precision of that value, not of subsequent
arithmetic. On Julia 1.13.2, square rooting that value under a 16-bit ambient
precision produced a 16-bit result. Ordinary conversion from a BigInt rational
to Float64 also used ambient BigFloat arithmetic, and multiplying a rational
sign by Float64 could promote the result back to BigFloat.

The shared numerical helpers address these issues using integer/rational
arithmetic and final Float64 rounding, without changing global MPFR state.
The contracts regression varies ambient precision only in a serial test and
requires identical rank/z-score preprocessing, depth scores, normalized scores
and exported weights. Density and neighbor tests cover their numerical paths
separately. The huge-query regression requires finite standardized coordinates
and checks their known values.

Per-pair threads read shared coordinate vectors and write separate result
columns and parameter entries. Exact rational temporaries are local to each
calculation. These properties avoid a shared arithmetic-state requirement;
they do not authorize a caller to mutate inputs concurrently. Public field
arrays remain mutable by design, and fit-time validation is not a lock. The
`PairSet` documentation now explicitly states this concurrency boundary.

## Julia indexing compatibility

The detached-snapshot `EntitySet` indexing implementation previously iterated
Boolean masks as scalar Boolean indices. It now resolves collection indices
against the entity index range before constructing snapshots, restoring ordinary
Boolean-mask selection, repeated integer indices and bounds checks. Selected
entities remain detached from the authoritative matrix.

## Reproduction and acceptance

Run the focused test from the repository root:

```sh
JULIA_DEPOT_PATH=/tmp/platypus-review-julia-depot:/home/simon/.julia \
julia --startup-file=no --project=. -e \
'using Test, Platypus, DataFrames; include("test/numerical_contracts.jl")'
```

Focused validation passed **91/91 checks** on Julia 1.13.2 after integration of
the shared rational conversion and square-root helpers (18 October 2026).

Acceptance requires the structural rejection and passing permutation controls,
exact source-weight normalization controls, metadata identity joins and snapshot
checks, ordinary mask indexing, and precision/finite-range regressions to pass.
Integration results are recorded separately by the coordinating reviewer.

## Primary Julia references

- [Composite types](https://docs.julialang.org/en/v1/manual/types/#Composite-Types):
  immutable objects can contain mutable arrays, so field immutability does not
  remove validation or concurrent-mutation concerns.
- [Inner constructors](https://docs.julialang.org/en/v1/manual/constructors/#Inner-Constructor-Methods):
  constructor invariants and later mutation validation are separate concerns.
- [Logical indexing](https://docs.julialang.org/en/v1/manual/arrays/#Logical-indexing):
  Boolean masks select the positions of their true elements.
- [BigFloat precision](https://docs.julialang.org/en/v1/base/numbers/#Base.MPFR.setprecision):
  ambient precision changes affect operations and are not thread-safe.
- [Data races](https://docs.julialang.org/en/v1/manual/multi-threading/#Communication-and-data-races-between-threads):
  callers must synchronize conflicting accesses; collection validation cannot
  make unsynchronized concurrent edits safe.
