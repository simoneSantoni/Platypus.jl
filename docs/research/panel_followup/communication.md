# Communication panel: claims, evidence, and release status

Reviewed 2026-10-18 against the working-tree 0.2.0 source and the renewed panel
evidence on Julia 1.13.2. This is a review, not a repair. Only this report was
written; public documentation, source, tests, and prior evidence were unchanged.
Human-perception validation remains outside scope.

## Assessment

The public explanation of descriptive scores, reference selection, method
differences, and novelty direction remains useful. Its numerical assurances are
now too strong. The central README and user guide do not surface the confirmed
unresolved defects, while the earlier closure table still reads as a current
resolution record. A passing suite and an implemented API must not be presented
as numerical acceptance.

The renewed [panel synthesis](README.md) handles the main evidence distinctions
well: 630 dependent flags are not 630 bugs; the 525 small-raw flags are not 525
proved mathematical ties; and the 91-case intervention is an isolated prototype,
not a production correction. Keep those qualifications beside every shortened
summary of these numbers.

## Prioritized passages and proposed replacements

### C1. High: the primary status description omits known unresolved correctness findings

**Passages:** `README.md:22–27`; `docs/src/index.md:15–24`.
These state that the methods are implemented and that 0.2.0 is unreleased, but
give no numerical-status qualification. Neither statement is false alone;
together with the advertised guarantees below, they leave readers with an
unsupported impression of completed correctness work.

**Evidence chain:** `docs/research/simulation_consistency.md:18–22` records the
630 flags; `docs/research/panel_followup/README.md:39–59` links independent
counterexamples to source mechanisms. The unchanged-source record is
[environment.txt](../../../validation/consistency/panel/integration/environment.txt).

**Proposed addition immediately after the implementation status:**

> Numerical consistency remains unresolved in this development version. An
> expanded audit found changes in some scores and neighborhoods under operations
> that should preserve the stated measurement. Independent small examples confirm
> defects in arithmetic identities, distance decisions, and singular-support
> handling. The existing test suite passes but does not cover these failures.
> See the renewed panel review before treating the scores as release-validated.

“Descriptive” does not excuse a failure to compute the declared descriptor.
The caveat belongs in the entry points, not only several links into research notes.

### C2. High: “exact” conflates the retrieval algorithm with the numerical predicate

**Passages:** `README.md:99` and `:113–114`;
`docs/src/guide/estimation.md:32–37` and `:133–137`;
`CHANGELOG.md:19–20` and `:42`.

**Evidence chain:** `src/neighbors.jl:46–72` uses rounded `hypot` values in the
final membership comparison. Both backends share that predicate. The
[HPC review](hpc.md), lines 37–106, distinguishes intended rank-lattice geometry,
an exact integer graph fixture, and exact rational distances on frozen Float64
coordinates. It shows both missing true boundary neighbors and included
genuinely unequal neighbors. No candidate-pruning failure was observed; a
universal pruning guarantee was not established.

**Proposed replacement for the search paragraph:**

> LOF/kNN use KD-tree retrieval or a streaming brute-force alternative, without
> approximate-neighbor retrieval or a dense distance matrix. Both currently use
> floating-distance boundary predicates with confirmed membership errors on
> selected tied and close-distance fixtures. Their agreement does not certify
> exact mathematical neighborhoods. The intended LOF definition includes every
> kth-distance tie; that contract is not yet fully met.

Retain the engineering claim about storage. Do not rename the implementation
“approximate search”: retrieval policy and floating arithmetic are different
issues. After repair, any stronger exactness claim needs its input/transform
semantics and predicate/fallback evidence stated explicitly.

### C3. High: intended identity and support rules read as achieved guarantees

**Passages:** `README.md:90–93` (“all-equal pair scores”);
`docs/src/index.md:56` and `:81–84`;
`docs/src/guide/estimation.md:73–77` and `:107–110`.

**Evidence chain:** `src/estimation.jl:80–101` correctly ranks equal *computed*
values. The depth triangle, weighted constant KDE, and beta-polynomial oracle
show that equal mathematical scores can arrive unequal at that function.
The [statistics review](statistics.md), lines 99–135, establishes that zero-weight
augmentation can change the existing reference's computed Gaussian scores.
The [mathematics review](mathematics.md), section 2, identifies exact comparison
for axis-constant support versus a tolerance tube for oblique singular support.

**Proposed qualifications:**

> Exactly equal computed raw scores map to 0.5. This conditional normalization
> rule does not guarantee that the current estimators preserve mathematical ties;
> independently confirmed cases produce unequal percentiles from equal targets.

> Zero-weight rows are intended to contribute no reference mass. Known differences
> between preprocessing paths can nevertheless change existing scores when such
> rows are appended.

> The declared depth convention assigns zero outside the reference affine span.
> Current singular-support branches apply inconsistent numerical tolerances, so
> that convention is not yet satisfied across coordinate representations.

Do not replace these with a vague “floating-point noise” disclaimer: one support
counterexample changes raw depth from 1 to 0, and some distance defects change
graph edges. Conversely, do not call the empirical-CDF formula itself incorrect.

### C4. High: the historical closure record still asserts current resolution

**Passages:** `docs/research/validation.md:11–29`, especially rows at lines
17, 19, 22, 23 and 26; `docs/research/README.md:80–84`;
`docs/research/communication.md:77`.

The closure record currently has a runtime supersession notice at lines 7–9 but
no comparable notice that subsequent numerical findings reopen its conclusions.
“Resolution and evidence,” “exact kth-distance neighborhoods,” “off-support
depth zero,” and “all ties map to 0.5” should not remain unqualified completion
claims. The research index's earlier warning at lines 67–75 helps, but does not
travel with a direct link to the closure document.

**Proposed prominent historical notice:**

> Historical implementation record. These changes and test runs occurred, but
> the subsequent simulation audit and renewed panel review reopened numerical
> correctness findings. This table records attempted resolutions and their
> original evidence, not current release acceptance. See the follow-up for
> unresolved cases and required validation.

Preserve old test counts and runtime details; do not rewrite history or discard
passing evidence. Replace the research index's “distinguishes fixed defects”
with “records the earlier implementation changes and tests; subsequent findings
qualify which defects can be considered resolved.” Replace “completed
implementation” in the earlier communication note with “implementation at that
stage.” The latest 911/911 pass establishes only the existing assertions.

### C5. Medium: exported provenance is useful but insufficiently bound to observations

**Passages:** `README.md:50–52` and `:73`;
`docs/src/guide/estimation.md:120–123`;
`docs/src/guide/research_workflow.md:42–47`;
`docs/src/api/estimation.md:68–70`.

**Evidence chain:** `src/estimation.jl:188–194` stores ordered weights in
configuration and IDs separately in the result; exported metadata has no
ID-to-weight binding. The [sociology review](sociology.md), lines 108–131,
demonstrates the gap. The
[Julia probe summary](../../../validation/consistency/panel/julia/summary.txt)
independently shows four incorrect positional matches after sorting a nonuniform
example's score table. It also shows requested KD-tree metadata despite actual
brute fallback: compare `src/neighbors.jl:35–43` with
`src/estimation.jl:270–272`.

**Proposed interim export advice:**

> Save ordered entity IDs alongside the metadata's weight vector before sorting
> or joining tables; the current metadata does not bind weights to IDs. Preserve
> input and source fingerprints and the resolved environment separately. Pairwise
> backend metadata records the requested backend and may not reveal a fallback.

Detached metadata is a valid ownership claim; it is not a completeness or
execution-attestation claim. Acceptance should require ID-bound weights,
requested/effective backend distinctions, and recorded numerical rank/support
policies once those contracts are decided.

## Concise user-facing explanation

> The package still has numerical consistency defects despite all 911 existing
> tests passing. The expanded audit recorded 630 overlapping flags; independent
> examples show that rounding and inconsistent numerical decisions can change
> scores and neighborhoods under equivalent representations. A small prototype
> removed 91 feature-order failures in isolation, which identifies one cause but
> does not repair the package. The review recommends targeted numerical fixes and
> independent end-to-end regression checks before numerical release acceptance.

Avoid describing all 525 small-raw cases as ties or all 630 flags as distinct
bugs. The original report's maximum normalized-output changes are not necessarily
headline atypicality changes. Where downstream effects are discussed, use the
aligned replay in the
[computational social science review](computational_social_science.md), lines
61–105, which separates headline changes, order reversals, and illustrative
threshold switches. It does not establish real-world prevalence or outcome
coefficient bias.

## Release acceptance language and procedure

A defensible future claim would name the repaired production revision, numerical
contract, tested domain, regression/oracle cases, unchanged comparison rules,
remaining limitations, and preserved before/after audit results. Require separate
checks of transformations, parameters, raw identities, geometric membership,
support decisions, normalized scores, and consequential ordering. Retain the
genuine-close-score and unequal-distance controls against blanket rounding.
A successful script exit or a larger test count alone is insufficient.

For the present tree, say “implemented, reviewed, unresolved numerical defects.”
Do not say “all issues fixed,” “numerically consistent,” or “exact mathematical
ties preserved.” The prototype does not change the unresolved status.

Procedure for this review:

1. Read the follow-up synthesis and mathematics, statistics, HPC, sociology and
   computational-social-science notes; inspect the original audit and historical
   closure/communication records.
2. Compare their claims with current README, public guides/API pages, changelog,
   normalization, neighbor predicates, pair validation, and metadata construction.
3. Independently sum the saved triage CSV: 630 total, 525 within the raw diagnostic
   threshold, 105 exceeding it. Check the feature-exchange CSV: 91 cases, with
   zero nonzero prototype raw or normalized discrepancies. These are artifact
   checks, not reruns or independent proofs of the numerical oracles.
4. Read the integration environment and Julia probe summaries, retaining the
   unchanged source fingerprint and distinction between numerical errors,
   metadata gaps, and successful snapshot controls.
5. Propose wording changes here only. No public pages or package files were edited.

No additional external authority was needed for this communication audit: its
claims concern local source, preserved observations, and the specialist evidence
chain linked above. Primary numerical sources are recorded in those specialist
reports. This panelist did not rerun their numerical probes or assign a proven
root cause to every audit flag.
