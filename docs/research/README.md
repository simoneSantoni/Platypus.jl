# Correctness and algorithm research log

Work started 2026-10-18 against commit `5ab596f` (v0.1.0).

## Mandate and procedure

The user requested implementation of the multidisciplinary review findings,
primary-source online research, comparison and adoption of suitable advanced
algorithms, and a record of the procedure. Correctness and internal consistency
take precedence over speed. Agreement with human perceptions is explicitly out
of scope. The panel comprises statistics, mathematics, Julia engineering,
computational social science, sociology, science communication, and Julia HPC.

1. Preserve the review's counterexamples as regression tests.
2. Define each estimator's target, reference-weight semantics, normalization,
   input requirements, and singular-case policy before optimizing it.
3. Search primary papers and official implementation documentation. Compare
   suitability, mathematical assumptions, numerical stability, cost and testable
   invariants; authorship or novelty alone is not an adoption criterion.
4. Implement in disjoint modules, cross-review, and compare numerical kernels
   with independent reference calculations.
5. Run the package suite, deterministic simulation validation, executable
   documentation and reproducible benchmarks. Record measured results, including
   failures and corrections, without presenting smoke tests as validation.

## Baseline evidence

- Existing suite: 109/109 checks passed on Julia 1.10.12 in the initial review.
- Confirmed failures: nonfinite inputs hidden by ranks; row-order-sensitive LOF
  ties; KDE weight-rescaling and zero-weight exclusion; insufficient reference
  support; mutable cache/metadata inconsistencies; unreliable tutorial recovery.
- Structural limitations: pairwise projection cannot identify all higher-order
  patterns; percentile scores are descriptive, not significance probabilities;
  centrality, local density contrast and copula density have different targets.

## Initial decisions

- Preserve endogenous reference collections and nonstructural labels.
- Preserve rank preprocessing as the default, but apply reference weights to
  preprocessing as well as estimation and final score normalization.
- Use exact neighbor search and original tie-inclusive LOF; approximate search
  must not silently change the statistic.
- Keep normalized summaries on a common explicit scale and expose raw method
  scores separately. Reject ambiguous legacy raw-summary requests with guidance.
- Adopt advanced methods only when their assumptions and implementation can be
  tested. In particular, marginal-tail methods such as ECOD/COPOD are not
  automatically substitutes for unusual feature-combination scoring.
- Planned release version is 0.2.0 because validation and normalization contracts
  change intentionally; no release or deployment is performed by this work.

## Primary-source screening started

- [Breunig et al., LOF (2000)](https://cs.fit.edu/~pkc/ml/related/breunig00sigmod.pdf)
  — exact kth-distance neighborhoods and local reachability density.
- [Segers, Sibuya and Tsukahara, empirical beta copula](https://arxiv.org/abs/1607.04430)
  — genuine smooth empirical copula under continuous untied ranks; do not assume
  arbitrary weighted or tied extensions preserve its margins.
- [Geenens, Charpentier and Paindaveine, probit copula estimation](https://arxiv.org/abs/1404.4414)
  — boundary-aware local-likelihood estimation; naive probit KDE is not the
  published improved estimator.
- [Li et al., COPOD](https://arxiv.org/abs/2009.09463) and
  [Li et al., ECOD](https://arxiv.org/abs/2201.00382) — computationally attractive
  tail scoring, requiring scrutiny of whether the target matches this package.

## Results and navigation

The current [numerical repair and validation record](numerical_repair.md)
documents implemented fixes, algorithm comparisons, new independent controls,
seven-specialist acceptance reports and a separate before/after audit. Earlier
review and benchmark records below describe their recorded source versions.
The final repaired audit passes all 14,637 checks across 84 datasets: 630 original
flags become zero, with the KDE normalization-contract change explicitly
documented. The final package suite passes 3,701 assertions and the independent
numerical comparator suite passes 976 assertions.

The [renewed panel investigation](panel_followup/README.md) traces the simulation
failures to numerical causes, adds independent counterexamples and isolated
interventions, and reviews the package's contracts and validation gaps. It is a
review, not a repair; its findings qualify the earlier closure record.

The subsequent [simulation consistency audit](simulation_consistency.md) found
additional unresolved numerical defects using 84 datasets and independent small
counterexamples. Its 630 violation flags qualify the earlier closure record;
the simulation audit preserves the findings without modifying estimators.

The subsequent [Julia 1.13.2 migration](julia_1_13_2.md) records the new default
runtime, dependency-environment repair and validation on that runtime.

The panel implemented the corrections and selected algorithm additions in the
working tree. The [closure and validation record](validation.md) maps review
concerns to changes, records integrated checks, and distinguishes fixed defects
from deliberate modeling limits. [CHANGELOG.md](../../CHANGELOG.md) describes
the intentional 0.2 API and score changes.

- [Statistics](statistics.md): weighting, reflected KDE, Gaussian/beta copulas,
  screened alternatives and numerical derivations.
- [Mathematics](mathematics.md): independent oracles, invariants, singular
  supports, normalization and higher-order limits.
- [Julia engineering](julia.md): data integrity, shared storage, snapshots and
  final cross-review overflow corrections.
- [Julia HPC](hpc.md) and [graph methods](graph.md): exact search, tie handling,
  sparse graph design, memory costs and reproducible measurements.
- [Computational social science](computational_social_science.md): primary-source
  method screening, deterministic recovery comparisons, negative controls and
  retained failures, including the separate graph comparison.
- [Science communication](communication.md): executed examples, API explanations,
  documentation procedure and migration.
- [Sociology](sociology.md): supplied-reference estimands, bibliography repairs,
  temporal and descriptive interpretation boundaries.

Procedure refinements during integration included deterministic weighted sums,
explicit Gaussian endpoint handling, representable extreme-value arithmetic,
and a strict documentation build that revealed an undeclared direct DataFrames
dependency in the docs environment and an omitted module docstring. Both docs
issues were corrected rather than disabling the checks. Human-perception
validation remained outside scope throughout.
