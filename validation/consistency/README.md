# Simulation consistency audit

This directory contains 84 reproducible datasets across 28 families, the audit
scripts, and the observed results on Julia 1.13.2. The objective is to identify
internal inconsistencies, not to validate agreement with human judgments or
select a method by its ability to recover planted observations.

**The repaired replay passes all 14,637 checks with zero flags**, compared with
630 original flags. It has zero unexpected errors, 1,464 successful baseline
fits and 48 documented model rejections. All eight small oracle replay checks
also pass. The final source passes 3,701 package assertions and 976 independent
mathematical comparator assertions; the sociology review adds 33 controls.

The before/after comparison verifies identical input bytes, check keys, seeds
and tolerances. **The KDE normalization contract changes explicitly:** the
weighted-mid-CDF oracle now uses exposed `normalization_scores`, which retain
kernel ordering when rounded physical log densities collide. Raw log-density
unit checks and representation-invariance checks remain in place. See the
[repair rationale and limits](../../docs/research/numerical_repair.md).

**The original audit found numerical defects.** That run recorded 14,007 passing
checks and 630 violation flags; these are correlated checks, not 630 independent
bugs. Small examples independently confirm incorrect score distinctions,
unstable LOF neighborhood membership, and loss of depth accuracy after adding
a large exactly representable offset. See the
[findings and procedure](../../docs/research/simulation_consistency.md).

The subsequent [panel investigation](../../docs/research/panel_followup/README.md)
adds root-cause evidence, independent controls, and package-review findings.
Its probes and outputs are preserved separately under `panel/`; each specialist
report gives its reproduction command. The subsequent
[numerical repair record](../../docs/research/numerical_repair.md) documents the
implemented fixes and their validation. Repaired outputs are saved separately
under `repaired_results/`; acceptance probes are under `repair_panel/`.

## Run

From the repository root, using its Julia project:

```sh
julia --threads=2 --startup-file=no --project=. validation/consistency/run.jl /tmp/platypus-consistency
julia --threads=2 --startup-file=no --project=. validation/consistency/reproduce.jl /tmp/platypus-consistency
julia --threads=2 --startup-file=no --project=. validation/consistency/verify_artifacts.jl /tmp/platypus-consistency
python3 validation/consistency/plot_results.py /tmp/platypus-consistency
python3 validation/consistency/compare_repairs.py validation/consistency/results /tmp/platypus-consistency
```

The optional output directory defaults to `validation/consistency/results`.
Use a different directory to preserve the original observations. Julia scripts
use dependencies already declared by the package. Plotting additionally needs
Python, NumPy, and Matplotlib. The numerical audit does not require Python.

The scripts **record violations and finish normally** so that every scenario is
examined. Successful process exit means the investigation completed; it does
not mean the package passed. `reproduce.jl` asserts the mathematical premises of
its independent oracles, then writes the package's discrepancies to CSV.
These scripts are deliberately separate from the package's normal CI suite.

## Coverage

Families include independent uniform/Gaussian data, heavy tails and skew,
positive and negative dependence, a parabola, ring and moons, unequal clusters
and bridges, disrupted dependence, ordinal and zero-inflated data, duplicates,
constant features/clouds, exact and near collinearity, a Cartesian grid, parity,
mixed measurement units, large offsets, tiny/huge scales, twelve features, and
two/three-row datasets. Seeds are 17, 101, and 2026. Nominal sizes are 32, 63,
and 96 respectively; structured and very small families override these sizes.

Methods: LOF, kNN distance, KDE, Mahalanobis depth, Gaussian copula, Gaussian
copula with regularization 0.05, empirical beta copula, and graph typicality.
Pair and graph methods use rank, raw, and z-score coordinates; copulas use rank
coordinates. Nonuniform reference checks target methods that support them.

Checks cover row/feature permutations, labels/IDs, monotone recoding, reflection,
unit scaling, reference-weight scaling, exclusion of zero-weight queries,
exact search backends, serial/threaded execution, bounds, an independent
weighted-mid-CDF calculation, constant clouds, and graph stationarity.
Changes to `k` are recorded as sensitivity, not treated as inconsistencies.

## Saved outputs

| File | Contents |
|---|---|
| `datasets.jl` | Seeded generators; no package estimators used in generation |
| `run.jl` | Full metamorphic/property audit |
| `reproduce.jl` | Small fixtures with exact rational and high-precision oracles |
| `plot_results.py` | Figures derived from saved CSVs |
| `verify_artifacts.jl` | Input round-trip, table-count, identity and source-fingerprint checks |
| `compare_repairs.py` | Identical-input/check/tolerance assertions and before/after flag comparison |
| `repaired_results/comparison.csv` | Before/after flags, resolved flags and new flags by method |
| `repaired_results/package_tests.log` | Complete repaired-source package test results |
| `repaired_results/audit.log` | Final full replay execution log |
| `repaired_results/provenance.txt` | Final input/table integrity and frozen source verification |
| `repaired_results/figures/before_after.*` | Before/after flags by method; PNG/PDF |
| `results/datasets.csv` | Shapes, seeds, duplicates, constants, numerical ranges |
| `results/datasets/*.csv` | All 84 input matrices with IDs, weights and diagnostic labels |
| `results/baselines.csv` | Successful fits and documented model rejections |
| `results/checks.csv` | Every attempted check, tolerance metrics and status |
| `results/counterexamples.csv` | A representative row/pair for every flagged comparison |
| `results/scores.csv` | Baseline per-entity atypicality scores |
| `results/sensitivity.csv` | Descriptive changes at `k=1` and `k=n-1` |
| `results/summary.csv` | Counts by method, transform, property and status |
| `results/minimal_*.csv` | Small inputs, expected/observed scores, neighbors and checks |
| `results/environment.txt` | Runtime, tolerances and package-source fingerprint |
| `results/provenance.txt` | Platform, dependencies and script fingerprints |
| `results/figures/*` | Dataset atlas, consistency heatmap, small counterexamples; PNG/PDF |

Each input CSV contains `entity_id`, `planted`, `label`, `reference_weight`, then
`f1`, `f2`, … . `planted` is only a generator diagnostic; the two alternating
labels are nonstructural metadata. Neither supplies a correctness target.
In `minimal_examples.csv`, `NaN` in `score_expected` means only the raw-score
oracle was evaluated. It is not a nonfinite package result.

Counterexample rows maximize the pair-percentile discrepancy; a row therefore
need not exhibit the maximum raw or summary discrepancy reported in
`checks.csv`. Graph entries use stationary mass as the raw value. Equal signed
infinities are allowed for zero log density. No package numerical implementation
was changed during the original audit; its failures remain as historical evidence.

![Small counterexamples](results/figures/minimal_counterexamples.png)
