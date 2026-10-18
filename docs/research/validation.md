# Implementation closure and validation record

Revision date: 2026-10-18. Baseline: commit `5ab596f`, version 0.1.0.
Working package version: **0.2.0, unreleased**. This work edits the checkout;
it does not publish a release, deploy documentation, or commit changes.

The subsequent [Julia 1.13.2 default migration](julia_1_13_2.md) supersedes the
development-runtime choice below. Earlier run versions and results are retained
as historical evidence.

## Finding-to-resolution record

| Review concern | Resolution and evidence |
|---|---|
| Nonfinite features hidden by ranks; invalid identity/schema metadata | Validate source data and metadata before transformation and again at analysis boundaries. `test/data_integrity.jl`. |
| Entity vectors diverging from cached matrix; mutation leaking into results | One canonical entity matrix; detached entity/accessor snapshots; source-aware pair validation; copied result metadata. `test/data_integrity.jl` and integration tests. |
| Weights changing under rescaling or leaking through excluded rows | Stable relative membership normalization; weighted preprocessing, fitting and mid-CDF normalization; explicit unsupported-method rejection. Density and mathematical-property tests. |
| Empty leave-one-out reference or hidden epsilon repairs | Require sufficient positive reference support; explicit constant/singular support conventions and exact zero-density log scores. Density and numerical-limit tests. |
| Row-order-dependent neighbor ties; arbitrary duplicate treatment | Exact kth-distance neighborhoods include boundary ties; original LOF reachability with documented duplicate extension; independent brute-force comparisons. `test/neighbors.jl`. |
| Dense quadratic distance matrices; repeated pair coordinates | KD-tree search with streaming brute-force fallback; shared per-feature vectors; optional independent-pair threads; sparse graph with explicit edge cap. Benchmarks and tree/brute/thread parity tests. |
| Unbounded kernels applied to bounded rank coordinates | Exactly normalized reflected compact kernels on unit support; Gaussian kernels on real coordinates; bandwidth heuristic and fixed-reference kernel LOO documented. Quadrature and analytic density oracles. |
| Scale-sensitive covariance regularization | Dimensionless spectral Mahalanobis calculation, explicit affine-support tolerance, off-support depth zero. Independent covariance, affine and unit-scaling checks. |
| Mixed raw-summary scales; unclear novelty direction | Mandatory weighted mid-percentiles for summaries, separate raw scores; lower novelty means more unusual; all ties map to 0.5. Independent summary calculations. |
| Unjustified redundancy weights | Removed the suggested weight column; retain diagnostic absolute correlations and undefined constant-feature correlations. |
| Exported estimator and plotting placeholders | Implement Gaussian/beta copulas, full-feature graph, descriptive QR outcome fitting, and all nine recipe constructors. Family-specific assumptions and distinct graph result type documented. |
| Missing reproducibility information and unreliable examples | Detached analysis metadata; explicit feature selection/cohort workflow; deterministic executed tutorial and scenario experiments; recorded environments. |
| Overstated methodological lineage and interpretation | Corrected citations, Uzzi's two summaries versus added IQR, estimand distinctions, shared-reference limits and descriptive rather than inferential meaning. Measurement-scope guide and sociology notes. |
| Broken or permissive documentation workflow | Portable docs project, strict exported-docstring checks, executed examples, explicit deployment opt-in and unreleased migration notes. |
| Overflow found in final independent cross-review | Scale outcomes before QR; stable standardized Gaussian differences; validate bandwidth conversion to Float64. `test/numerical_limits.jl` includes 28 public-API regressions. |

## Algorithm selection

The seven specialists screened primary papers and official implementations.
Selection depended on matching the package's numerical target and being able to
verify the implementation, not on author reputation or recency alone.

| Candidate | Outcome |
|---|---|
| Exact LOF and kNN | Repair and optimize existing estimators without approximate neighbors or arbitrary tie truncation. |
| Reflected rank-space KDE | Adopt a normalized compact-kernel construction; distinguish it from the complete estimator/bandwidth optimizer in the screened 2025 reflection paper. |
| Gaussian copula | Adopt explicit normal-score correlation estimation; reject singularity unless regularization is requested. |
| Empirical beta copula | Adopt the complete equal-weight mixture with untied ranks; reject unsupported weighted/tied extensions. |
| Self-tuning affinity graph | Adopt as a separate full-feature connectivity estimator, with explicit disconnected-graph stationary-mass selection. |
| Probit local-likelihood copula density | Defer rather than substitute a plain transformed KDE for the published local-likelihood method. |
| ECOD, COPOD, Isolation Forest, LoOP, INFLO | Screen and document target/assumption differences; no unsupported implementation or superiority claim. These external algorithms were not empirically benchmarked here. |

See [statistics](statistics.md), [HPC](hpc.md), [graph](graph.md), and
[computational social science](computational_social_science.md) for primary
sources, search procedures, derivations and adoption/defer decisions.

## Comparative evidence

The deterministic pairwise study retains 180 rows: six methods, five scenarios,
three seeds and two summaries. The graph supplement retains 15 headline recovery
runs and separate stability diagnostics. All are generated-data experiments,
with 160 observations and 12 designated anomalies per collection, not a broad
real-world benchmark. See `validation/results/` for the CSVs and environments.

Pairwise methods recover the constructed multi-pair dependence violation well,
but recovery weakens after irrelevant features are added. Graph connectivity
performs better in the selected higher-order and irrelevant-feature cases. All
methods fail to recover the designated bridge cluster adequately. Changing the
feature schema materially changes rankings, even when row permutations and exact
backends agree. No method is declared universally best; numerical agreement is
not evidence of agreement with human judgments. The latter is explicitly outside
this phase.

## Reproduction and environment

Normal project commands are:

```sh
julia --threads=2 --project -e 'using Pkg; Pkg.test()'
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'
julia --project validation/experiments.jl
julia --project validation/graph_experiments.jl
julia --project benchmark/neighbors.jl
julia --project benchmark/graph.jl
```

Local validation uses a writable temporary depot over the installed depot,
`JULIA_DEPOT_PATH=/tmp/platypus-review-julia-depot:/home/simon/.julia`, and
`--startup-file=no --compiled-modules=no`. The initial Julia 1.13.0 environment
failed during dependency precompilation; numerical and documentation validation
therefore uses Julia 1.10.12 and the supported minimum Julia 1.9.4. A subsequent
normal compiled-module package load and fit passed on Julia 1.10.12. This is not
a claim that every Julia version has been verified locally. CI retains the
minimum/current Julia matrix and now uses two threads.

The Julia 1.9.4 Linux binary was downloaded from the official archive and checked
against SHA256 `07d20c4c2518833e2265ca0acee15b355463361aa4efdab858dad826cf94325c`.
The hardened host rejected that old runtime's executable-stack flag on
`libopenlibm.so`; the flag was cleared only in the temporary extracted runtime.
No package source or dependency implementation was modified for that adjustment.
The minimum-version project resolved compatible dependencies independently
(including DataFrames 1.7.1 and NearestNeighbors 0.4.24); Julia 1.10 used
DataFrames 1.8.2 and NearestNeighbors 0.4.29.

Focused suites, executable examples, comparative experiments and benchmark
measurements are recorded in the specialist notes. Timings are workload
observations on a shared host; allocation counts are cumulative allocations,
not measured peak memory.

## Executed integrated checks

`Pkg.test` passed **898/898 checks on both Julia 1.9.4 and Julia 1.10.12**, with
two Julia threads. The testsets contain 148 package/API, 60 data-integrity,
150 neighbor, 133 density/copula, 137 mathematical-property, 219 graph,
15 workflow, 28 extreme-numerical-value and 8 serial/thread agreement checks.
The final integration review added the 28 numerical-limit regressions to the
earlier 870-check suite.

After the visualization fixes, the complete current suite passed **911/911
checks on Julia 1.10.12 with normal compiled modules and two threads**. This adds
13 label, overlay, user-override and color-key regressions to the 898 checks
above. The final 13 plotting regressions were not separately run on Julia 1.9.4.

Normal compiled-module loading and a small default fit also passed on Julia
1.10.12. The 19 named tutorial/guide examples passed independently. The strict
Documenter build passed without warnings, checking the actual isolated docs
environment, executed examples, cross-references and exported docstrings.
Its initial failures revealed missing direct DataFrames visibility and an
unincluded module docstring, both fixed in the documentation configuration.

All nine plotting helper names, plus the graph distribution overload, rendered
to SVG using Plots 1.41.7 and GR 0.73.27 in an isolated temporary environment
resolved from cached packages. Visual inspection then exposed axis-label/title
overrides in the outcome identity line and an implicit heatmap colorbar setting;
both were corrected and covered by the 13 plotting regressions. The preserved
`validation/render_plots.jl` then passed **23/23 actual-backend assertions**,
rendered all ten SVG figures, and checked processed subplot labels and color
keys. Final outcome and heatmap images were visually inspected: axis labels,
title and numeric colorbar are present. Generated figures stayed in a temporary
directory rather than becoming repository artifacts.

## Deliberate limits

- Pairwise densities do not identify the full joint distribution. This is not a
  blanket impossibility theorem for every aggregation of pair functions.
- No calibrated anomaly probability, randomization null, confidence interval,
  causal estimator or frozen historical-reference prediction API is claimed.
- Nonuniform weights remain unsupported for neighbor/graph methods and empirical
  beta copulas. Tied margins are unsupported by the empirical beta family.
- Bandwidth, k, feature schema and reference collection remain explicit modeling
  choices; adding sophisticated algorithms does not eliminate them.
- Exact KDE/beta evaluation requires quadratic kernel work, with additional
  sorting in deterministic reductions (up to O(N² log N) work). Score storage is
  proportional to entities times feature pairs, and exact neighbor ties can be
  expensive. Sparse graph storage has an enforced edge limit.
- Direct public array fields remain mutable; copied accessors and validation
  protect normal workflows but do not turn Julia arrays into immutable objects.
