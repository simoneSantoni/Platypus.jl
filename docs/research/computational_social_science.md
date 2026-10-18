# Computational social science: design, comparison and validation

## Review response

The review found that the original examples did not establish recovery, numeric
column discovery could admit outcomes, and historical comparisons and provenance
were underspecified. The revised research workflow uses explicit `feature_cols`,
dated collections constructed before ranking, detached `analysis_metadata`, and
an explicit interpretation of lower-tail scores. `test/workflow.jl` exercises
these contracts. Statistical measurement and numerical reliability are the scope;
human-perception validation is excluded.

The endogenous-reference design is preserved. Chronological prefix construction
prevents later observations from entering a fit; it is not an external-reference
prediction interface. Historical novelty and causal impact remain distinct
research questions, not capabilities supplied by a descriptive score or OLS fit.

## Primary-source search and decisions

Searches covered ECOD, COPOD, Isolation Forest, anomaly benchmark protocols and
the temporal design of Askin–Mauskapf. Sources were inspected online during the
implementation review; algorithm selection below is our inference about the
package's target, not a performance claim from those papers.

| Candidate | Relevant primary evidence | Decision for this package |
|---|---|---|
| ECOD | Estimates coordinate-wise empirical distributions and aggregates tail scores. [Li et al.](https://arxiv.org/abs/2201.00382) | Do not substitute it for pairwise dependence scoring: rare marginal positions and unusual dependence combinations are distinct targets. No claim of numerical superiority or inferiority; no ECOD implementation benchmark was run. |
| COPOD | Scores multivariate extremeness using empirical-copula-inspired tail probabilities. [Li et al.](https://arxiv.org/abs/2009.09463) | Tail probability is not local copula density or local density contrast. Keep it as an external comparator for a future matched-target benchmark, rather than treating the shared word “copula” as algorithmic equivalence. No COPOD implementation benchmark was run. |
| Isolation Forest | Uses random partition trees and isolation path length to detect anomalous instances. [Liu, Ting and Zhou](https://cs.nju.edu.cn/zhouzh/zhouzh.files/publication/icdm08b.pdf) | Potential full-dimensional comparator, but changes the score's meaning and loses the native pairwise explanation. No automatic adoption and no empirical performance claim from this review. |
| Exact LOF/kNN, reflected KDE, global depth, Gaussian/beta copulas | Implemented methods expose different structural targets and assumptions. See specialist mathematical/statistical notes. | Compare these actual implementations on the same generated collections; do not rank algorithms solely by speed or a single recovery example. |

Emmott et al. vary anomaly difficulty, frequency, clusteredness and feature
relevance, and report both ROC AUC and average precision. Their findings motivate
our limited fixed-design checks and explicit negative controls, rather than a
claim to replicate their much broader benchmark. [Primary paper](https://arxiv.org/abs/1503.01158)

Askin–Mauskapf's yearly measurement uses the prior 52 weeks, illustrating why a
pooled retrospective collection is a different temporal design. [Primary paper](https://gwern.net/doc/culture/2017-askin.pdf)

## Executable comparison

Run `julia --project validation/experiments.jl`. It writes
`validation/results/recovery.csv`, `summary.csv`, `stability.csv` and
`environment.txt`. Seeds 101, 202 and 303 define independent synthetic
collections of 160 entities with 12 designated anomalies. No settings are tuned
using those labels. All six methods use rank coordinates and their defaults;
Gaussian and beta copulas are evaluated only on untied continuous data.

- **Independent features:** designated labels have no feature signal. The
  expected AUROC is 0.5; low scores in a fitted collection are not discoveries.
- **Dependence violations:** four features share a noisy latent component;
  reversing two features in designated rows violates four of six pairs, so the
  headline median has sufficient affected pairs to respond.
- **Unequal clusters:** two clusters differ in size and dispersion, with bridge
  observations designated anomalous. Global centrality, density contrast and
  density estimation can disagree legitimately. Rank preprocessing also changes
  raw geometry, so this is not a benchmark of raw-space density estimation.
- **Irrelevant features:** the dependence case receives four independent noise
  features. This tests dilution by the all-pairs summary without changing the
  planted labels.
- **Higher-order structure:** `U`, `V` and `(U+V) mod 1` have independent
  bivariate margins but a joint constraint. Replacing the third coordinate for
  designated rows removes that constraint while preserving their population
  pairwise distributions. The population bivariate densities have no contrast;
  finite-sample scoring and global centrality can still behave differently.

Recovery is reported for headline atypicality and `1 - novelty`, separately.
AUROC gives tied scores half credit. Average precision groups equal-score
thresholds. Precision at the known anomaly count averages over ties at the
boundary. That count is known only in the simulation; it is not a production
threshold. Prevalence (0.075) is recorded for interpreting precision. Three
seeds give an illustrative range, not a confidence interval.

Stability checks independently permute rows, compare exact tree/brute neighbor
backends and correlate rankings before/after adding irrelevant features. The
benchmark also records allocation and elapsed time after one method warm-up;
these timings are indicative, and controlled scalability belongs to the HPC
benchmark. No held-out or human judgment data are used. Every generated case is
retained, including poor recovery. Additional datasets and sampling designs are
needed before choosing a method for a substantive application.

## Executed results

Executed with Julia 1.10.12, Platypus 0.2.0, one thread, seeds 101/202/303.
The artifacts contain 180 recovery rows and 18 stability rows. Mean headline
AUROC over the three seeds is shown below; per-seed values, min/max, average
precision and the separate tail summary are retained in the CSV files.

| Scenario | LOF | kNN | KDE | Depth | Gaussian copula | Beta copula |
|---|---:|---:|---:|---:|---:|---:|
| Independent features | 0.502 | 0.532 | 0.545 | 0.520 | 0.463 | 0.495 |
| Dependence violations | 0.963 | 0.980 | 0.969 | 0.866 | 0.992 | 0.825 |
| Unequal clusters / bridge observations | 0.343 | 0.154 | 0.504 | 0.090 | 0.320 | 0.096 |
| Dependence plus irrelevant features | 0.561 | 0.620 | 0.636 | 0.574 | 0.737 | 0.575 |
| Higher-order constraint | 0.502 | 0.557 | 0.588 | 0.461 | 0.532 | 0.539 |

Recovery is strong for the constructed multi-pair dependence violation. It
deteriorates when independent features dilute the median. The bridge observations
are not recovered: all six methods have mean precision at 12 equal to zero in
that scenario. They form their own small cluster, and rank preprocessing changes
the original gap geometry. The experiment therefore does **not** support a broad
claim that these defaults detect arbitrary anomalies in heterogeneous data.
Gaussian copula leads this particular dependence construction, which is favorable
to its assumptions; it does not win across all cases or establish a universal
replacement for LOF. Three replications are insufficient to establish an exact
chance-level claim for the higher-order case.

All 18 row-permutation comparisons had maximum headline-score difference 0 and
Spearman correlation 1. All six exact neighbor backend comparisons had maximum
raw-score difference 0. These support implementation consistency on the tested
continuous datasets. Correlations after adding irrelevant features ranged from
0.211 to 0.795, demonstrating material schema sensitivity despite numerical
consistency. `test/workflow.jl` passed 15/15 checks. Separate metric reference
checks verified perfect ordering, reversed ordering and all-tied scores.

Artifacts: [recovery](../../validation/results/recovery.csv),
[summary](../../validation/results/summary.csv),
[stability](../../validation/results/stability.csv),
[environment](../../validation/results/environment.txt).

## Full-feature graph supplement

`julia --project validation/graph_experiments.jl` reuses exactly the same five
generators, seeds and recovery metrics. It leaves the original 180 recovery rows
unchanged and writes a separate 15-row `graph_recovery.csv` and its summary.
The graph's headline atypicality measures **full-feature connectivity** under
self-tuning affinities; it is not a seventh pairwise density method. No graph
novelty or conventionality measure is fabricated. Defaults use rank coordinates
and `k=13` on these 160-row collections; no graph settings were tuned to labels.

Executed on Julia 1.10.12 with one thread:

| Scenario | Mean graph AUROC | Seed range | Mean average precision | Mean precision at 12 |
|---|---:|---:|---:|---:|
| Independent features | 0.556 | 0.526–0.606 | 0.105 | 0.056 |
| Dependence violations | 0.940 | 0.856–0.985 | 0.841 | 0.750 |
| Unequal clusters / bridge observations | 0.035 | 0.021–0.051 | 0.043 | 0.000 |
| Dependence plus irrelevant features | 0.882 | 0.867–0.891 | 0.623 | 0.528 |
| Higher-order constraint | 0.793 | 0.727–0.896 | 0.515 | 0.417 |

The full-feature graph recovers the planted higher-order violations better in
these three simulations than the six pairwise summaries (mean AUROC 0.461–0.588).
It also retains more recovery after the addition of irrelevant features in this
construction. However, it ranks the designated bridge cluster as highly typical
and performs worse there than the pairwise alternatives. Local scale adaptation
can make a small dense group well connected, even when that group was designated
anomalous by the generator. These results support using the graph as a distinct
connectivity target, not as a universal improvement or an automatic repair of
feature selection. The independent-feature control remains a negative control;
three seeds do not establish calibrated false-positive rates.

The graph supplement also saves row-permutation and irrelevant-feature ranking
sensitivity in `graph_stability.csv`. Edge and component counts are retained in
the recovery artifact so disconnected-graph behavior can be inspected instead of
implicitly assuming a unique stationary distribution.
Row permutations had zero score differences and correlation 1 in all three
checks. Adding noise features changed the full ranking substantially (Spearman
0.113–0.161), despite retaining high recovery of the designated anomalies.
The unequal-cluster case at seed 303 had two graph components; the other 14
graphs were connected. Recorded times are indicative: the four-feature warm-up
does not remove dimension-specific compilation on the first three/eight-feature
fits, so these rows must not be used for a controlled speed comparison.

Artifacts: [graph recovery](../../validation/results/graph_recovery.csv),
[graph summary](../../validation/results/graph_summary.csv),
[graph stability](../../validation/results/graph_stability.csv),
[graph environment](../../validation/results/graph_environment.txt).

## Remaining design limits

Reference selection, missingness and feature engineering remain analyst choices.
Metadata records those choices, but does not validate their truth or remove
selection bias. Paired measurements share variables and fitted references;
independent-row inferential assumptions require separate justification. The
package provides descriptive outcome fits without claiming causal or calibrated
inferential conclusions. These are documented limits, not silently deferred bug
fixes.
