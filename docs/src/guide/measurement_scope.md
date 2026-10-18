# Measurement scope

Platypus measures statistical atypicality relative to a supplied collection,
feature schema, transformation and scoring rule. These choices define the
quantity being measured. The collection supplies an endogenous reference;
the algorithm does not discover a unique population boundary or reconstruct
observations missing from the collection.

## Match the algorithm to the numerical target

There is no single algorithm for every meaning of atypicality. Even with the
same data, local isolation, global centrality, dependence density and graph
connectivity can give different orderings. Agreement between these orderings
is not an implementation requirement.

| Method | Raw numerical target | Interpretation boundary |
|---|---|---|
| `fit_typicality(...; method=:knn)` | Negative exact kth-neighbor distance in each pair | Neighborhood crowding at the selected k; not a calibrated density or probability |
| `method=:lof` | Reciprocal tie-inclusive Local Outlier Factor in each pair | Density contrast with neighbors; a small coherent cluster need not be anomalous |
| `method=:kde` | Leave-one-out product log density in each pair | Reflected compact kernels on rank support; Gaussian kernels on unbounded coordinates; depends on bandwidth |
| `method=:depth` | Weighted global Mahalanobis depth in each pair | Covariance-based centrality with exact algebraic rank and strict affine support; not local density |
| `fit_copula(...; family=:gaussian)` | Fitted log Gaussian-copula density in each pair | Normal-score correlation model; singularity requires explicit regularization |
| `family=:beta` | Log empirical beta-copula density in each pair | Complete in-sample mixture; untied margins and equal positive weights required |
| `fit_graph_typicality` | Selected degree-proportional stationary mass on a full-feature similarity graph | Connectivity under locally scaled affinities; low density need not mean low connectivity |

Rank preprocessing reduces the influence of marginal scales. Ties remain tied;
finite empirical margins need not be continuously uniform. A rank transform
does not itself specify a randomization null, and a bounded KDE on rank
coordinates need not be a copula density. Distances, boundary effects and the
chosen scoring rule still matter. The Gaussian and empirical beta copula
families have explicit dependence models; neither fits every possible joint
distribution.

The pairwise methods inspect all two-feature projections. Identical bivariate
distributions can accompany different higher-order constraints, leaving no
population pair-density contrast. This does not prove that every aggregation of
pair scores must be powerless, but pairwise density does not identify the full
joint distribution. The graph retains all selected features jointly, but its neighborhood and affinity
choices still determine which configurations receive low connectivity.
Adding a feature also changes the collection of pairs; related or irrelevant
features can change summary rankings. Method and feature sensitivity are
numerical questions to examine for the intended target.

## The summary contract

For each pair, `normalization_scores` are mapped to reference-weighted
mid-percentiles in `pair_scores`. Larger values mean more typical. Equal
normalization values share the midpoint of their reference mass; if all are
equal, all percentiles are 0.5.

Normally `normalization_scores` equals the method-specific `raw_scores`.
Gaussian KDE on real coordinates (`:none` or `:zscore`) is the explicit
exception: physical log densities remain in `raw_scores`, while
`normalization_scores` holds ordinal labels of dimensionless kernel masses.
Those labels preserve computed order and equality before a physical-unit log
adjustment can erase a small contrast. Equal rounded raw logs can therefore
have different percentiles. The labels are not density measurements; their
magnitudes have no interpretation across pairs or fits. Rank-space KDE uses
raw values for normalization.

Export physical/method-specific values with `scoretable(result; raw=true)` and
ordering values with `scoretable(result; normalization=true)`. The options are
mutually exclusive. Metadata records `normalization_basis=:raw_score` or
`:scale_free_kernel_order`. This is an explicit normalization-contract change,
not a claim that the original raw-score mid-CDF contract was always satisfied.

Real-coordinate KDE compares direct and shifted kernel masses on a common
basis, preserving original reference-weight ratios outside the exponential.
Underflowed mixtures are not assigned to a separate low-score group.
The numerical comparison certifies the represented mass expressions, not exact
original Gaussian densities; terms can still round or underflow after shifting.
All 976 independent mathematical comparator checks passed, including all eight
mixed-underflow oracle rows. The final 84-dataset audit passed 14,637 checks with
zero violations and zero unexpected comparison errors. These results support
implementation correctness in the tested scenarios, not construct validity. See the
[estimation guide](estimation.md). No blanket tolerance merges nearby scores.

For pairwise results:

- `conventionality` is the median of an entity's pair percentiles.
- `novelty` is their selected quantile, normally 0.1: **smaller means more
  unusual**. It describes the lower tail, not the single minimum or historical
  first appearance.
- `dispersion` is their interquartile range. It measures spread across pairs,
  not uncertainty about the estimate.
- `atypicality` is `1 - conventionality`: larger means more atypical under the
  selected rule. One unusual pair can leave this median-based score unchanged.

The median and lower-tail design is inspired by Uzzi et al.; IQR is an
additional package diagnostic. Their randomized co-citation baseline is not
implemented here. Percentiles are neither anomaly probabilities nor p-values,
and their common scale does not make different methods interchangeable.
See the [source distinctions](../references.md).

Graph results have their own `typicality`, `atypicality`, degree and stationary
mass. They do not fabricate pairwise conventionality or novelty. For a
disconnected graph the reported stationary distribution is an explicit
degree-based selection; it is not unique. Read the recorded component and
isolation conventions when interpreting that result.

## Reference membership and labels

`ReferenceWeights` defines relative reference mass. KDE, depth and the Gaussian
copula use supported weights throughout marginal transformation, reference
fitting and percentile normalization. Zero-weight rows are query observations
without reference mass; positive common rescaling of weights leaves the target
unchanged. This is a membership-weight contract, not a survey-design estimator,
an inverse-probability correction, or a count of independent replicates.
The empirical beta copula requires equal positive weights; kNN, LOF and the
graph reject nonuniform weights rather than silently ignoring them.

Entity identity and the unit represented by a row are inputs to the analysis.
Splitting one row into several records is not guaranteed to preserve scores,
even if their weights sum to the original mass: self-exclusion, neighborhoods
and effective-sample-size calculations can change. Resolve entities before
fitting and retain the ordered entity IDs and source reference weights in the
exported metadata.

Each fit uses one shared reference-weight vector. Labels are descriptive
metadata and never automatically define reference groups. A string label does
not encode graded multiple memberships, and the API does not infer separate
reference systems for different labels. If distinct collections are required,
construct and record them explicitly. Their percentiles describe different
reference distributions and are not a fixed absolute scale.

## Time and outcomes

Construct dated cohorts before transformation and fitting. Appending later
positive-weight observations can change historical ranks, neighborhoods and
scores even if the earlier entities' features do not change. A chronological
prefix limits the available observations, but a fit including the focal entity
is not automatically a strictly prior comparison. Leave-one-out KDE omits the
focal kernel; it does not refit every transformation and bandwidth without that
entity. There is no general fitted-reference/predict API.

`fit_outcome` performs descriptive least squares on selected pairwise summaries.
Its optional interaction uses conventionality times `1 - novelty`, so greater
tail unusualness enters with the intended direction. The fit supplies neither
causal identification nor uncertainty estimates and does not establish a
historical reference design. It is not a reproduction of the cited scientific
impact or music-success models.

Use the [research workflow](research_workflow.md) for explicit feature selection,
cohort construction, provenance, and deterministic recovery experiments. Save
`analysis_metadata(result)` with exported tables so these measurement choices
remain identifiable.
