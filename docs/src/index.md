# Platypus.jl

**Descriptive atypicality relative to a chosen collection.**

Platypus.jl scores how unusual an entity's measured features are within the
collection you supply. Pairwise estimators score every two-feature projection
and summarize each entity's scores. A separate graph estimator describes
connectivity in the full selected feature space.

Scores depend on the features, reference population, and method. They are not
anomaly probabilities or validated measures of human judgments of novelty.
The observed collection is a reference chosen by the analyst. Calling it
endogenous does not resolve its substantive boundaries.

## What is implemented

The 0.2 development code implements data loading, exploratory diagnostics,
LOF/kNN/KDE/depth scoring, Gaussian and empirical beta copulas, sparse graph
scoring, descriptive outcome regression, and plotting recipes. A working-tree
version is not evidence that a matching release tag has been published.

The final 84-dataset numerical audit passed 14,637 checks with zero violations
and zero unexpected comparison errors. The package suite passed 3,701 assertions,
and the independent mathematical comparator checks passed all 976 assertions.
Repairs preserve selected algebraic and reflection identities, use exact
neighbor-boundary decisions, and replace spectral depth truncation with exact
algebraic rank. Real-coordinate KDE exposes a separate normalization ordering;
its weighted mid-CDF contract is explicit. The
[Estimation guide](guide/estimation.md) explains these contracts, computational
costs, and remaining finite-precision limits. Validation covers the tested
scenarios rather than every possible input.

Start with [Getting Started](getting_started.md), then choose a method using the
[Estimation guide](guide/estimation.md). Plotting functions return recipe objects;
[Visualization](api/visualization.md) explains rendering.

## The pairwise workflow

For ``N`` entities and ``p`` selected features, let
``P = \binom{p}{2}`` denote the number of feature pairs.

1. Construct an [`EntitySet`](@ref) with a common feature schema.
2. Use [`build_pairs`](@ref) to form its projections.
3. Score with [`fit_typicality`](@ref) or [`fit_copula`](@ref).
4. Inspect individual pairs, summaries, and saved configuration.

The default transformation replaces each feature with its empirical mid-rank:

```math
u_{ej} = \frac{\operatorname{tiedrank}(x_{\bullet j})_e - 1/2}{N}.
```

Untied continuous features receive a uniform grid of ranks. Ties retain masses;
constant features remain constant. Ranking removes marginal units and changes
geometry. It is a transformation, not a randomized null model, and does not
guarantee that every score responds exclusively to dependence violations.

Each pair has physical or method-specific `raw_scores`, explicit
`normalization_scores`, and normalized `pair_scores`. Larger ordering values mean
more typical. Let ``r_{ep}`` denote the value in `normalization_scores`.
Normalization uses its reference-weighted mid-CDF:

```math
s_{ep} =
\frac{\sum_i w_i\,\mathbf{1}(r_{ip}<r_{ep})+
\tfrac12\sum_i w_i\,\mathbf{1}(r_{ip}=r_{ep})}{\sum_i w_i}.
```

Equal normalization values share the midpoint of their reference mass; if all
are equal, all percentiles are 0.5. Normally `normalization_scores` equals
`raw_scores`. Real-coordinate Gaussian KDE instead supplies ordinal labels
ordered by dimensionless kernel mass, before physical-unit log normalization can
erase a small distinction. Raw log-density ties can therefore have different
percentiles. These labels have no density units or comparable magnitude across
pairs or fits.

Exact rational calculations and canonical symmetry rules preserve selected
identities; no blanket tolerance groups nearby scores. Represented exponential
terms and final arrays remain floating, with finite precision and underflow
limits. The percentile is descriptive, not a probability of normality or anomaly.

## Reading the summaries

| Summary | Definition | Interpretation |
| --- | --- | --- |
| Conventionality | Median over pairs | Higher means more typical overall |
| `novelty` | Lower quantile, default ``Q_{0.1}`` | **Lower** means more unusual |
| Dispersion | IQR over pairs | Variation across pairs, not a standard error |
| Atypicality | One minus conventionality | Higher means more atypical overall |

The lower quantile is not the minimum. Neither a median nor a lower quantile
guarantees detection of one unusual pair among many. Pairwise projections cannot
identify every higher-order relationship.

The design draws inspiration from Uzzi et al.'s median and lower-tail summaries
of randomized journal co-citation scores. Platypus measures different continuous
feature configurations; its rank normalization is not Uzzi et al.'s null-model
standardization. Dispersion is an additional package summary.
See [References](references.md) and the
[Research workflow](guide/research_workflow.md) for interpretation and validation.

## Reference and representation

Labels are retained for inspection; they do not select peers. Supported methods
accept [`ReferenceWeights`](@ref) for graded membership. With those methods,
zero-weight rows are queries and do not contribute to transformations, reference
fits, or normalization.

LOF measures local relative density. KDE measures density, Mahalanobis depth
measures global centrality, and graph scores describe full-feature connectivity.
These quantities can produce different orderings.

Save [`analysis_metadata`](@ref) with exported scores. It binds reference weights
to ordered entity IDs and records configuration and caller-supplied provenance.
Pairwise metadata retains the original represented weights as well as their
normalized values. Neighbor fits record requested and effective backends; both
requests currently execute a filtered exact scan with quadratic candidate work.
Retain the dataset, source revision, code, Julia environment, and substantive
selection decisions separately.

## Reading order

- [Getting Started](getting_started.md): an executed deterministic example.
- [Entities and Data](guide/entities.md): feature selection and integrity.
- [Building Feature Pairs](guide/pairs.md): transforms and snapshots.
- [Exploratory Analysis](guide/eda.md): diagnostics and their limits.
- [Estimation](guide/estimation.md): methods, assumptions, weights, and cost.
- [Research workflow](guide/research_workflow.md): reporting and validation.
- [Estimation API](api/estimation.md) and [Visualization](api/visualization.md):
  signatures and output types.

## Citation

If you use Platypus.jl in your research, please cite the software. See
[Citing Platypus.jl](api/citation.md) for the citation, BibTeX entry, and guidance
on recording the version and methods used.
