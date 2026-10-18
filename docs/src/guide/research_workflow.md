# Reproducible research workflow

Platypus measures relative statistical structure in a supplied collection. Its
scores describe that collection and chosen features. They are not probabilities
of anomaly, hypothesis tests, historical first appearances, or causal effects.
The workflow below concerns numerical and statistical validation; it does not
attempt to measure human perceptions.

## Declare the collection and features

Record the data source, version, inclusion rules, observation period and missing
data decisions. Charting songs, indexed papers and surviving firms are selected
populations: local scoring cannot recover entities absent from the source.
Reference membership weights, recency weights and sampling weights encode
different substantive choices. State what the chosen weights mean.

Always select measurement features explicitly in outcome studies. Numeric
discovery includes numeric years, outcomes and controls unless excluded by an
explicit `feature_cols` list.

```julia
using Platypus, DataFrames

df = DataFrame(id=1:8, year=[2020,2020,2020,2021,2021,2021,2022,2022],
    energy=[.1,.2,.4,.3,.6,.7,.8,.5],
    acousticness=[.8,.7,.5,.6,.3,.2,.1,.4],
    streams=[100,120,200,180,300,350,400,250])
measurement_features = [:energy, :acousticness]
cohort = df[df.year .<= 2021, :]
# The package rejects missing and nonfinite features. Here the decision is
# explicit complete-case selection; record and inspect the exclusions.
complete = dropmissing(cohort, measurement_features)
es = load_entities(complete; id_col=:id, feature_cols=measurement_features)
result = fit_typicality(build_pairs(es); method=:kde,
    provenance=(dataset="illustrative-v1", cohort="year <= 2021", cutoff=2021,
                preprocessing="complete cases on energy/acousticness",
                excluded_ids=setdiff(cohort.id, complete.id)))
scores = summarytable(result)
metadata = analysis_metadata(result)
```

Save `metadata` with the exported score table and the script, input checksum and
resolved dependency environment. Metadata includes the transformation, normalized
reference weights, effective sample size, method parameters, version and the
supplied provenance. It records your declarations; it does not verify the source
data or the accuracy of a declared cutoff. Join by `entity_id`, verify the join is
one-to-one, and check outcome alignment after filtering or sorting.

Dropping incomplete observations changes the reference collection. Imputation
also changes its geometry; document the method and fit preprocessing within the
permitted time window. Compare reasonable preprocessing choices when missingness
is appreciable. There is no built-in missing-data model.

## Keep time explicit

Feature ranks, neighborhood structure and final percentiles are computed from
the supplied collection. Pooling years changes the historical comparison by
including later entities. Recompute separate dated cohorts or chronological
prefixes when the research question requires information available by a cutoff.
Construct the cohort **before** `load_entities` and `build_pairs` so future values
cannot influence marginal ranks.

A prefix ending at a focal entity includes that entity in the endogenous
collection; depth and copula fitting can include it, and transformations use the
declared reference weights. This is not automatically a strictly prior reference
or a model for scoring previously unseen observations. Platypus deliberately has
no external-reference API. A strict historical prediction design needs additional
methodological work; `fit_outcome` does not establish such a design from scores.
Askin and Mauskapf's yearly measure explicitly compares against the preceding
52 weeks, a distinct temporal specification. [Primary paper](https://gwern.net/doc/culture/2017-askin.pdf)

Percentiles in independently fitted cohorts are relative to different reference
distributions. A score change need not mean an entity's features changed. Retain
cohort IDs and avoid treating these values as a fixed absolute measurement scale.

## Interpret summaries and outcome fits

`atypicality` is one minus median pairwise typicality. `novelty` is a lower-tail
typicality score: **smaller means more unusual**. It is not a minimum, and one
unusual pair can leave the median almost unchanged. `dispersion` is variation
across feature-pair scores, not a standard error.

LOF measures density contrast with neighbors; kNN and KDE measure neighborhood
crowding or density; Mahalanobis depth measures global centrality. Copula fits
impose their own dependence assumptions. Choosing between them changes the
measurement target, not merely its computational implementation.

`fit_outcome` is a descriptive least-squares fit. It does not provide causal
identification, automatic control selection, historical validity or calibrated
uncertainty. Pair rows share features, and entity scores share the fitted
reference collection. Do not treat `scoretable` rows as independent observations.
For inferential work outside this package, choose resampling and clustering from
the research design (for example artist, author or firm and time), and refit the
measurement pipeline within each valid resample. Resampling only fixed score
rows omits uncertainty from estimating the reference collection.

## Check recovery and stability

The executable `validation/experiments.jl` compares all four typicality methods
and both copula families on fixed seeded collections. It reports AUROC,
threshold-based average precision, tie-aware precision at the planted anomaly
count, row-permutation stability and sensitivity to added irrelevant features.
The labels are generated from known numerical structure. Independent-feature
and higher-order-only negative controls are retained even when methods fail.

```sh
julia --project validation/experiments.jl
```

The separate graph supplement uses those same generated collections and reports
only graph atypicality, a full-feature connectivity target. It does not rerun
the pairwise comparisons or invent a graph novelty summary:

```sh
julia --project validation/graph_experiments.jl
```

This is a small diagnostic experiment, not a universal method ranking or a test
of predictive generalization. Different methods share a percentile output scale
but target different structures. Bivariate density targets do not identify a
joint constraint when every bivariate margin is identical to the comparison population.
Inspect individual pairs and compare feature schemas, neighborhood sizes,
bandwidths and reference definitions. A useful validation report explains both
successes and deliberate failure cases.

The design follows the methodological lesson of varying anomaly frequency,
cluster structure and feature relevance rather than reporting one convenient
example; the supplied experiment covers a small subset of those dimensions.
[Emmott et al.](https://arxiv.org/abs/1503.01158)
