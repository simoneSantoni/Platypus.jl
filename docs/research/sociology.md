# Measurement scope and scholarly attribution

Review date: 2026-10-18. This implementation pass addresses numerical correctness,
internal consistency and faithful attribution. Agreement with human perceptions
is outside the requested scope and is not an acceptance gate.

## Procedure

The initial review compared source, examples and tests against the package's
claims about reference populations, novelty, categories and its cited studies.
The implementation pass inspected `estimation.jl`, `density.jl`, `graph.jl`,
`outcome.jl`, and the reference-weight contract. It separates properties of the
implemented statistic from contextual ideas borrowed from the literature.

Targeted online searches used the titles and authors of Uzzi et al. (2013),
Askin and Mauskapf (2017), Hsu et al. (2009), and Rosch and Mervis (1975), plus
the empirical beta copula, statistical depth and self-tuning graph papers.
Full primary papers were used where accessible; publisher records verified
metadata. Repository mirrors of primary papers remain primary texts. This is
a focused claim audit, not an exhaustive systematic literature review.

## Evidence and adopted corrections

| Source inspected | Evidence relevant to the package | Adopted distinction |
|---|---|---|
| [Uzzi et al., pp. 468–470](https://gwern.net/doc/psychology/novelty/2013-uzzi.pdf) | Journal-pair counts standardized against randomized citation networks; median and tenth percentile are the two summaries | Describe aggregation as inspired by this study; IQR is a package addition; rank percentiles provide no null calibration |
| [Askin and Mauskapf, pp. 9–11](https://gwern.net/doc/culture/2017-askin.pdf) | Full-feature cosine similarity with a preceding 52-week reference and pair-specific genre weights | A shared endogenous collection with one reference-weight vector is a different specification; OLS is not a replication of their outcome models |
| [Hsu et al., publisher record and abstract](https://journals.sagepub.com/doi/abs/10.1177/000312240907400108) | Category spanning incorporates producer- and audience-side processes; published pages are 150–169 | Correct pagination; distinguish descriptive labels and reference mass from category-specific mechanisms |
| [Rosch and Mervis, publisher abstract](https://doi.org/10.1016/0010-0285(75)90024-9) | Shared attributes relate to category prototypicality | Retain as historical conceptual context without claiming a derivation of density or reference weights |
| [Breunig et al., definition of k-distance neighborhood](https://cs.fit.edu/~pkc/ml/related/breunig00sigmod.pdf) | LOF compares local reachability densities with boundary ties included | Distinguish local contrast from global rarity and honor the original neighborhood definition |
| [Segers, Sibuya and Tsukahara](https://arxiv.org/abs/1607.04430) | Rank-based beta-mixture copula construction | Name the actual restricted estimator; no arbitrary weighted/tied extension is implied |
| [Zelnik-Manor and Perona](https://papers.neurips.cc/paper/2619-self-tuning-spectral-clustering.pdf) and [Lovász](https://www.cs.cmu.edu/afs/cs/academic/class/15859n-f16/RelatedWork/random-walks-on-graphs.pdf) | Local scaling defines affinities; undirected walks admit degree-proportional stationary mass | Call the result full-feature connectivity, specify disconnected-graph selection, and remove reinforcement/ant-colony claims |

The bibliography's incomplete “Kim, H., Lee, J., et al.” entry was removed
because no title, year or identifiable work supported its attribution. Broad
background citations that did not identify implemented methods were pruned;
the remaining annotations explain what is and is not adopted.

## Operational decisions checked against implementation

1. The supplied collection, features, transform and score rule define the
   estimand. An endogenous reference removes the need for an external benchmark;
   it does not identify a unique population boundary.
2. Pairwise percentile summaries retain median conventionality, lower-quantile
   typicality named `novelty`, IQR dispersion and complementary atypicality.
   Lower `novelty` means greater unusualness. Historical newness, minimum score
   and inferential significance are different quantities.
3. Reference weights mean relative membership mass. Supported estimators apply
   them to transformation, fitting and score normalization. They do not implement
   survey sampling, effective independent replicate counts or missing-population
   reconstruction. Unsupported weighting is rejected explicitly.
4. Pairwise distance, local contrast, smooth density, global depth and copula
   density remain distinct targets. The graph adds a joint-feature representation
   with its own connectivity target. No method universally covers all notions
   of atypicality; method adoption depends on assumptions, analytical identities,
   invariants and computational feasibility, not a source's prominence.
5. Descriptive OLS uses the supplied summaries and outcome alignment. An optional
   conventionality-by-tail-atypicality interaction has explicit sign semantics.
   It does not implement either cited study's full design or uncertainty model.
6. Cohorts must be defined before preprocessing. A dated prefix is an endogenous
   cohort, not automatically a strictly prior fitted reference. KDE's omitted
   focal kernel is not complete pipeline leave-one-out fitting. The API has no
   general reusable fitted-reference/predict separation.

The resulting public contract is in
[`measurement_scope.md`](../src/guide/measurement_scope.md), with the corrected
[`references.md`](../src/references.md). The computational-social-science
panel owns the executable cohort/feature safeguards and diagnostic experiments
in [its notes](computational_social_science.md) and the
[research workflow](../src/guide/research_workflow.md); these are not duplicated
here. This pass changes documentation only. Numerical test results belong to
the individual implementation records and the integrated validation report.
