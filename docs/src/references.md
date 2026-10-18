# References

These sources establish distinct measurement ideas and numerical methods.
Platypus combines selected ideas under its own [measurement contract](guide/measurement_scope.md);
it does not reproduce every cited study or claim that the methods measure the
same quantity. Algorithm selection follows the intended numerical target,
assumptions and testable properties, rather than the prominence of a source.

## Combination atypicality and cultural products

- **Uzzi, B., Mukherjee, S., Stringer, M., & Jones, B.** (2013).
  [Atypical combinations and scientific impact](https://doi.org/10.1126/science.1240474).
  *Science*, 342(6157), 468–472.
  Compares journal-pair co-citation counts in Web of Science with randomized
  citation networks. The two paper-level summaries are the median and tenth
  percentile of null-standardized scores. Platypus borrows this aggregation
  idea; its IQR is an additional diagnostic, and its percentiles have no
  randomized-null interpretation.
  [Full paper](https://gwern.net/doc/psychology/novelty/2013-uzzi.pdf).

- **Askin, N., & Mauskapf, M.** (2017).
  [What makes popular culture popular? Product features and optimal differentiation in music](https://doi.org/10.1177/0003122417728662).
  *American Sociological Review*, 82(5), 910–944.
  The yearly specification uses genre-weighted cosine similarity to songs on
  the preceding 52 weeks' charts. This motivates explicit feature and reference
  choices, but differs from Platypus's pairwise density scores and shared
  reference weights. Descriptive OLS in `fit_outcome` is not a replication of
  the paper's outcome models.
  [Full paper](https://gwern.net/doc/culture/2017-askin.pdf).

## Categories and reference-relative descriptions

- **Rosch, E., & Mervis, C. B.** (1975).
  [Family resemblances: Studies in the internal structure of categories](https://doi.org/10.1016/0010-0285(75)90024-9).
  *Cognitive Psychology*, 7(4), 573–605.
  Studies relations between category prototypicality and shared attributes.
  Historical conceptual context; it does not derive a local-density estimator
  or equate typicality with membership weights.

- **Hsu, G., Hannan, M. T., & Koçak, Ö.** (2009).
  [Multiple category memberships in markets: An integrative theory and two empirical tests](https://doi.org/10.1177/000312240907400108).
  *American Sociological Review*, 74(1), 150–169.
  Relates category spanning to audience-side and producer-side processes.
  Context for distinguishing category membership from feature-space position;
  Platypus's descriptive labels and shared reference weights do not implement
  that theory's category-specific mechanisms.

## Density, local outliers and centrality

- **Breunig, M. M., Kriegel, H.-P., Ng, R. T., & Sander, J.** (2000).
  [LOF: Identifying density-based local outliers](https://cs.fit.edu/~pkc/ml/related/breunig00sigmod.pdf).
  *Proceedings of ACM SIGMOD 2000*, 93–104.
  Source of kth-distance neighborhoods including boundary ties, reachability
  density and local density ratios. Platypus records explicit conventions for
  zero-distance cases and exposes reciprocal LOF as typicality.

- **Silverman, B. W.** (1986).
  [Density Estimation for Statistics and Data Analysis](https://ned.ipac.caltech.edu/level5/March02/Silverman/Silver_contents.html).
  Chapman and Hall. Background on kernel and nearest-neighbor estimation and
  bounded domains. The package's dimensional effective-size bandwidth and
  compact-kernel scale adjustment are documented heuristics, not a claim of
  optimal bandwidth for every supplied distribution.

- **Zuo, Y., & Serfling, R.** (2000).
  [General notions of statistical depth function](https://doi.org/10.1214/aos/1016218226).
  *The Annals of Statistics*, 28(2), 461–482.
  Framework for depth functions, including Mahalanobis depth. Platypus uses
  weighted mean/covariance centrality with an explicit singular-support policy;
  it implements neither Tukey halfspace depth nor simplicial depth.

- **Muia, M. N., Atutey, O., & Hasan, M.** (2025).
  [Kernel smoothing for bounded copula densities](https://arxiv.org/abs/2502.05470).
  Preprint. A source screened for the mirror-reflection construction. Platypus
  adopts bounded-support reflection with its own compact product kernel and
  bandwidth rule, not the complete estimator-selection procedure in the paper.
  Its rank-space KDE is not guaranteed to have exactly uniform fitted margins.

## Copulas

- **Nelsen, R. B.** (2006).
  [An Introduction to Copulas](https://doi.org/10.1007/0-387-28678-0), second edition.
  Springer. Background on copulas, marginal distributions and dependence.
  The package's Gaussian family estimates normal-score correlation; it does
  not claim maximum-likelihood estimation or universal dependence recovery.

- **Segers, J., Sibuya, M., & Tsukahara, H.** (2016 preprint).
  [The empirical beta copula](https://arxiv.org/abs/1607.04430).
  Defines the rank-based beta mixture. Platypus requires untied margins and
  equal positive reference weights to retain its copula construction, and
  explicitly evaluates the complete in-sample mixture.

## Full-feature graph connectivity

- **Zelnik-Manor, L., & Perona, P.** (2004).
  [Self-tuning spectral clustering](https://papers.neurips.cc/paper/2619-self-tuning-spectral-clustering.pdf).
  *Advances in Neural Information Processing Systems*, 17.
  Source of locally scaled affinities. Platypus adopts the affinity construction
  on a symmetrized neighbor graph; it does not implement the paper's spectral
  clustering procedure.

- **Lovász, L.** (1993).
  [Random walks on graphs: A survey](https://www.cs.cmu.edu/afs/cs/academic/class/15859n-f16/RelatedWork/random-walks-on-graphs.pdf).
  *Combinatorics, Paul Erdős is Eighty*, 2, 1–46.
  Source for degree-proportional stationary mass of an undirected walk.
  Disconnected graphs require a stated selection among stationary distributions;
  stationarity alone does not imply convergence from every initial distribution.
  Platypus computes its selected mass directly, without reinforcement or ants.

Additional algorithm screening and adoption decisions are recorded in the
repository's `docs/research/` notes. The executable study-design examples are
in [Reproducible research workflow](guide/research_workflow.md).
