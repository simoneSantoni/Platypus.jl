# References

Platypus.jl sits at the intersection of three literatures: the
*atypicality-and-impact* tradition that gives the pairwise design, the
*outlier- and depth-based density estimation* tradition that gives the
scoring rules, and the *categorisation* tradition that frames typicality
as a graded, reference-relative property. Citations on the API pages
refer to entries in this list.

## Atypical combinations and impact

- **Uzzi, B., Mukherjee, S., Stringer, M., & Jones, B.** (2013).
  Atypical combinations and scientific impact. *Science*, 342(6157),
  468–472. The originating paper. Measures a publication's atypicality
  from the *pairs* of journals in its reference list — how often each
  pair is co-cited across all of science, relative to a randomised null —
  and summarises the resulting distribution by its median
  (*conventionality*) and 10th percentile (*novelty*). The finding that
  the highest-impact papers are *both* highly conventional *and*
  tail-novel is the design Platypus.jl generalises from discrete
  co-citation counts to continuous features.

- **Askin, N., & Mauskapf, M.** (2017). What makes popular culture
  popular? Product features and optimal differentiation in music.
  *American Sociological Review*, 82(5), 910–944. Applies a typicality
  measure to songs, using audio features, and finds a non-monotone
  ("optimal differentiation") relationship with chart success. The
  closest prior application of the present idea to a music market, and
  the natural comparison for the planned
  [`fit_outcome`](@ref Platypus.fit_outcome) estimator.

- **Kim, H., Lee, J., et al.** treatments of recombinant novelty in
  patents and science more broadly motivate measuring novelty as a
  property of *combinations* rather than of single attributes; the Uzzi
  et al. design is the canonical statement and is taken as the reference
  point here.

## Categorisation and typicality

- **Rosch, E., & Mervis, C. B.** (1975). Family resemblances: Studies in
  the internal structure of categories. *Cognitive Psychology*, 7(4),
  573–605. The foundational statement that category membership is graded
  and that typicality is a function of an item's *configuration of
  features* relative to other members — the cognitive-science root of
  treating atypicality as endogenous and feature-based.

- **Hsu, G., Hannan, M. T., & Koçak, Ö.** (2009). Multiple category
  memberships in markets: An integrative theory and two empirical tests.
  *American Sociological Review*, 74(1), 150–166. Develops typicality and
  category spanning for markets with *fuzzy, overlapping* category
  boundaries — the organisational-theory counterpart of the "no clear
  reference set" problem, and the motivation for graded
  [`ReferenceWeights`](@ref Platypus.ReferenceWeights).

- **Hannan, M. T., Pólos, L., & Carroll, G. R.** (2007). *Logics of
  Organization Theory: Audiences, Codes, and Ecologies*. Princeton
  University Press. Formal treatment of categories, typicality, and
  grade-of-membership in organisational populations.

## Density, outliers, and data depth

- **Breunig, M. M., Kriegel, H.-P., Ng, R. T., & Sander, J.** (2000).
  LOF: Identifying density-based local outliers. *Proceedings of the
  2000 ACM SIGMOD International Conference on Management of Data*,
  93–104. Introduces the Local Outlier Factor — a density score defined
  *relative to a point's own neighbours' density* — which is the
  `:lof` method and the default scoring rule, chosen because it needs no
  global density scale and so suits a heterogeneous collection.

- **Tukey, J. W.** (1975). Mathematics and the picturing of data.
  *Proceedings of the International Congress of Mathematicians*, 2,
  523–531. Introduces halfspace (Tukey) depth, the origin of the
  data-depth view of centrality used by the `:depth` method.

- **Liu, R. Y.** (1990). On a notion of data depth based on random
  simplices. *The Annals of Statistics*, 18(1), 405–414. Develops
  simplicial depth and the depth-based centre-outward ordering of a point
  cloud.

- **Zuo, Y., & Serfling, R.** (2000). General notions of statistical
  depth function. *The Annals of Statistics*, 28(2), 461–482. The
  unifying treatment of depth functions, including the Mahalanobis depth
  implemented by the `:depth` method.

- **Silverman, B. W.** (1986). *Density Estimation for Statistics and
  Data Analysis*. Chapman & Hall. Source of the rule-of-thumb bandwidth
  used by the `:kde` method.

## Copulas

- **Sklar, A.** (1959). Fonctions de répartition à n dimensions et leurs
  marges. *Publications de l'Institut de Statistique de l'Université de
  Paris*, 8, 229–231. Sklar's theorem: any joint distribution factors
  into its marginals and a copula. The theoretical licence for the
  copula-rank transform in [`build_pairs`](@ref Platypus.build_pairs),
  which isolates dependence structure from marginal shape.

- **Nelsen, R. B.** (2006). *An Introduction to Copulas* (2nd ed.).
  Springer. Standard reference for copula families and estimation; the
  basis for the planned [`fit_copula`](@ref Platypus.fit_copula)
  estimator.

## Graphs, random walks, and emergence

- **Lovász, L.** (1993). Random walks on graphs: A survey.
  *Combinatorics, Paul Erdős is Eighty*, 2, 1–46. The reference for the
  stationary distribution of a random walk as a density-like measure on
  a graph — the principled form of the "emergent reference set" idea
  reserved for [`fit_graph_typicality`](@ref Platypus.fit_graph_typicality).

- **Coifman, R. R., & Lafon, S.** (2006). Diffusion maps. *Applied and
  Computational Harmonic Analysis*, 21(1), 5–30. The diffusion-process
  view of density and geometry on a data graph; the continuous-density
  counterpart of the graph random-walk estimator.

- **Dorigo, M., & Stützle, T.** (2004). *Ant Colony Optimization*. MIT
  Press. The canonical statement of ant-colony optimization. Listed
  precisely to mark the boundary: ACO is a *combinatorial optimiser*, not
  a density estimator, and is therefore *not* used here — see the
  discussion on the [Home](index.md) page. The reusable kernel of the
  idea (pheromone reinforcement with evaporation) is subsumed by the
  graph random-walk model above.

## Sibling software

- **VacancyChains.jl** (Santoni, ongoing). Vacancy chain analysis in
  Julia. Used here as the structural template for package layout, naming
  conventions, and documentation. See
  [simoneSantoni/VacancyChains.jl](https://github.com/simoneSantoni/VacancyChains.jl).
