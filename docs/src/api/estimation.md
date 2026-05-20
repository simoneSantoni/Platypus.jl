# Estimation API

Platypus.jl provides one implemented estimator and three reserved names.
The implemented estimator, [`fit_typicality`](@ref Platypus.fit_typicality),
scores every entity's local typicality on each feature pair and
aggregates the ``n(n-1)/2`` scores per entity. The reserved names —
[`fit_graph_typicality`](@ref Platypus.fit_graph_typicality),
[`fit_outcome`](@ref Platypus.fit_outcome),
[`fit_copula`](@ref Platypus.fit_copula) — currently raise
`not implemented`; they are documented here so callers can plan against a
stable public API.

## `fit_typicality`: pairwise local-density atypicality

### Statistical model

Let an entity carry ``n`` continuous features. Each feature is first
replaced by its empirical copula rank,

```math
u_i = \frac{\operatorname{rank}(x_i) - \tfrac{1}{2}}{N},
```

over the ``N`` entities of the collection, so every feature is marginally
uniform on ``(0,1)``. The collection is then projected onto each of the
``\binom{n}{2}`` feature pairs; pair ``p = (i, j)`` is a cloud of ``N``
points ``\{(u_i^{(e)}, u_j^{(e)})\}_e``.

For each pair, every entity ``e`` receives a **local-density typicality
score** ``s_{ep}``. With the default rank-normalisation, the raw density
estimates of a pair are replaced by their ranks rescaled to ``(0, 1]``,
so that scores are comparable across pairs before they are aggregated —
the continuous-feature analogue of the ``z``-scoring that
[Uzzi et al. (2013)](../references.md) applied to discrete co-citation
counts.

The four scoring rules differ in how the local density is estimated:

- **`:knn`** — let ``d_k(e)`` be the distance from entity ``e`` to its
  ``k``-th nearest neighbour in the pair cloud. The density estimate is
  ``\widehat{f}(e) \propto 1 / d_k(e)``.

- **`:lof`** — the Local Outlier Factor of
  [Breunig et al. (2000)](../references.md). Writing ``\mathrm{lrd}_k(e)``
  for the local reachability density of ``e`` and ``N_k(e)`` for its
  ``k``-neighbourhood,

  ```math
  \mathrm{LOF}_k(e) = \frac{1}{|N_k(e)|}
      \sum_{o \in N_k(e)} \frac{\mathrm{lrd}_k(o)}{\mathrm{lrd}_k(e)},
  ```

  which compares ``e``'s density to that of *its own neighbours*.
  ``\mathrm{LOF} \approx 1`` marks an inlier; ``\mathrm{LOF} \gg 1`` an
  outlier. Typicality is taken as ``1 / \mathrm{LOF}``. Because LOF is
  defined relative to local density, it needs no single global density
  scale — the property that makes it the default for a heterogeneous
  collection.

- **`:kde`** — a weighted two-dimensional Gaussian kernel density,

  ```math
  \widehat{f}(e) = \frac{1}{\sum_{o \ne e} w_o \, 2\pi h_u h_v}
      \sum_{o \ne e} w_o \,
      \exp\!\Big(\!-\tfrac{1}{2}\big[(\Delta u/h_u)^2 + (\Delta v/h_v)^2\big]\Big),
  ```

  with per-dimension Silverman bandwidths ``h_u, h_v`` and reference
  weights ``w_o``.

- **`:depth`** — the Mahalanobis data depth
  ([Liu 1990](../references.md); [Zuo & Serfling 2000](../references.md)).
  With weighted mean ``\mu`` and weighted covariance ``\Sigma`` of the
  pair cloud,

  ```math
  \mathrm{depth}(e) = \frac{1}{1 + (x_e - \mu)^\top \Sigma^{-1} (x_e - \mu)},
  ```

  a centre-outward ordering that assumes no particular distribution.

### Aggregating the pairwise scores

Each entity ends with a distribution of ``\binom{n}{2}`` scores, which is
collapsed into the [Uzzi et al. (2013)](../references.md) summary triple:

```math
\textit{conventionality}_e = \operatorname{median}_p s_{ep}, \quad
\textit{novelty}_e = Q_{q}\big(s_{e\bullet}\big), \quad
\textit{dispersion}_e = \operatorname{IQR}_p s_{ep},
```

where ``q`` is `novelty_quantile` (default ``0.1``, matching Uzzi et al.'s
10th-percentile "novelty"). The headline quantity is
``\textit{atypicality}_e = 1 - \textit{conventionality}_e``.

Keeping the three summaries separate is deliberate: Uzzi et al. found
that the highest-impact papers scored *high on conventionality and high
on tail novelty at once*, a pattern invisible to any single atypicality
number. The three can be entered separately into a model of an outcome
(see [`fit_outcome`](@ref Platypus.fit_outcome)).

### Theoretical context

The pairwise decomposition is inherited from
[Uzzi et al. (2013)](../references.md), who measured the typicality of
*pairs* of journals in a paper's reference list against a randomised
null. Platypus.jl carries the design from discrete co-citation counts to
continuous features, which requires two substitutions: a *null model* —
here the empirical copula transform, which removes marginal rarity so
only the combination is scored — and a *density estimator* in place of a
co-occurrence ``z``-score.

The use of *local* density (rather than a global baseline) is what makes
the measure robust when the collection is heterogeneous and the reference
set has no clean boundary; [Breunig et al. (2000)](../references.md) is
the source of the local-density-ratio idea. The applied target —
typicality of cultural products and its non-monotone effect on success —
is the subject of [Askin & Mauskapf (2017)](../references.md).

### Reference-set weighting

When the reference population is fuzzy, [`ReferenceWeights`](@ref Platypus.ReferenceWeights)
let each entity contribute partially to the reference cloud. Weights are
honoured by `:kde` (weighted kernel sum) and `:depth` (weighted mean and
covariance). The neighbour-based rules `:knn` and `:lof` ignore weights —
they have no natural weighted form — which the docstring states
explicitly.

```@docs
fit_typicality
```

## `TypicalityResult` and accessors

`TypicalityResult` exposes the full pairwise score matrix and the
per-entity summary vectors. The `atypicality` / `conventionality` /
`novelty` / `dispersion` accessors return the summary vectors;
`summarytable` returns a per-entity `DataFrame`; `scoretable` returns the
long-form entity × feature-pair table, suitable for joining against
entity-level covariates or for plotting.

```@docs
TypicalityResult
atypicality
conventionality
novelty
dispersion
summarytable
scoretable
```

## Reserved estimators (not yet implemented)

These names are exported so downstream code can bind against them today;
they currently raise `ErrorException` with a "not implemented" message.
The intended specifications are below.

### `fit_graph_typicality`

The endogenous, graph-based estimator for the genuinely fuzzy-reference
case. It builds a similarity graph over entities and takes typicality to
be the stationary distribution of a random walk on it: dense,
well-connected regions accumulate probability mass on their own, so the
reference "density" is *emergent* rather than imposed. This is the
rigorous form of the self-organising intuition discussed on the
[Home](../index.md) page — and the reason an ant-colony *optimiser* is
the wrong tool, since a reinforced random walk, not a path search, is
what the problem calls for. [Lovász (1993)](../references.md) is the
reference for random walks on graphs;
[Coifman & Lafon (2006)](../references.md) for the diffusion-map view.

```@docs
fit_graph_typicality
```

### `fit_outcome`

A regression of an entity-level outcome — song popularity, scientific
impact, firm performance — on the typicality summaries
`conventionality`, `novelty`, and `dispersion`. This is the estimator
that tests the [Uzzi et al. (2013)](../references.md) hypothesis
directly: that the outcome is highest when conventionality and tail
novelty are *both* high, a non-monotone pattern also reported for music
by [Askin & Mauskapf (2017)](../references.md). Because success
influences what later becomes typical, the intended specification
measures typicality against a strictly prior window.

```@docs
fit_outcome
```

### `fit_copula`

A parametric alternative to the non-parametric scoring rules: fit a
copula to each feature pair and score typicality from the fitted copula
density. Because the copula is by construction marginal-free
([Sklar 1959](../references.md); [Nelsen 2006](../references.md)), this
gives a smooth combination-atypicality score without the bandwidth or
neighbourhood-size choices the non-parametric rules require.

```@docs
fit_copula
```

## References

The estimator design follows [Uzzi et al. (2013)](../references.md) for
the pairwise decomposition and the conventionality/novelty summary,
[Breunig et al. (2000)](../references.md) for the local-density-ratio
score, [Liu (1990)](../references.md) and
[Zuo & Serfling (2000)](../references.md) for data depth,
[Sklar (1959)](../references.md) for the copula transform, and
[Lovász (1993)](../references.md) for the graph random-walk model. Full
citations are on the [References](../references.md) page.
