# Platypus.jl

**Entity atypicality for Julia.**

How unusual is one entity, relative to a large and heterogeneous set of
its peers? A song among hundreds of thousands of songs; a firm among the
firms of a loosely bounded industry; a paper among a sprawling, fuzzily
delimited literature. Platypus.jl gives a measure of *atypicality* for
exactly this setting — and, crucially, for the setting in which the peer
group has **no clean boundary**.

The package takes its name from its subject. The platypus is the textbook
atypical entity: a mammal that lays eggs, has a bill, and detects prey
electrically. It is unusual not because any one feature is strange, but
because of the *combination* of features it carries. That is the idea
Platypus.jl operationalises.

## The problem: typicality without a reference set

Typicality is, by definition, relative. To say an entity is atypical is
to say it is unlike *the others* — and that immediately raises the
question: which others?

In a clean problem the reference set is given. Uzzi, Mukherjee, Stringer
& Jones (2013), in [*Atypical Combinations and Scientific
Impact*](references.md), had all of science to compare a paper against.
But many interesting populations do not come pre-bounded:

- A *music market* has fuzzy, overlapping genre boundaries. Is a song
  atypical "for pop", "for the market", or "for 2026"? There is no single
  right answer, and no crisp set to compute against.
- An *industry* is a category with graded, contested membership.
- A *literature* has no edge — only a denser core and a thinning
  periphery.

When the reference set is itself fuzzy, Platypus.jl makes two moves.

1. **Endogenous reference.** The observed collection *is* its own
   reference set. Typicality is read off the cloud of entities you have,
   not imported from an external baseline. Where membership is graded
   rather than binary, [`ReferenceWeights`](@ref) lets each entity
   contribute partially.
2. **Local reference.** Typicality is measured against *local* density —
   how crowded the neighbourhood around an entity is, compared to the
   neighbourhoods around nearby entities — not against one global
   baseline. A heterogeneous, lumpy collection then stops being a problem
   to be corrected and becomes the structure the measure adapts to.

## The approach: typicality of pairs of features

Platypus.jl does not score an entity in the full ``n``-dimensional
feature space. High-dimensional density estimates are fragile, and a
single "this entity is weird" number is hard to interpret. Instead,
following Uzzi et al. (2013), it decomposes the entity into the
``\binom{n}{2}`` **pairs of features**.

Uzzi et al. studied the combinations of *journals* a paper cited: for
each pair of journals co-cited by a paper, they measured how often that
pair is co-cited across all of science, relative to chance. A paper then
has a whole *distribution* of pairwise typicalities, and they summarised
it by its **median** (how conventional the paper is overall) and its
**10th percentile** (how striking its single most novel combination is).

Platypus.jl carries this design over from discrete co-citation counts to
**continuous features**. For a song with features danceability,
acousticness, energy, valence, …, each *pair* of features defines a
two-dimensional cloud — one point per song. The focal song's position in
that cloud is scored for local typicality; do this for every pair; then
summarise the per-entity distribution of pairwise scores.

```
    feature pair (danceability, acousticness)
                                                 v (acousticness rank)
    1.0 ┤                          • focal song  │  ← off the dependence
        │        • • • •                         │     ridge: atypical pair
        │      • • • • • •                       │
        │    • • • • • •                         │
        │      • • • •                           │
    0.0 ┤  • •                                   │
        └────────────────────────────────────────
        0.0                              1.0  u (danceability rank)
```

### Step 1 — copula-rank each feature

Each feature is replaced by its rank, rescaled to ``(0, 1)``:

```math
u_i \;=\; \frac{\operatorname{rank}(x_i) - \tfrac{1}{2}}{n}.
```

This is the empirical *copula* transform. It makes every feature
marginally uniform, so the only thing left in a feature pair's cloud is
the **dependence structure** of the two features. An entity then scores
as atypical because of an unusual *combination* — danceable yet acoustic
— not merely because one feature is extreme. [`build_pairs`](@ref)
applies this transform by default.

### Step 2 — score local typicality on each pair

For each feature pair, every entity gets a local-density typicality
score. Platypus.jl offers four rules, selected by the `method` keyword of
[`fit_typicality`](@ref):

| Method | Idea | Best when |
|--------|------|-----------|
| `:knn`   | inverse distance to the ``k``-th nearest neighbour | a fast, assumption-free baseline is wanted |
| `:lof`   | density *relative to neighbours'* density (Breunig et al. 2000) | the cloud's density varies from region to region |
| `:kde`   | smooth Gaussian kernel density | a smooth estimate is wanted; supports weights |
| `:depth` | Mahalanobis data depth — centre-outward ordering | a distribution-free central-tendency score is wanted |

`:lof` is the natural default for a genuinely heterogeneous set: it
compares an entity's local density to that of its own neighbours, so it
has no single global density scale to be misled by.

### Step 3 — summarise the pairwise distribution

Each entity now has ``\binom{n}{2}`` pairwise scores. [`fit_typicality`](@ref)
collapses them into the Uzzi-style triple:

```math
\textit{conventionality}_e = \operatorname{median}_p\, s_{ep}, \qquad
\textit{novelty}_e = Q_{0.1}\big(s_{e\bullet}\big), \qquad
\textit{dispersion}_e = \operatorname{IQR}_p\, s_{ep},
```

where ``s_{ep}`` is entity ``e``'s typicality on pair ``p``.
*Conventionality* is how typical the entity is overall; *novelty* is how
striking its most unusual single combination is; *dispersion* tells you
whether the entity mixes very conventional and very unconventional pairs.
The headline quantity, [`atypicality`](@ref), is ``1 - \textit{conventionality}``.

Because the three summaries are separate, they can be entered separately
into a model of an outcome (song popularity, scientific impact) — which
is how Uzzi et al. found that the highest-impact papers are *both* highly
conventional *and* carry a tail of novelty.

## A worked example

```julia
using Platypus

# An n_entities × n_features matrix of continuous features.
es = EntitySet(rand(500, 8),
               [:danceability, :energy, :acousticness, :valence,
                :tempo, :loudness, :speechiness, :liveness])

# Project onto all 8*7/2 = 28 feature pairs (copula-rank transform).
pairs = build_pairs(es)

# Score and aggregate.
result = fit_typicality(pairs; method=:lof)

summarytable(result)        # conventionality / novelty / dispersion / atypicality
atypicality(result)         # the headline per-entity vector
```

See [Getting Started](getting_started.md) to reproduce this end to end.

## Why not an ant-colony algorithm?

A reasonable instinct, when the reference set is fuzzy and self-defining,
is to reach for an *emergent, self-organising* method — an ant-colony
algorithm, say. The instinct is sound; ant-colony **optimization** is
not the right instrument.

Ant Colony Optimization is a *combinatorial optimiser*: it searches for
good paths in a graph by having ants reinforce edges with pheromone.
Typicality estimation has no objective function to optimise. Strip ACO
down to the mechanism that *is* relevant here — pheromone deposition with
evaporation — and what remains is a decaying frequency count, i.e. an
online kernel density estimate, useful only if the population drifts over
time.

The genuinely useful version of the "emergent, no imposed reference"
intuition is a **reinforced random walk on a similarity graph**: take
typicality to be the stationary distribution of a random walk over a
graph of entities, so dense, well-connected regions accumulate
probability mass on their own. That is the design reserved for
[`fit_graph_typicality`](@ref) — see the [Estimation API](api/estimation.md).

## Notation used in this manual

| Symbol | Meaning |
|--------|---------|
| ``n`` | Number of features per entity |
| ``\binom{n}{2}`` | Number of feature pairs, ``n(n-1)/2`` |
| ``x_i`` | Raw value of feature ``i`` |
| ``u_i`` | Copula-rank transform of feature ``i`` |
| ``s_{ep}`` | Typicality of entity ``e`` on feature pair ``p`` |
| ``Q_q(\cdot)`` | The ``q``-quantile of a set of values |
| ``\textit{atypicality}`` | ``1 - \textit{conventionality}`` |

## Where to next

- [Getting Started](getting_started.md) — install the package and run the
  example above end to end.
- [Entities and Data](guide/entities.md) — how a collection is
  represented as `Entity` and `EntitySet`.
- [Building Feature Pairs](guide/pairs.md) — the copula-rank transform
  and [`build_pairs`](@ref).
- [Estimation](api/estimation.md) — the full statistical specification of
  [`fit_typicality`](@ref), the four scoring rules, and the reserved
  estimators.
- [References](references.md) — annotated bibliography, starting from
  Uzzi et al. (2013).
