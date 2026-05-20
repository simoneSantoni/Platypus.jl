# Types

The data model has three layers:

1. **Atomic records** — an [`Entity`](@ref) is one member of a
   collection: an id, a vector of continuous feature values, and an
   optional category label.
2. **Collections** — an [`EntitySet`](@ref) is the full set of observed
   entities sharing one feature schema. It is both the unit of analysis
   *and*, in the absence of an external population, the reference set.
3. **Projections** — a [`FeaturePair`](@ref) is the collection projected
   onto two of its features; a [`PairSet`](@ref) is the set of all
   ``n(n-1)/2`` such projections, the object [`fit_typicality`](@ref)
   consumes.

A fourth helper, [`ReferenceWeights`](@ref), attaches a graded
membership weight to each entity, for the case where the boundary of the
reference set is fuzzy rather than crisp.

The split into "entities" (raw observations) and "pairs" (projections)
follows the conceptual move made by
[Uzzi et al. (2013)](../references.md): atypicality is not read off a
single high-dimensional point but assembled from the typicality of
*pairs* of attributes. [Rosch & Mervis (1975)](../references.md) is the
older source for treating typicality as a graded, feature-configuration
property of category members rather than a binary one.

## Why the reference set is endogenous

Platypus.jl never takes a reference population as an external input. The
`EntitySet` you provide is the reference set: every estimator scores an
entity relative to the other members of the same collection. When
membership is graded — overlapping genres, contested industry
boundaries, a literature with no edge — [`ReferenceWeights`](@ref) lets
each entity contribute *partially* to the reference cloud instead of
being forced in or out. This is the design choice that lets the package
operate when, in the language of the [Home](../index.md) page, "there is
no clear reference set". The categorisation literature
([Hsu, Hannan & Koçak 2009](../references.md)) reaches the same
conclusion from the side of markets with fuzzy category boundaries.

## Data types

```@docs
Entity
EntitySet
FeaturePair
PairSet
ReferenceWeights
TYPICALITY_METHODS
```

## Construction and projection

```@docs
load_entities
build_pairs
```

## Exploratory analysis

These return summary scalars, matrices, or `DataFrame`s suitable for
direct plotting. The redundancy diagnostics — [`pair_correlations`](@ref)
and [`pair_redundancy`](@ref) — address the caveat that the ``n(n-1)/2``
pairwise scores are not an i.i.d. sample when features are correlated.

```@docs
feature_summary
rank_matrix
pair_correlations
pair_redundancy
nearest_entities
compute_entity_statistics
pair_table
```

## Predicates and accessors

```@docs
features
label
feature_names
n_entities
n_features
n_pairs
```

## Internal helpers

```@docs
Platypus.parse_feature
```

## References

The conceptual choices encoded in these types — atypicality as an
aggregate over feature *pairs*, the reference set as endogenous and
graded, typicality as local density — come from the works listed on the
[References](../references.md) page. The canonical entry points are
Uzzi et al. (2013) and Breunig et al. (2000).
