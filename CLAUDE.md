# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What the package is

Platypus.jl measures **entity atypicality** in large, heterogeneous sets — sets whose reference population has no clean boundary. It decomposes each entity into the `n(n-1)/2` *pairs* of its features and scores the local typicality of each pairwise combination, following Uzzi et al. (2013), *Atypical Combinations and Scientific Impact*. The package is scaffolded but partial: data loading, pair construction, EDA, and the core estimator `fit_typicality` are implemented; `fit_graph_typicality`, `fit_outcome`, `fit_copula`, and all `plot_*` functions are exported stubs that raise `ErrorException("... is not implemented yet")` so the public API in [README.md](README.md) is stable from v0.1.0 onward.

## Development commands

```bash
# Run the full test suite from the shell
julia --project -e 'using Pkg; Pkg.test()'

# Or, faster while iterating, run runtests.jl directly
julia --project -e 'include("test/runtests.jl")'

# Build the docs (separate environment under docs/)
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'
```

The test suite is a single [test/runtests.jl](test/runtests.jl) with named `@testset` blocks. To run a single testset, copy its body into the REPL with the package loaded — there is no per-file test entry point.

## Architecture

The module entry point [src/Platypus.jl](src/Platypus.jl) re-exports the public surface and `include`s files in order: `types.jl`, `entities.jl`, `pairs.jl`, `eda.jl`, `estimation.jl`, `visualization.jl`. This ordering matters — later files use types and helpers from earlier ones.

### Core types ([src/types.jl](src/types.jl))

- **`Entity`** — one collection member: `id`, `features::Vector{Float64}`, optional `label`. The `label` is a soft category tag (genre, sector) and is **never** used to define the reference set; it is carried through only so results can be grouped. The inner constructor enforces at least one feature.

- **`EntitySet`** — entities sharing one feature schema, with a cached `n_entities × n_features` matrix. It is both the unit of analysis and, since the reference set is endogenous, the reference set itself. Constructor validates that every entity has `length(feature_names)` features.

- **`FeaturePair`** — the projection of an `EntitySet` onto two features: a 2-D cloud with coordinates `u`, `v` (one point per entity).

- **`PairSet`** — the set of all `n(n-1)/2` `FeaturePair`s; the object `fit_typicality` consumes.

- **`ReferenceWeights`** — per-entity, non-negative graded membership in the reference set, for the fuzzy-boundary case. Honoured by the `:kde` and `:depth` methods only.

- **`TYPICALITY_METHODS = (:knn, :lof, :kde, :depth)`** — the closed set of scoring rules accepted by `fit_typicality`.

### Data loading ([src/entities.jl](src/entities.jl))

`load_entities(df; ...)` builds an `EntitySet` from a `DataFrame` (or CSV path). By default every numeric column other than `id_col`/`label_col` becomes a feature. `parse_feature` coerces cells to `Float64` and **rejects `missing`** — the density estimators have no missing-value semantics, so impute or drop before loading.

### Pair construction ([src/pairs.jl](src/pairs.jl))

`build_pairs(es; transform=:rank)` applies a per-feature transform, then enumerates the `n(n-1)/2` pairs. `:rank` is the empirical copula transform (`(tiedrank(x) - 0.5)/n`) and is the default — it strips marginal rarity so pair scores reflect the *combination*. Other transforms: `:zscore`, `:none`. Unknown transforms `throw(ArgumentError)`. Fewer than two features also throws.

### EDA ([src/eda.jl](src/eda.jl))

`feature_summary`, `rank_matrix`, `pair_correlations`, `pair_redundancy`, `nearest_entities`, `compute_entity_statistics`, `pair_table`. `pair_correlations`/`pair_redundancy` exist to expose the caveat that the `n(n-1)/2` pairwise scores are not an i.i.d. sample when features are correlated.

### Estimation ([src/estimation.jl](src/estimation.jl))

- **`fit_typicality(ps; method, k, novelty_quantile, normalize, weights)`** — scores each entity's local typicality on every feature pair (`:knn`/`:lof`/`:kde`/`:depth`), optionally rank-normalises per pair so scores are comparable, then aggregates per entity into `conventionality` (median), `novelty` (low quantile), `dispersion` (IQR), `atypicality` (`1 - conventionality`). This is the Uzzi et al. (2013) summary triple, generalised from discrete co-citation z-scores to continuous-feature densities.
- **`TypicalityResult`** — exposes `atypicality`, `conventionality`, `novelty`, `dispersion`, `summarytable` (per-entity `DataFrame`), `scoretable` (long-form entity × pair).
- **`fit_graph_typicality`, `fit_outcome`, `fit_copula`** — exported but raise `ErrorException` until implemented. Reserve their existing signatures.

### Visualization ([src/visualization.jl](src/visualization.jl))

All `plot_*` functions are exported stubs that raise. Implement them as Plots.jl recipes (or CairoMakie equivalents) consuming the `DataFrame`/matrix shapes that EDA and estimation already produce. The names and arities are the public API — don't rename them when implementing.

### Documentation site ([docs/](docs/))

Docs are built with `Documenter.jl` from [docs/make.jl](docs/make.jl). The docs environment has its own `Project.toml` and pulls the package in via a `[sources]` path entry. The page tree is committed in `make.jl`; if you add a new page, both add the file under `docs/src/` and register it in the `pages = [...]` block.

## Working conventions

- **The reference set is endogenous.** Every estimator scores an entity against the other members of the same `EntitySet`. Do not add an API that takes an external reference population; graded membership goes through `ReferenceWeights`.
- **Labels are never structural.** `Entity.label` is descriptive only. Code that computes typicality must not branch on it.
- **The copula-rank transform is the default** and the reason pair scores measure *combination* atypicality. Don't change `build_pairs`'s default away from `:rank`.
- **`TYPICALITY_METHODS` is a closed set.** Adding a method means updating the tuple, `_pair_typicality`, and the `fit_typicality` docstring.
- **Stub `fit_*` and `plot_*` are intentional.** Implement by replacing the body; do not change the exported name or arity.

## Template

The package layout, naming conventions (`load_*`/`build_*`/`fit_*`/`plot_*`, per-feature files, single `runtests.jl`, separate `docs/Project.toml`), and documentation structure are modeled on **VacancyChains.jl** by the same maintainer. Mirror its conventions when extending Platypus.jl.
