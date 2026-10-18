# Communication and research artefact decisions

Date: 2026-10-18. Scope: documentation, reproducible examples, public contracts,
release communication, and practical interpretation. Human-perception validation
is explicitly outside this phase; no claims of such validation are made.

## Procedure

1. Read README, all user guides and API pages, docstrings, Project.toml, docs
   configuration, CI workflows, and the panel review. Compare advertised behavior
   with implementation and tests rather than treating documentation as evidence.
2. Re-read the integrated 0.2 source for raw/normalized scales, weighting,
   transformations, graph component conventions, copula support restrictions,
   OLS assumptions, plotting recipes, and snapshot ownership.
3. Search primary official documentation using “Documenter @example blocks”
   and “RecipesBase recipes user recipe”. Inspect the relevant official pages,
   recorded below. These informed executable examples and the rendering contract.
4. Rewrite package positioning, introductory explanations, user/API pages, and
   migration notes around implemented behavior. Distinguish method definitions,
   package conventions, and empirical validation.
5. Replace the unseeded promised-outlier example with a deterministic 40-row,
   three-feature fixture. Execute the default LOF pipeline and assert that row 40
   has maximal headline atypicality. Convert core guide code to named Documenter
   examples so later builds execute assertions.
6. Execute tutorial/guide blocks against the shared Julia 1.10.12 environment,
   then integrate with the root agent's strict Documenter build. The build owner
   controls docs deployment routing and export coverage.

## Primary documentation and decisions

Documenter's official syntax manual states that `@example` evaluates code and
inserts its result. Named example blocks share a module on a page. We use those
blocks for executable workflows and assertions, reserving ordinary Julia fences
for installation and optional rendering instructions that should not modify the
build environment. A successful syntax-only code fence is not an execution test.
[Documenter syntax manual](https://documenter.juliadocs.org/stable/man/syntax/).

RecipesBase supports recipes that can be consumed by plotting packages. We
document `TypicalityPlot` as a data-bearing recipe object, and show the explicit
rendering call `using Plots; plot(plot_...(args))`. This keeps the package's
backend-independent construction contract visible instead of suggesting that
each helper immediately displays a figure.
[RecipesBase recipe types](https://docs.juliaplots.org/stable/RecipesBase/types/),
[Plots recipes overview](https://docs.juliaplots.org/stable/recipes/).

## Communication changes

- Position the package as descriptive measurement relative to a selected
  collection, not a solution that removes reference-population choices.
- Distinguish local relative density, density, global centrality, fitted copula
  dependence, and full-feature graph connectivity.
- Define `novelty` as lower-tail typicality with the direction prominently stated.
  Quantiles are not minima; IQR is not sampling uncertainty.
- Explain that rank transforms and mid-CDF score normalization are not null-model
  calibration. Attribute median/lower-tail inspiration without attributing the
  package's IQR to Uzzi et al.
- Explain weighted transformations, zero-weight queries, unsupported weighted
  methods, singularity choices, and raw log-density zeros.
- Explain plotting dependencies and ID-versus-row selectors.
- Remove placeholder/stub claims after implementation; describe the 0.2 working
  tree without claiming an unpublished release exists.
- Document migration, algorithmic costs, snapshot ownership, and metadata needed
  to interpret exported tables.

## Verification and limits

The default-LOF fixture was executed on Julia 1.10.12: row 40 had atypicality
0.9875; the largest score among rows 1–39 was 0.95. This verifies the specific
worked example, not broad performance. Assertions remain in the documentation.

All 19 named example blocks across Getting Started and the four user guides
executed successfully as separate page modules on Julia 1.10.12. An identifier
collision with DataFrames.metadata was removed by naming the saved configuration
fit_metadata. The integrated strict Documenter build provides the final check
of rendering, links, and public API coverage.

The completed implementation's numerical tests, synthetic studies, and benchmark
records are reported by their respective owners. None measures perceived
atypicality in human observers. The documentation therefore makes no claim to
that construct validation. Stable release routing and dependency environments
are controlled by the root integration changes; check the final build record
before publishing the documentation.
