# Julia 1.13.2 default migration

Requested on 2026-10-18 after the multidisciplinary implementation review.
The user's “Julia 13.2” is interpreted as the installed Julia 1.13.2 release.

## Configuration

- `.julia-version` records 1.13.2 for development.
- The primary CI job and documentation build use the exact 1.13.2 version.
  The minimum-version CI job and package compatibility floor remain Julia 1.9.
  Choosing a default runtime does not require dropping tested older support.
- A Juliaup directory override selects 1.13.2 for `/home/simon/Platypus.jl`
  and its subdirectories. Each additional checkout needs its own override;
  `.julia-version` alone does not configure Juliaup. The global default channel
  was not changed.
- README, getting-started instructions, repository guidance and changelog now
  explain the default and setup procedure.

The Juliaup commands and precedence were checked against the
[official Juliaup documentation](https://github.com/JuliaLang/juliaup#overrides).
Exact workflow version selection follows the
[setup-julia documentation](https://github.com/julia-actions/setup-julia#julia-versions).

## Dependency migration and investigation

1. Confirmed the installed binary reports Julia 1.13.2 and can load its standard
   libraries in the default environment.
2. Attempted resolution and normal tests with the previous Julia 1.10 manifests.
   These failed during Markdown/Pkg precompilation with annotated-string stack
   overflow; a non-precompiled diagnostic also found Statistics pinned to 1.10.0,
   conflicting with the available 1.11 series.
3. Preserved the old root and docs manifests as `Manifest-v1.10.toml` and resolved
   fresh environments under 1.13.2. Fresh root resolution selected Statistics
   1.11.5, RecipesBase 1.4.0, PrecompileTools 1.3.4 and the runtime's own
   StyledStrings 1.11.0. Pkg and Markdown then loaded normally in the project.
   This supports stale environment pins as the cause; no standard-library code
   was patched and no global cache was deleted.
4. Stored the new root and docs manifests as `Manifest-v1.13.toml`. Generated manifests
   remain ignored by git, consistent with this package's existing convention.
   The procedure follows [Pkg's version-specific manifest documentation](https://pkgdocs.julialang.org/v1/toml-files/#Different-Manifests-for-Different-Julia-versions).
5. Ran full package tests with normal compilation and two Julia threads,
   plus resolution and strict documentation validation under 1.13.2.

Older Julia 1.9/1.10 benchmark and simulation artifacts were not relabeled as
new-runtime measurements. This migration changes runtime selection and
environment resolution, not the definitions of statistical scores.

## Validation

- `julia --version` reports **1.13.2** in both the repository root and `docs/`.
- `Pkg.test()` passed **911/911 checks** on Julia 1.13.2 with two threads and
  normal compiled modules. No numerical implementation changes were required.
- The strict Documenter build passed, including executed examples, exported
  docstrings and cross-references, with normal compiled modules.
- These full checks used the writable validation depot layered over the installed
  depot. Resolution succeeded from cached packages despite a registry-network
  warning in the restricted shell. The ordinary installed depot was then
  instantiated, and package loading plus a default fit passed using plain
  `julia --startup-file=no --project=.` with no `JULIA_DEPOT_PATH` override.
  Daily use does not require the validation depot.

Reproduction from this configured checkout:

```sh
julia --threads=2 --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```
