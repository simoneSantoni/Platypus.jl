# Julia engineering decisions

The review found two distinct storage defects: entity vectors could disagree with a cached matrix, and each feature was copied into every pair containing it. The revised design makes `EntitySet.matrix` the authoritative numerical store. Constructors copy caller data; `es[i]`, iteration, `es.entities`, `features`, and `feature_names` return detached snapshots. Direct field arrays remain available for advanced mutation, with dimensions, identities, and finiteness checked at analysis boundaries. A `struct` alone cannot guarantee immutable data: its array fields remain mutable in Julia. [Julia types manual](https://docs.julialang.org/en/v1/manual/types/#Composite-Types)

`build_pairs` creates one contiguous `Vector{Float64}` per transformed feature and shares it among feature pairs. A detached raw matrix preserves preprocessing provenance for weighted fitting. Coordinate storage therefore grows with entities times features rather than entities times feature pairs; pair descriptors still grow quadratically in the number of features. Concrete vector fields preserve estimator dispatch and specialization. Julia documents both the copying behavior of ordinary slices and the tradeoff between contiguous copies and views. Sharing each column once avoids repeated copies without introducing a new abstract coordinate field. [Julia performance tips](https://docs.julialang.org/en/v1.9/manual/performance-tips/#Consider-using-views-for-slices)

Weighted rank preprocessing uses the midpoint of each tied group's reference probability mass. Rows with zero reference weight remain scoreable queries but contribute no mass; queries outside the support may map to 0 or 1. Weighted z-scoring computes location and spread using only positive support, with the reliability-weight correction that recovers sample standard deviation for uniform weights. Weight normalization divides by the maximum before summation to avoid overflow; an unrepresentable positive support weight raises an error instead of silently disappearing. Finite-value validation is repeated after direct field mutations and before fitting.

The prior six-argument `PairSet` constructor remains available with `transform=:unknown` and no raw source matrix. A fit cannot reconstruct weighted marginal transformations from such coordinates. New `build_pairs` output stores both `transform` and `source_matrix`.

Compatibility changes are intentional: duplicate entity identifiers/schema names, invalid feature values, and mismatched metadata lengths now raise errors. Entity/accessor mutation no longer edits a collection indirectly. `pair_redundancy` now reports only descriptive absolute correlation; its previous suggested `1 - abs(correlation)` weight had no established interpretation as a correction for dependence between pair scores. Constant-feature correlations remain undefined (`NaN`).

Validation is in `test/data_integrity.jl`: input rejection, mutation isolation, finite-weight normalization, shared-storage identity, uniform-weight equivalence, weighted-transform scaling and zero-row exclusion, and EDA boundary cases. Statistical estimator validation is documented separately.

## Final integration review: representable numerical limits

A finite response of `1e308` produced infinite coefficients even in an intercept-only outcome regression: the orthogonal projection overflowed before dividing by the design norm. The outcome solver now scales the response before pivoted QR, forms fitted values and residuals in scaled units, and rescales final outputs. It rejects results whose coefficients, predictions or residuals cannot be represented as finite Float64 values. R² uses scaled residuals throughout.

The Gaussian KDE formerly subtracted opposite-sign values before dividing by bandwidth. For points `(-1e308,-1)` and `(1e308,1)` with bandwidths `(1e308,1)`, subtraction overflowed even though each standardized difference is exactly two. The kernel now falls back to dividing before subtracting only when ordinary subtraction overflows, preserving accuracy for nearby ordinary inputs. Explicit bandwidths are also checked after conversion to Float64; finite BigFloat values outside the representable positive Float64 range are rejected. Public-API regressions are in `test/numerical_limits.jl`.

Actual rendering was validated in an isolated temporary environment resolved offline from cached packages (Julia 1.10.12, Plots 1.41.7, GR 0.73.27). All nine plotting helpers and the graph-result atypicality overload produced SVG files successfully. RecipesBase recipe expansion is also exercised by package tests. Visual inspection of outcome and heatmap output identified label and color-scale problems that recipe-expansion success alone did not reveal.

The rendering issues were corrected at the recipe layer: empty subplot labels are omitted, and the outcome identity-line overlay explicitly repeats the main series labels because Plots applies defaults to every series. Heatmaps request a numeric colorbar despite suppressed series legends. User-supplied labels and colorbar settings remain overridable. `test/visualization.jl` checks these contracts, including detached plot data.


The optional integration script `validation/render_plots.jl` preserves the actual-backend check without adding Plots to the package dependencies. It renders all helpers to SVG and checks processed subplot labels and numerical color keys. From the repository root, reproduce in an isolated environment:

```sh
GKSwstype=100 julia --startup-file=no -e 'using Pkg; repo=pwd(); Pkg.activate(mktempdir()); Pkg.develop(path=repo); Pkg.add("Plots"); include(joinpath(repo,"validation","render_plots.jl"))'
```

The script prints its retained temporary output directory. To choose a directory with an already prepared environment, run `GKSwstype=100 julia --project=/path/to/environment validation/render_plots.jl /path/to/output`. `Pkg.offline(true)` before resolution reuses installed packages without network access when all dependencies are cached.

Final validation: the preserved script passed 23/23 actual-backend assertions and wrote all ten SVG files. Visual inspection confirmed that the outcome plot retained its title and Observed/Fitted axes and that score heatmaps displayed their numerical color scales. The numerical-limit regressions passed 28/28, and recipe-level plotting regressions cover 13 checks.
