# Statistical methods: search and adoption notes

Review date: 2026-10-18. Scope: pairwise descriptive atypicality, weighting,
bounded-support smoothing, and a feasible parametric copula implementation.
These notes distinguish literature results from implementation proposals.
Human perception validation is outside this revision's scope. No claim of
universal anomaly calibration or validated social meaning follows from the
numerical methods below.

## Procedure and sources

The screen started with the review's reproduced defects: row-dependent LOF
ties; zero-weight leakage into bandwidths and ranks; ordinary Gaussian KDE
boundary depression on an independent uniform grid; and an empty weighted
leave-one-out reference silently converted to zero density. Searches included
“empirical beta copula Segers Sibuya Tsukahara”, “Geenens Charpentier Paindaveine
copula density local likelihood”, “mirror reflection copula kernel estimator”,
“LoOP Local Outlier Probabilities”, and “INFLO symmetric neighborhood”.
Primary papers and official implementations/documentation were preferred to
search summaries. This is a targeted engineering screen, not a systematic
literature review or an exhaustive benchmark.

| Candidate | Evidence inspected | Decision for this revision |
|---|---|---|
| Empirical beta copula | Segers, Sibuya and Tsukahara, sections 2.1–2.2 | Adopt strictly for equally weighted, positive, untied ranks; reject unsupported inputs. |
| Probit/local-likelihood copula KDE | Geenens, Charpentier and Paindaveine, transformation identity, sections 3–4 | Defer full method: requires actual local likelihood and bandwidth selection; a plain probit KDE is not equivalent. |
| Reflected KDE | Muia, Atutey and Hasan, definition 2.1; direct support-integral derivation below | Adopt the boundary-reflection principle for rank coordinates, with explicit kernel/support choices. |
| Gaussian copula | Official `copula` density documentation and normal-copula source | Recommend a narrowly named parametric family, with transparent parameter estimation and limitations. |
| Corrected LOF | Breunig et al., k-distance neighborhood definition | Retain and repair ties; establish oracle and permutation tests first. |
| LoOP | Authors' publication record and indexed abstract | Defer: probability interpretation would require additional assumptions and validation; full PDF retrieval redirected. |
| INFLO | Jin et al., influence-space construction | Defer: reverse-neighbor neighborhoods address a different modeling choice, not the current numerical defects. |

### Empirical beta copula

For untied ranks `R_ij` in `1:n`, its density is

```text
c_beta(u,v) = (1/n) sum_i beta_pdf(u; R_i1,n+1-R_i1)
                            beta_pdf(v; R_i2,n+1-R_i2).
```

The equal mixture has uniform margins because each rank occurs once in each
coordinate. The paper explicitly identifies ties as breaking that argument
and discusses random tie breaking. Its finite-sample comparisons concern
copula estimation; they do not establish superior in-sample anomaly rankings.
[Segers, Sibuya and Tsukahara, original paper](https://arxiv.org/pdf/1607.04430).

Engineering assessment: replacing `1/n` with arbitrary membership weights
does not preserve the rank-mixture identity. Removing one mixture component
without recomputing ranks likewise destroys it. Random tie breaking would
add artificial distinctions to rounded or ordinal measurements. The final
implementation therefore restricts this family to equal strictly positive
weights and untied marginal ranks, and evaluates the complete in-sample
mixture. It rejects ties, unequal weights, and zero-weight queries. This is
a restricted advanced alternative, not a weighted general-purpose estimator.

### Transformation and local likelihood

With `s=Phi^-1(u)` and `t=Phi^-1(v)`, the transformation identity is
`c(u,v)=f_ST(s,t)/(phi(s)*phi(t))`. The improved estimators use local log-linear
or log-quadratic density estimation, and investigate adaptive bandwidths.
Their paper distinguishes these from a naive transformed kernel estimator.
[Geenens, Charpentier and Paindaveine, original paper](https://arxiv.org/pdf/1404.4414).

Engineering assessment: the Jacobian can amplify numerical tail error.
Endpoints, singular local fits, bandwidth selection, and comparison with an
independent implementation all need coverage. Defer a faithful implementation
instead of calling a simpler transformation the published improved method.

### Reflection: concrete implementable option

The mirror estimator augments each coordinate `x` with `x`, `-x`, and `2-x`,
giving nine bivariate images. The bounded-density assumptions in the screened
paper do not cover arbitrary corner-singular copulas.
[Muia, Atutey and Hasan, definition 2.1](https://arxiv.org/pdf/2502.05470).

The following is an engineering derivation, not a claim that the complete
paper's estimator or bandwidth optimizer has been adopted. For a symmetric
compact kernel `K` supported on `[-1,1]` and `0<h<=1`, define

```text
R_h(u,x) = [K((u-x)/h) + K((u+x)/h) + K((u-2+x)/h)] / h
K(z) = 3/4 * (1-z^2) * indicator(abs(z)<=1)     # Epanechnikov option
f_-i(u,v) = sum_{j!=i} p_j R_hu(u,u_j) R_hv(v,v_j) / (1-p_i).
```

Reflection folds exactly the kernel mass outside either endpoint into the
unit interval, so `integral_0^1 R_h(u,x) du = 1` for `x` in `[0,1]`.
The product mixture therefore integrates to one if `1-p_i>0`. It need not
have exactly uniform estimated margins, so it is a bounded rank-space KDE,
not automatically a genuine finite-sample copula density. Compact support
can yield true zero density for isolated queries; define ranking ties and
avoid fabricating positive evidence with an arbitrary denominator epsilon.

For an unbounded Gaussian kernel, three images leave a small but nonzero
unrepresented tail. Either sum further periodic reflections to a stated
tolerance, or explicitly normalize each reflected kernel by its exact
unit-interval integral. Do not claim exact mass conservation for an
unnormalized three-image Gaussian. Keep ordinary real-line smoothing for
`:none` and `:zscore`; support treatment must use stored transform metadata.

Use positive reference weights only when computing location, scale and
effective sample size. With normalized weights `p`, `n_eff=1/sum(p.^2)`.
The dimensional bandwidth factor for a two-dimensional Gaussian rule of
thumb is `n_eff^(-1/6)`, not the univariate `n^(-1/5)` factor.
[SciPy's official weighted KDE documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html).
If using the compact kernel above, specify its bandwidth constant as an
explicit heuristic or perform justified kernel-specific selection; borrowing
a Gaussian constant is not an optimal-bandwidth result. Expose a multiplier
for sensitivity checks. Degenerate coordinates need a stated convention or
clear rejection, not a hidden scale of `1e-6`.

## Weighting contract

The following contract is a package design proposal. Weights express relative
reference membership, not survey design correction or independent replicate
counts. Validate finite, nonnegative entries and a positive total; normalize
stably before arithmetic. A common positive rescaling must change nothing.

For each feature and query value `x`, use the weighted mid-distribution

```text
U(x) = sum_j p_j indicator(x_j<x) + 0.5 sum_j p_j indicator(x_j==x).
```

Zero-weight rows may still receive scores, but must not affect any reference
rank, bandwidth, fitted parameter, or normalization. In particular, adding
zero-weight rows must leave all pre-existing scores unchanged. Unweighted
midranks recover `(tiedrank-.5)/n`. Ties remain ties; the transformed marginal
is not literally continuous-uniform when the original feature is discrete.
Zero-weight queries outside positive-weight support can map to 0 or 1;
bounded KDE can evaluate these directly. A probit-based model needs an
explicit endpoint convention derived from positive reference support, not
from total row count.

For raw typicality scores `s`, use the same weighted mid-distribution
`T_i=sum_j p_j 1{s_j<s_i} + .5 sum_j p_j 1{s_j=s_i}`. This yields `.5` for
an entirely tied reference and is invariant to excluded rows. It is a
relative percentile, not an anomaly probability. Relative to the old
`tiedrank/n`, uniform mid-percentiles shift results by `1/(2n)`; document
this numerical compatibility change. Compute ranks in sorting order with
ties pooled, avoiding a quadratic rank-normalization loop.

Weighting only `fit_typicality` after an unweighted `build_pairs` cannot
repair the marginal geometry. Store transform and reference provenance in
`PairSet`, accept weights during construction, and either reconstruct from
original features or reject incompatible weights at fit time. Never silently
reuse coordinates calculated under a different reference distribution.
Neighbor methods without a defined weighted algorithm should reject
nonuniform membership weights rather than ignore them.

For kernel-contribution leave-one-out scoring, require positive remaining
weight for every scored reference member. State whether ranks and bandwidths
are held fixed at their full-reference values; that is not full pipeline
leave-one-out fitting. If an entity has all reference mass, rejecting the
fit is preferable to declaring it maximally atypical from an empty cloud.

## Narrow parametric copula proposal

Implement `fit_copula(...; family=:gaussian)` only for rank-based pairs with
documented reference weights. The Gaussian family has standard-normal latent
margins; it does not assert normally distributed observed features.
The official `copula` implementation evaluates a multivariate-normal density
divided by the marginal normal densities.
[Official density documentation](https://search.r-project.org/CRAN/refmans/copula/html/Copula.html),
[normal-copula implementation](https://github.com/cran/copula/blob/master/R/normalCopula.R).

For two coordinates, derive the log-density directly:

```text
z1 = Phi^-1(u); z2 = Phi^-1(v)
log c = -0.5*log1p(-rho^2)
        - (rho^2*(z1^2+z2^2)-2*rho*z1*z2)/(2*(1-rho^2)).
```

A feasible parameter estimator is the weighted Pearson correlation of the
normal scores. Call this a normal-score correlation estimate, not maximum
likelihood. Handle constant reference margins explicitly. Require
`abs(rho)<1`; any regularization away from perfect dependence must be
explicit and recorded. Work in log density for comparison/normalization so
underflow or overflow does not create artificial ties. Under independence
`rho=0`, log density is exactly zero everywhere: this is a valuable oracle
that Mahalanobis depth cannot satisfy.

This family cannot represent every nonlinear, asymmetric or multimodal
dependence pattern. It adds a clear model-based comparator; it is not a
general replacement for LOF or bounded KDE. Fitted-density scores remain
descriptive and can be optimistic when evaluated on the fitting collection.

## Neighbor alternatives screened

The original LOF neighborhood contains every entity at or within the kth
distance. Repair this definition and zero-distance cases before extending
the method menu.
[Breunig et al., original paper](https://cs.fit.edu/~pkc/ml/related/breunig00sigmod.pdf).

INFLO adds reverse nearest neighbors to its influence space. This changes
what local context means, and requires its own efficiency and tie treatment;
it does not remove rank-boundary or reference-weight problems.
[Jin et al., original paper](https://hanj.cs.illinois.edu/pdf/pakdd06_outlier.pdf).

LoOP presents a probability-oriented outlier score. The historical PDF URL
redirected during inspection, so only the author record and indexed abstract
were screened. Its probabilistic interpretation should not be adopted merely
because its output lies in `[0,1]`. Neither neighbor alternative was selected
for implementation in this revision.
[Author publication record](https://brava.dbs.ifi.lmu.de/publications/891).

## Implementation and verification record

The implementation is in `src/density.jl`, with focused tests in
`test/density.jl`. `_density_scores` returns scores and support, bandwidth,
constant-coordinate, effective-size and leave-one-out metadata. It exposes
log-density output, uses sorted reductions for permutation reproducibility,
and rejects a reference with fewer than two positive-weight entities. Its
bandwidth heuristic uses weighted population variance and the dimensional
effective-size factor. For the compact kernel, a `sqrt(5)` factor matches
kernel variance to a unit Gaussian; bandwidths are capped at one. Explicit
bandwidths are absolute coordinate units, scalar or a two-coordinate tuple.

Constant reference coordinates are dropped from the kernel product. Queries
off that constant support receive zero density; two constant coordinates give
one on their joint support. These are reduced-support scoring conventions,
not a two-dimensional Lebesgue density for degenerate data.

`_copula_scores` returns log densities and fitted parameters. Gaussian
normal-score correlations are shrunk only when the caller supplies positive
`regularization`; a numerically singular unregularized fit raises an error.
Endpoint queries raise an error by default because a corner copula-density
limit can depend on the path of approach. The caller may specify
`query_clip` strictly between zero and one half to clip zero-weight endpoint
queries only; this setting is recorded in metadata. Positive-reference
endpoints always raise an error. Beta scores use the full equal-weight empirical beta mixture and
are explicitly marked in-sample. Neither is claimed to use full-pipeline
leave-one-out fitting.

The focused tests cover the independent-grid boundary check, analytic
Gaussian-kernel values, actual reflected-kernel quadrature, raw/log
consistency, weight rescaling, zero-weight extension, row permutations,
constant support, analytic Gaussian-copula identities, the independent
Gaussian fit, singular-fit handling, and the exact two-observation beta
density. Integration tests additionally need to verify the public API's
weighted transformation and result normalization.

The first focused standalone run passed 128 checks under Julia 1.13 with
compiled modules disabled. The default development-runtime precompilation
failed in Julia's Markdown/Test dependency stack, so verification was moved
to the installed Julia 1.10 runtime. Endpoint-policy and metadata tests were
added after that initial run; the final run is recorded in the parent
implementation report. The final standalone density/copula suite passed
**133/133 checks on Julia 1.10.12** in 3.2 seconds after compilation, using
the temporary review depot and the repository project. This includes the
explicit endpoint policy and raw/log-density metadata checks.

## Required numerical oracles and limitations

1. Permute entity rows, reverse ties, and permute feature columns; restored
   pair scores and summaries must agree within declared tolerances.
2. Append zero-weight points both inside and far outside reference support;
   existing transformed coordinates and normalized scores must not change.
3. Rescale all weights by small and large finite constants; test normalization
   without introducing epsilon-dependent estimates.
4. Integrate each reflected one-dimensional kernel numerically for centers
   at 0, near 0, .5, near 1, and 1, across allowed bandwidths. Verify unit mass.
5. On a dense independent uniform grid, inspect corner/interior KDE density
   ratios before ranking. Avoid claiming perfect finite-sample equality,
   particularly for leave-one-out evaluation.
6. For a Gaussian copula, test `rho=0` against density 1; coordinate exchange;
   simultaneous reflection; `c_rho(u,v)=c_-rho(u,1-v)`; and selected interior
   density values against an independent implementation or analytic constants.
7. Reject or deliberately handle singleton effective reference support,
   constant features, NaN/Inf, invalid bandwidths, and singular dependence.
8. Use seeded dependence-violation fixtures for every implemented estimator,
   alongside null fixtures. Compare sensitivity to `k`, bandwidth and feature
   subsets rather than asserting every method must induce one ordering.

These checks validate stated numerical behavior, not construct validity.
Population bivariate densities do not identify every higher-order interaction;
this is not a theorem that every aggregation of pair scores is powerless. Scores
sharing features are dependent; IQR across pairs is not a standard error.
Any subsequent bootstrap must resample entities or appropriate study-design
clusters and repeat preprocessing, rather than resampling pair columns as
independent observations. No inferential confidence level or human judgment
agreement is claimed by these proposals.
