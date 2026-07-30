# bml 0.11.0

## Compact posterior storage

* `monitor` is now capability-based and defaults to `"summary"`, which retains
  estimates and convergence diagnostics without retaining posterior draws.
* Draws requested through `"parameters"`, `"random_effects"`, `"weights"`,
  `"fitted"`, `"predictive"`, or `"log_lik"` are stored once in a compact
  `posterior::draws_array`. `"full"` selects every family-supported capability.
* `monitor = "jags"` is the explicit escape hatch for retaining the complete
  raw R2jags object. Logical `monitor` values are no longer accepted.
* Summary-only fits retain rank-normalized R-hat, bulk/tail ESS, and Monte Carlo
  standard errors for every reported parameter. `summary()`, `tidy()`,
  `glance()`, and `bmlCompare()` expose compact convergence information.
* Post-estimation methods now report the exact posterior capability needed when
  a fit did not retain the required draws.

# bml 0.10.0

This release redesigns the model-specification grammar around the micro-macro
framework (weights, aggregation functions, effects) and aligns the package API
with brms. **All breaking changes are clean breaks — old spellings error with a
migration message; there are no deprecation shims.**

## The aggregation function `fn()` (emergent terms)

* Every `mm()` block now selects an **aggregation function** `fn()` — the `f` in
  the framework's `theta^{micro,f}(M_it)`. `fn("sum")` (default) is the additive
  weighted mean; new **emergent** types read properties of the whole member set:
  - `fn("var", moment = p)` — weighted central moments (distributional shape)
  - `fn("hhi")`, `fn("effn")`, `fn("entropy")` — concentration of the weights
  - `fn("threshold", c =, kappa =)` — smooth critical-mass terms; the cutpoint
    can be **estimated**: `c = est()`
  - `fn("smax", kappa =)` — smooth max/min; `kappa = est()` makes the
    aggregation function itself the estimand (min <-> mean <-> max)
  - `fn("gmean", p = est())` — power/CES mean (substitutes vs complements)
  - `fn("cov")` — weighted covariance of two member attributes
* **Expression DSL:** `fn(~ E((x - E(x))^2))` — user-written reductions with the
  expectation operator `E()` over the normalized member weights; whitelisted
  operations transpile to JAGS. Free symbols become estimable shape parameters.
* Blocks stack: combine mean + variance + concentration features in one model.

## Weights: `fn()` -> `w()`

* Weights moved to their own component: `w(~ 1/n)`, `w(~ importance,
  scale = TRUE)`, `w(~ ilogit(b0 + b1*q))`. The old `fn(w ~ ...)` spelling and
  the `c =` normalization flag were **removed** (`c` -> `scale`, matching brms).
* **One parameter rule:** any symbol in a `w()` or `fn()` formula that is not a
  data column is a free parameter; the build messages
  `Estimating parameters: ...` so typos surface immediately.

## Named blocks and interactions

* Every block has a (auto-generated or user-set) `name =`; single-feature blocks
  can be referenced in main-formula interactions: cross-level
  (`Ax:education`) and block x block (`Ax:Vx`). Bare references error.

## Effects grammar: `re()` / `fe()`

* `RE = re(1 + x)` — random intercepts + slopes (lme4-style LHS; grouping
  implied by `id()`). Default **`cor = FALSE`** (independent; a deliberate
  divergence from brms's correlated default); `cor = TRUE` opts into the
  correlated 2x2 with an LKJ-style prior.
* `FE = fe(1 + x, showFE =)` — the no-pooling peer (unit dummies + unit-specific
  slopes; first unit is the reference).
* `hm()` carries **structure only**: `type =` and `vars =` were removed. Fixed
  effects of unit-level covariates go in the main formula
  (`Y ~ gdp + hm(id(cid), RE = re(1 + gdp))`), retiring `b.hm`.

## brms-aligned call syntax

* `bml(formula, data, family = gaussian(), prior = ...)` — `data` is second.
* **Family functions:** `gaussian()`, `bernoulli()`, `weibull()`,
  `cox(intervals =)` (absorbing `cox_intervals`); strings still accepted.
* **MCMC arguments renamed:** `iter`, `warmup`, `thin`, `chains`, `cores`
  (replacing `n.iter`, `n.burnin`, `n.thin`, `n.chains`, `parallel`).
* **`prior()` system:** class-targeted priors on the natural (SD) scale,
  translated to JAGS automatically — `prior(normal(0, 5), class = "b")`,
  classes `Intercept`, `b`, `sd`, `cor`, `sigma`, `shape`, `w`, `fn`;
  `get_prior()` lists every settable prior. Raw JAGS strings remain an escape
  hatch. The old `priors =` list argument was removed.

## Parameter naming

* One user-facing coefficient family **`b`**, labeled by term
  (`(Intercept)`, `majority`, `A_x`, `A_x:education`); internal arrays `b`,
  `b.fn.<k>` retire `b.mm`/`b.hm`. Weight parameters are `b.w.<k>` (class
  `"w"`); `fn` shape parameters are `fn.<name>.<k>` (class `"fn"`) and are
  reported in their own "Function parameters" section, never in the coefficient
  table. RE variances print as `sd(Intercept)`, `sd(x)`, `cor(Intercept,x)`.

## Post-estimation: the brms / posterior vocabulary

* `as_draws()` / `as_draws_df()` / `as_draws_matrix()` / `as_draws_array()`
* `fixef()`, `ranef()`, `coef()`
* `fitted()`, `predict()`, `posterior_predict()` (newdata not yet supported)
* `log_lik()`, `loo()`, `waic()` (gaussian/bernoulli; pointwise `log_lik` is
  monitored in the generated JAGS model)
* `pp_check()` (via bayesplot, Suggests), `conditional_effects()`
* `make_jagscode()` / `make_jagsdata()` (analogues of `make_stancode()`)

## Autoregressive random effects: `ar` moved into `re()`

* The random-walk option is part of the effects grammar now — it is a property
  of the random effect, not of the block. `mm(ar = TRUE)` / `hm(ar = TRUE)`
  were **removed** (migration error): write `RE = re(1, ar = TRUE)`.
* **New: time-indexed walks.** `RE = re(1, ar = year)` orders each member's
  participations by a numeric time variable and scales the step variance by
  the normalized time gaps, `u_t ~ N(u_{t-1}, sigma^2 * gap_t / mean(gap))` —
  long gaps allow more drift than short ones; equal spacing reproduces the
  `ar = TRUE` model. Requires an intercept-only `re(1)`; time values must be
  unique within each member/unit.

## Temporal helpers

* `bml_delta()`, `bml_realloc()`, `bml_turnover()` — aggregation over
  membership spells: the matched (shift-share style) decomposition of a
  composition change into reallocation + attribute change + entry/exit, with
  the adding-up identity holding exactly. Observed-panel case; the in-model
  spell path (estimated weights, reallocation toward high-RE members) is
  future work.

## Vignettes

* The examples vignette gained four new examples on the coalition data:
  mean vs. spread (`fn("var")`, stacked feature blocks), estimating the
  aggregation function itself (`fn("smax", kappa = est())`), cross-level
  interactions via named blocks (`Afin:majority`), and heterogeneous member
  effects (`vars(x + x:z)` + `re(1 + x)` in one block). The friendship example
  now adjudicates linear-in-means vs. best-friend with `loo_compare()` and
  shows `pp_check()`; the Boston example uses `prior(..., class = "w")`.

## Bug fixes

* **Canonical member-row ordering.** All member-level structures (design
  matrices, weights, slope data, and the member-id index vectors) now share one
  ordering, sorted by (group, member). Previously the member-id index used a
  different ordering than the design/weight matrices, so models with member
  random effects could associate rows with the wrong member's effect when the
  input data were not already sorted by member within group. Results of
  affected models change (even with the same seed); the refit examples
  vignette reflects this — most notably, the friendship example's model
  comparison now favors linear-in-means over the best-friend specification.

## Removed (clean breaks)

* `fn(w ~ ...)` weights spelling; `c =` normalization flag
* `hm(type =)`, `hm(vars =)`, `hm(name =)` (now `labels =`), `hm(showFE =)`
* `mm(ar =)` / `hm(ar =)` — moved into `re(1, ar = ...)`
* `bml(n.iter =, n.burnin =, n.thin =, n.chains =, parallel =, priors =,
  cox_intervals =)`; `data` as last argument
* `family = "Gaussian"`-style capitalized strings still resolve, but the
  canonical spellings are lowercase / family functions

---

# bml 0.9.0

## Major Changes

* **Package renamed from `rmm` to `bml`** (Bayesian Multiple-Membership Multilevel Models)

* **New syntax for weight functions:** The `ar` parameter has been moved from the `fn()` specification to the `mm()` block level for clearer API
  - Old: `fn(w ~ 1/n, c = TRUE, ar = FALSE)`
  - New: `fn(w ~ 1/n, c = TRUE)` with `ar = FALSE` at the `mm()` level

* **Support for multiple mmid groups:** The package now supports models with multiple membership identifiers, allowing more complex membership structures

* **Enhanced documentation:** Comprehensive documentation added for the `coalgov` dataset including:
  - Detailed variable descriptions organized by level (identifiers, government-level, country-level, party-level)
  - Statistical summaries for all variables
  - Clear explanation of multiple-membership structure
  - Updated references and examples

## New Features

* **Flexible weight function parameterization:** Enhanced support for parameterizing weight functions with covariates and group-specific structures

* **Per-group random effects:** Random effects can now be specified separately for different mmid groups

* **Improved JAGS code generation:** Optimized model string generation for better performance with complex multiple-membership structures

## Breaking Changes

* **`ar` parameter moved:** Existing code using `fn(w ~ ..., ar = TRUE)` must be updated to place `ar` in the `mm()` block instead

* **Dataset changes:**
  - Removed `schoolnets` dataset (including `nodedat` and `edgedat` objects)
  - Updated `coalgov` dataset with enhanced documentation and additional variables

## Bug Fixes

* Fixed issues with weight function constraints when using multiple `mm()` blocks

* Improved handling of group-level indices in JAGS variable creation

## Documentation

* Updated vignette examples to use new syntax
* Added comprehensive FAQ section
* Improved installation instructions with CRAN and GitHub options
* Updated references to reflect 2026 publication

---

# rmm 0.2.0

* Major syntax updates for weight function specification
* Improved prior specification interface
* Updated model string generation

# rmm 0.1.1

* Removed posterior predictive p-values and fixed variance at second level for Weibull AFT regression

# rmm 0.1.0

* Initial release
* Added a `NEWS.md` file to track changes to the package
