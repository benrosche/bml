# Frequently Asked Questions

## 1. What’s the difference between a multiple-membership model and a conventional multilevel model?

In a conventional hierarchical multilevel model, each lower-level unit
belongs to exactly one higher-level unit (e.g., students nested in
schools). In a *multiple-membership model (MMMM)*, lower-level units can
belong to multiple higher-level units simultaneously, or alternatively,
higher-level outcomes can be influenced by multiple lower-level units.

For example:

- **Conventional MLM:** Students (level 1) nested in schools (level 2) —
  each student attends exactly one school.
- **MMMM:** Coalition governments (level 2) composed of multiple
  political parties (level 1) — each party can participate in multiple
  governments over time, and each government comprises multiple parties.

The MMMM accounts for this complex membership structure by weighting the
contributions of each lower-level unit, allowing researchers to model
how multiple units jointly shape higher-level outcomes.

## 2. What’s the difference between a conventional MMMM and the extended MMMM implemented in `bml`?

The *conventional MMMM* (as implemented in MLwiN, brms, or other
software) uses *fixed, pre-specified weights* to aggregate lower-level
effects. For instance, you might use equal weights (1/n) or weights
based on time spent in each context.

The *extended MMMM* in `bml` allows you to:

1.  **Parameterize the weight function:** Instead of fixing weights, you
    can specify a functional form for weights (e.g., `w ~ 1/n^exp(b*X)`)
    and estimate the parameters that determine how lower-level units are
    aggregated.

2.  **Test alternative aggregation mechanisms:** Compare different
    weighting schemes (equal weights, proportional weights, functions of
    covariates) to determine which best fits the data.

3.  **Endogenize weight matrices:** Rather than imposing spatial or
    network weights externally, let the data determine connection
    strengths as functions of covariates.

This flexibility enables researchers to explicitly model the
*micro-to-macro link* — how lower-level characteristics aggregate to
produce higher-level outcomes.

## 3. When should I use `bml` instead of other multilevel modeling packages?

Use `bml` when:

- **You have multiple-membership structures:** Higher-level outcomes
  depend on multiple lower-level units (e.g., coalitions composed of
  parties, teams composed of individuals, neighborhoods influenced by
  surrounding areas).

- **You need flexible weight functions:** Your theory suggests weights
  should depend on covariates, group size, or other features, not be
  fixed in advance.

- **You’re studying aggregation processes / micro-to-macro
  relationships:** Your research question focuses on how lower-level
  units jointly shape higher-level outcomes, rather than how
  higher-level contexts shape lower-level outcomes. Rather than assuming
  a fixed aggregation (e.g., simple average), you want to test and
  estimate how lower-level effects combine.

For standard hierarchical models without multiple membership, packages
like `lme4`, `brms`, or `MCMCglmm` will be more efficient. For
conventional MMMMs with fixed weights, `brms` or MLwiN are excellent
alternatives.

## 4. What outcome types and distributions does `bml` support?

`bml` supports a variety of outcome distributions commonly used in
social science research:

- **Continuous outcomes:** Gaussian (normal) regression
- **Binary outcomes:** Logit (logistic) regression
- **Survival/duration outcomes:**
  - Cox proportional hazards model
  - Weibull accelerated failure time (AFT) model

You can also specify hierarchical random effects
([`hm()`](https://benrosche.github.io/bml/reference/hm.md) blocks) in
addition to multiple-membership structures, allowing for hierarchically
nested and cross-classified designs.

## 5. How do I specify the weights, and what are the `scale` and `ar` options?

The weights are specified in the
[`w()`](https://benrosche.github.io/bml/reference/w.md) container within
an [`mm()`](https://benrosche.github.io/bml/reference/mm.md) block (the
aggregation function lives in
[`fn()`](https://benrosche.github.io/bml/reference/fn.md), which
defaults to the additive `fn("sum")`):

``` r

mm(
  id   = id(member_id, group_id),
  vars = vars(X),
  w    = w(~ 1/n, scale = TRUE),
  fn   = fn("sum"),
  RE   = TRUE
)
```

**Weight function components:**

- **`w(~ ...)`**: Specifies the functional form for weights (one-sided
  formula). You can use:
  - Group size: `n` (number of members in each group)
  - Covariates: Any variable in your data
  - Mathematical functions: [`exp()`](https://rdrr.io/r/base/Log.html),
    [`log()`](https://rdrr.io/r/base/Log.html),
    [`sqrt()`](https://rdrr.io/r/base/MathFun.html), etc.
  - Examples: `w(~ 1/n)` (equal weights), `w(~ tenure)` (proportional to
    tenure), `w(~ 1/n^exp(b*similarity))` (similarity-weighted)
  - **One parameter rule:** any symbol that is not a data column (like
    `b` above) becomes a free parameter with a default prior; the build
    messages `Estimating parameters: ...`
- **`scale` parameter:** Controls weight normalization
  - `scale = TRUE` (default): Weights are normalized to sum to 1 within
    each group (row-standardization); this makes them a probability
    measure over the group’s members
  - `scale = FALSE`: Weights are not normalized (useful when aggregating
    sums rather than averages)
- **Autoregressive random effects** live inside
  [`re()`](https://benrosche.github.io/bml/reference/re.md) (they are a
  property of the random effect):
  - `RE = re(1)` (default): independent random effects for each member
  - `RE = re(1, ar = TRUE)`: member effects evolve as a random walk
    across repeated group participations
  - `RE = re(1, ar = year)`: random walk over a numeric time variable,
    with the step variance scaled by the (normalized) time gaps — long
    gaps allow more drift than short ones

**Important:** Only one
[`mm()`](https://benrosche.github.io/bml/reference/mm.md) block per
member-id can carry random effects.

## 6. How do I fix parameters to known values?

There are two ways to fix parameters in `bml`, depending on where in the
model they appear.

### Main equation and `hm()` blocks: Using `fix()`

Use the [`fix()`](https://benrosche.github.io/bml/reference/fix.md)
helper inside
[`vars()`](https://benrosche.github.io/bml/reference/vars.md) to hold a
coefficient at a specified value rather than estimating it. This works
both in the main equation and in
[`hm()`](https://benrosche.github.io/bml/reference/hm.md) blocks.

**Main equation:**

``` r

# Fix the coefficient of 'exposure' to 1 (i.e., use as an offset)
m <- bml(
  y ~ 1 + fix(exposure, 1) + x1 + x2,
  data = dat,
  family = gaussian()
)
```

This is equivalent to adding `exposure` as a predictor but constraining
its coefficient to 1 instead of estimating it. This is useful for
offsets or when prior theory dictates a specific coefficient value.

**Unit-level covariates:** with the
[`re()`](https://benrosche.github.io/bml/reference/re.md)/[`fe()`](https://benrosche.github.io/bml/reference/fe.md)
grammar, [`hm()`](https://benrosche.github.io/bml/reference/hm.md)
carries only the effect structure; unit-level covariates (and
[`fix()`](https://benrosche.github.io/bml/reference/fix.md) offsets of
them) go in the main formula:

``` r

m <- bml(
  y ~ 1 + x1 + fix(investiture, 0.5) + gdp +
    mm(id = id(pid, gid), w = w(~ 1/n), fn = fn("sum"), RE = TRUE) +
    hm(id = id(cid), FE = fe(1)),
  data   = coalgov,
  family = gaussian()
)
```

**In [`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks:**

``` r

m <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(
      id   = id(pid, gid),
      vars = vars(fix(cohesion, 1) + finance),
      w    = w(~ 1/n, scale = TRUE),
      fn   = fn("sum"),
      RE   = TRUE
    ),
  data   = coalgov,
  family = weibull()
)
```

Here, the coefficient of `cohesion` is fixed at 1 while `finance` is
freely estimated.

### Weights `w()`: Omitting parameters

In the weight formula, you fix weights by simply not including any free
symbols. When every symbol is a data column (or `n`), the weights are a
known, fixed transformation and are pre-computed in R.

``` r

# Fixed equal weights (no parameters to estimate)
w(~ 1/n, scale = TRUE)

# Fixed weights proportional to seat share (no parameters to estimate)
w(~ pseat, scale = TRUE)
```

Compare this with weight formulas that include free parameters:

``` r

# Estimable: b is a free symbol, so it is estimated from data
w(~ 1/n^exp(b * x), scale = TRUE)

# Estimable: b0 is estimated from data
w(~ 1 / (1 + (n - 1) * exp(-b0)), scale = FALSE)
```

The key distinction is the *one parameter rule*: any symbol that is not
a data column or `n` is a free parameter with a default prior. The build
always messages the detected parameters (`Estimating parameters: ...`),
so a misspelled column name surfaces immediately.

## 7. I get “Error in node w.1\[…\]: Invalid parent values” — what does this mean?

This error occurs when JAGS cannot evaluate the weight function for a
particular observation because the computed weight is numerically
invalid (e.g., `NaN`, `Inf`, or a value outside the domain of a
downstream function). It most commonly arises with **parameterized
weight functions** during the initialization phase.

**Why it happens:** Weight function parameters (`b0`, `b1`, …) are given
vague priors by default (`dnorm(0, 0.0001)` in JAGS precision, which
corresponds to SD = 100). When JAGS initializes the MCMC chains, it may
draw extreme starting values (e.g., `b0 = 80`). For weight functions
involving nonlinear transformations like `ilogit()` or
[`exp()`](https://rdrr.io/r/base/Log.html), extreme parameter values can
cause numerical issues downstream — even if the weight function itself
is mathematically well-defined, the resulting weights may produce
overflow in the likelihood (e.g., `exp(-mu * shape)` in the Weibull
model).

**Built-in safeguard:** `bml` initializes all weight parameters at 0 by
default. This ensures numerically stable starting values (e.g.,
`ilogit(0) = 0.5`). However, if the error persists or occurs during
sampling (not just initialization), consider the following steps.

**How to fix it:**

1.  **Use more informative priors.** Narrow the prior spread for weight
    parameters so that the sampler stays in a numerically stable region:

    ``` r

    m <- bml(
      ...,
      prior = prior(normal(0, 1), class = "w")   # SD = 1 instead of 100
    )
    ```

    This is especially important for parameters inside `ilogit()`,
    [`exp()`](https://rdrr.io/r/base/Log.html), or other functions that
    saturate or explode at extreme inputs.

2.  **Ensure the weight function is bounded.** Unbounded weight
    functions can produce extreme values that destabilize the
    likelihood. Common strategies:

    - Use `ilogit()` to bound weights between 0 and 1:
      `w(~ ilogit(b0 + b1 * x), scale = TRUE)`
    - Use [`exp()`](https://rdrr.io/r/base/Log.html) carefully — it
      grows rapidly, so pair it with `scale = TRUE` (normalization) or
      wrap the argument: `w(~ exp(b1 * x), scale = TRUE)` where `x` is
      standardized

3.  **Standardize covariates in the weight function.** If a weight
    variable has a large range (e.g., income in thousands), the product
    `b1 * x` can easily overflow. Standardize such variables before
    including them in
    [`w()`](https://benrosche.github.io/bml/reference/w.md).

4.  **Supply explicit initial values.** If the default initialization at
    0 doesn’t work for your model, provide custom starting values:

    ``` r

    m <- bml(
      ...,
      inits = list("b.w.1" = c(0.5, -0.1))
    )
    ```

5.  **Re-run the model.** Since the error can be seed-dependent
    (different chains draw different initial values), simply re-running
    may resolve it. However, this indicates a fragile parameterization —
    consider steps 1–3 for a robust solution.

## 8. I’m coming from brms — what maps onto what?

`bml`’s member random effects aggregated by weights are exactly brms’s
weighted multiple-membership term, and the `bml` API deliberately
follows brms where the semantics coincide:

| brms | bml | note |
|----|----|----|
| `brm(formula, data, family = gaussian())` | `bml(formula, data, family = gaussian())` | same positions; family functions or strings |
| `iter`, `warmup`, `thin`, `chains`, `cores`, `seed` | same | same meanings |
| `prior(normal(0, 5), class = "b")` | same | plus bml-specific classes `"w"` (weight parameters) and `"fn"` (aggregation shape parameters) |
| [`get_prior()`](https://benrosche.github.io/bml/reference/get_prior.md) | [`get_prior()`](https://benrosche.github.io/bml/reference/get_prior.md) | lists every settable prior |
| `(1 + x \| mm(g1, g2, weights = w, scale = TRUE))` | `mm(id = id(k, i), w = w(~ q, scale = TRUE), fn = fn("sum"), RE = re(1 + x))` | see the divergences below |
| `mmc(x1, x2)` | not needed | `bml`’s long format gives every member row its own covariate value |
| [`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html), [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html), [`coef()`](https://rdrr.io/r/stats/coef.html), [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), [`predict()`](https://rdrr.io/r/stats/predict.html), [`posterior_predict()`](https://benrosche.github.io/bml/reference/posterior_predict.md), [`log_lik()`](https://benrosche.github.io/bml/reference/log_lik.md), [`loo()`](https://mc-stan.org/loo/reference/loo.html), [`waic()`](https://mc-stan.org/loo/reference/waic.html), [`pp_check()`](https://benrosche.github.io/bml/reference/pp_check.md), [`conditional_effects()`](https://benrosche.github.io/bml/reference/conditional_effects.md), [`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html) | same | `newdata` prediction is not supported yet |
| `make_stancode()` / `make_standata()` | [`make_jagscode()`](https://benrosche.github.io/bml/reference/make_jagscode.md) / [`make_jagsdata()`](https://benrosche.github.io/bml/reference/make_jagscode.md) |  |
| `sd(Intercept)`, `sd(x)`, `cor(Intercept,x)` output labels | same |  |

**Deliberate divergences:**

- **`id(member, group)` is long format** (one row per membership, any
  number of members per group), not brms’s wide `mm(g1, g2)` (fixed
  number of membership columns). Long format is strictly more general
  for variable set sizes.
- **[`re()`](https://benrosche.github.io/bml/reference/re.md) defaults
  to `cor = FALSE`** (independent intercept and slope), where brms’s `|`
  is correlated. The correlated draw is a real JAGS cost, and in
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks the
  member effects are only read through the weighted aggregate, so the
  correlation is rarely of interest. Opt in with
  `re(1 + x, cor = TRUE)`.
- **[`bml::mm()`](https://benrosche.github.io/bml/reference/mm.md) is
  more than brms’s
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md):** brms’s
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md) is only the
  grouping constructor inside `( | )`; `bml`’s
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md) bundles the
  grouping, the weights
  [`w()`](https://benrosche.github.io/bml/reference/w.md), the
  aggregation function
  [`fn()`](https://benrosche.github.io/bml/reference/fn.md) (including
  emergent features like `fn("var")`), and the member effects
  `RE = re(...)`.
