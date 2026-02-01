# Bayesian Multiple-Membership Multilevel Models with Parameterizable Weight Functions Using JAGS

The **bml** package provides a user-friendly interface for fitting
Bayesian multiple-membership multilevel models with parameterizable
weight functions via JAGS.

JAGS must be installed separately:
<https://sourceforge.net/projects/mcmc-jags/>.

## Usage

``` r
bml(
  formula,
  family = "Gaussian",
  priors = NULL,
  inits = NULL,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = max(1, floor((n.iter - n.burnin)/1000)),
  n.chains = 3,
  seed = NULL,
  run = TRUE,
  parallel = FALSE,
  monitor = TRUE,
  modelfile = FALSE,
  cox_intervals = NULL,
  data = NULL
)
```

## Arguments

- formula:

  A symbolic model formula. See 'Formula Components' section for
  details. The general structure is:
  `outcome ~ 1 + predictors + mm(...) + hm(...)`. For survival models,
  use `Surv(time, event)` on the left-hand side.

- family:

  Character string specifying the outcome distribution and link
  function. Options:

  - `"Gaussian"`: Normal distribution with identity link (continuous
    outcomes)

  - `"Binomial"`: Binomial distribution with logit link (binary
    outcomes)

  - `"Weibull"`: Weibull survival model (requires `Surv(time, event)`
    outcome)

  - `"Cox"`: Cox proportional hazards model (requires
    `Surv(time, event)` outcome)

- priors:

  Named list or character vector of JAGS prior specifications. Parameter
  names follow the pattern: `b.mm.k` (mm block k coefficients), `b.w.k`
  (weight function parameters), `tau.mm` (mm random effect precision),
  etc. Example:
  `list("b.mm.1 ~ dnorm(0, 0.01)", "tau.mm ~ dgamma(2, 0.1)")`. Default
  priors are weakly informative.

- inits:

  List of initial values for MCMC chains. Applied to all chains. If
  `NULL`, JAGS generates initial values automatically.

- n.iter:

  Total number of MCMC iterations per chain. Default: 1000. Increase for
  better convergence (e.g., 10000-50000 for production models).

- n.burnin:

  Number of burn-in iterations to discard at the start of each chain.
  Default: 500. Should be sufficient for chains to reach stationarity.

- n.thin:

  Thinning rate: save every k-th iteration to reduce autocorrelation.
  Default: `max(1, floor((n.iter - n.burnin) / 1000))` (targets ~1000
  samples). Increase if posterior samples show high autocorrelation.

- n.chains:

  Number of MCMC chains. Default: 3. Use 3-4 chains to assess
  convergence via Gelman-Rubin diagnostics.

- seed:

  Integer random seed for reproducibility. If `NULL`, results will vary
  across runs.

- run:

  Logical; if `TRUE` (default), JAGS is executed and the model is
  fitted. If `FALSE`, returns the model specification without fitting
  (useful for inspecting generated JAGS code or data structures).

- parallel:

  Logical; if `TRUE`, run MCMC chains in parallel using multiple cores.
  Requires parallel backend setup. Default: `FALSE`.

- monitor:

  Logical; if `TRUE`, store full MCMC chains and additional outputs for
  diagnostic plots. Required for
  [`monetPlot`](https://benrosche.github.io/bml/reference/monetPlot.md)
  and
  [`mcmcDiag`](https://benrosche.github.io/bml/reference/mcmcDiag.md).
  Default: `TRUE`.

- modelfile:

  Logical or character path:

  - `FALSE` (default): JAGS code generated internally

  - `TRUE`: Save generated JAGS code to `modelstring.txt` in working
    directory

  - Character path: Read JAGS code from specified file instead of
    generating

- cox_intervals:

  For Cox models only: controls baseline hazard flexibility and
  computational efficiency.

  - `NULL` (default): Non-parametric baseline hazard using all unique
    event times (maximum flexibility, slower for large datasets)

  - Integer k: Piecewise constant baseline hazard with k intervals
    (faster, suitable for datasets with many unique event times).
    Recommended: k = 10-20 for most applications.

- data:

  Data frame in member-level (long) format where each row represents a
  member-level observation. Must contain all variables referenced in the
  formula, including identifiers specified in
  [`id()`](https://benrosche.github.io/bml/reference/id.md).

## Value

A list of class `"bml"` containing:

- `reg.table`: Data frame of posterior summaries (means, SDs, credible
  intervals) for main parameters. Access columns via `$Parameter`,
  `$mean`, `$sd`, `$lb`, `$ub`.

- `jags.out`: Full JAGS output object (if `monitor = TRUE`)

- `modelstring`: Generated JAGS model code

- `jags.data`: Data passed to JAGS

- `jags.inits`: Initial values used

- `formula`: Original formula

- `family`: Model family

- Additional components for random effects, weights, predictions (when
  applicable)

## Details

In addition to hierarchical and cross-classified multilevel models, the
**bml** package allows users to fit Bayesian multiple-membership models.
Unlike tools such as `brms` or
[MLwiN](https://www.bristol.ac.uk/cmm/software/mlwin/), **bml** lets
users specify and estimate models in which membership weights are
parameterized through flexible formula syntax. This enables a more
nuanced examination of how effects from member-level units aggregate to
group level (the micro-macro link).

The package automatically generates JAGS code to fit the model and
processes the output to facilitate interpretation of model parameters
and diagnostics.

The package and modeling framework are introduced in: Rosche, B. (2026).
*A Multilevel Model for Coalition Governments: Uncovering Party-Level
Dependencies Within and Between Governments*. *Political Analysis*.

For accessible introductions to multiple-membership models, see Fielding
and Goldstein (2006) and Beretvas (2010). Advanced treatments include
Goldstein (2011, Ch. 13), Rasbash and Browne (2001, 2008), Browne et al.
(2001), and Leckie (2013).

## Formula Components

- **Outcome (Y):** The dependent variable. For survival models, use
  `Surv(time, event)`.

- **Intercept:** Follows standard R formula conventions (like
  [`lm()`](https://rdrr.io/r/stats/lm.html)):

  - `y ~ x`: Includes intercept by default

  - `y ~ 1 + x`: Explicitly includes intercept (same as default)

  - `y ~ 0 + x` or `y ~ -1 + x`: Excludes intercept

- **Main-level predictors (X.main):** Variables defined at the main
  (group) level, separated by `+`.

- **HM-level predictors (X.hm):** Variables defined at the nesting
  level, separated by `+`.

- **Multiple membership object
  ([`mm()`](https://benrosche.github.io/bml/reference/mm.md)):** Defines
  how member-level units are associated with group-level constructs
  using a user-specified weighting function. Multiple
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md) objects can
  be specified with different weight functions.

- **Hierarchical membership
  ([`hm()`](https://benrosche.github.io/bml/reference/hm.md)):**
  Specifies nesting of main-level units within higher-level entities.
  Cross-classified structures can be modeled by including multiple
  [`hm()`](https://benrosche.github.io/bml/reference/hm.md) objects.

**Formula Features:** The main formula and
[`vars()`](https://benrosche.github.io/bml/reference/vars.md)
specifications support standard R formula syntax:

- **Interactions:** Use `*` for main effects plus interaction, or `:`
  for interaction only. Example: `y ~ a * b` expands to
  `y ~ a + b + a:b`.

- **Transformations:** Use [`I()`](https://rdrr.io/r/base/AsIs.html) for
  arithmetic operations. Example: `y ~ I(x^2)` or `y ~ I(a + b)`.

These features work in:

- **Main formula:** `y ~ 1 + a * b + I(x^2)`

- **mm() vars:** `vars(a * b)` or `vars(I(x^2))`

- **hm() vars:** `vars(a:b)` or `vars(I(log(x)))`

**Note on weight functions:** The
[`fn()`](https://benrosche.github.io/bml/reference/fn.md) weight
function in [`mm()`](https://benrosche.github.io/bml/reference/mm.md)
does NOT support interactions or
[`I()`](https://rdrr.io/r/base/AsIs.html) transformations. Users must
pre-create any needed transformed variables in their data before using
them in weight functions. For example, instead of `fn(w ~ b1 * x^2)`,
first create `data$x_sq <- data$x^2` and use `fn(w ~ b1 * x_sq)`.

**Note on intercepts:** Intercept syntax (`1`, `0`, `-1`) only applies
to the main formula. Numeric literals in
[`vars()`](https://benrosche.github.io/bml/reference/vars.md) are
ignored (e.g., `vars(1 + x)` is equivalent to `vars(x)`).

## Multiple Membership Object [`mm()`](https://benrosche.github.io/bml/reference/mm.md)


    mm(
      id   = id(mmid, mainid),
      vars = vars(X.mm),
      fn   = fn(w ~ 1/n, c = TRUE),
      RE   = TRUE,
      ar   = FALSE
    )

**Components:**

- `id(mmid, mainid)`: Specifies identifiers linking each member-level
  unit (`mmid`) to its corresponding group-level entities (`mainid`).

- `vars(X.mm)`: Specifies member-level covariates aggregated across
  memberships. Use `+` to include multiple variables. Supports
  interactions (`*`, `:`) and transformations
  ([`I()`](https://rdrr.io/r/base/AsIs.html)). Set to `NULL` for RE-only
  blocks.

- `fn(w ~ ..., c)`: Defines the weight function (micro-macro link). The
  `c` parameter controls weight normalization: when `c = TRUE`
  (default), weights are normalized to sum to 1 within each group
  (\\\tilde{w}\_{ik} = w\_{ik} / \sum\_{k} w\_{ik}\\). Set `c = FALSE`
  for unnormalized weights (e.g., when aggregating sums). Note: Does not
  support interactions or [`I()`](https://rdrr.io/r/base/AsIs.html) -
  pre-create transformed variables.

- `RE`: Logical; if `TRUE`, include random effects for this block.
  Automatically `TRUE` if `vars = NULL`.

- `ar`: Logical; if `TRUE`, member-level random effects evolve as a
  random walk across repeated participations in groups. This captures
  dynamics where a member's unobserved heterogeneity changes over time.
  Default: `FALSE`.

**Multiple mm() blocks:** You can specify multiple
[`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks with
different weight functions. However, `RE = TRUE` can only be specified
for one [`mm()`](https://benrosche.github.io/bml/reference/mm.md) block.


    mm(id = id(pid, gid), vars = vars(X1), fn = fn(w ~ 1/n), RE = FALSE) +
    mm(id = id(pid, gid), vars = vars(X2), fn = fn(w ~ tenure), RE = FALSE) +
    mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE)

## Hierarchical Membership Object [`hm()`](https://benrosche.github.io/bml/reference/hm.md)


    hm(id = id(hmid), vars = vars(X.hm), name = hmname, type = "RE", showFE = FALSE)

**Components:**

- `id = id(hmid)`: Variable identifying nesting-level groups.

- `vars = vars(X.hm)`: Nesting-level variables, or `NULL`. Supports
  interactions (`*`, `:`) and transformations
  ([`I()`](https://rdrr.io/r/base/AsIs.html)).

- `name = hmname`: Optional labels for nesting-level units.

- `type`: `"RE"` (default) or `"FE"`.

- `showFE`: If `TRUE` and `type = "FE"`, report the fixed effects.

## Supported Families / Links

- Gaussian (continuous): `family = "Gaussian"`

- Binomial (logistic): `family = "Binomial"`

- Weibull survival: `family = "Weibull"`, outcome: `Surv(time, event)`

- Cox survival: `family = "Cox"`, outcome: `Surv(time, event)`

## Priors

Priors can be specified for parameters. With multiple mm() blocks, use
indexed names:


    priors = list(
      "b.mm.1 ~ dnorm(0, 0.01)",
      "b.w.1 ~ dnorm(0, 0.1)",
      "tau.mm ~ dscaled.gamma(25, 1)"
    )

## References

Rosche, B. (2026). A Multilevel Model for Coalition Governments:
Uncovering Party-Level Dependencies Within and Between Governments.
*Political Analysis*.

Browne, W. J., Goldstein, H., & Rasbash, J. (2001). Multiple membership
multiple classification (MMMC) models. *Statistical Modelling*, 1(2),
103-124.

## See also

[`summary.bml`](https://benrosche.github.io/bml/reference/summary.bml.md)
for model summaries,
[`monetPlot`](https://benrosche.github.io/bml/reference/monetPlot.md)
for posterior visualization,
[`mcmcDiag`](https://benrosche.github.io/bml/reference/mcmcDiag.md) for
convergence diagnostics,
[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`hm`](https://benrosche.github.io/bml/reference/hm.md) for model
specification helpers

## Author

Benjamin Rosche \<benrosche@nyu.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
data(coalgov)

# Basic multiple-membership model
# Parties (pid) within governments (gid), nested in countries (cid)
m1 <- bml(
  Surv(govdur, earlyterm) ~ 1 + majority +
    mm(
      id   = id(pid, gid),
      vars = vars(fdep),
      fn   = fn(w ~ 1/n, c = TRUE),
      RE   = TRUE
    ) +
    hm(id = id(cid), type = "RE"),
  family  = "Weibull",
  n.iter  = 10000,
  n.burnin = 5000,
  data    = coalgov
)

# View results
summary(m1)
monetPlot(m1, "b[2]")  # Plot for majority coefficient

# Multiple mm() blocks with different weight functions
m2 <- bml(
  Surv(govdur, earlyterm) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE) +
    mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE),
  family  = "Weibull",
  data    = coalgov
)

# Cox model with piecewise baseline hazard (faster for large datasets)
m3 <- bml(
  Surv(govdur, earlyterm) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE),
  family  = "Cox",
  cox_intervals = 10,  # Use 10 intervals instead of all unique times
  data    = coalgov
)

# Parameterized weight function
m4 <- bml(
  Surv(govdur, earlyterm) ~ 1 + majority +
    mm(
      id   = id(pid, gid),
      vars = vars(fdep),
      fn   = fn(w ~ b0 + b1 * govmaxdur, c = TRUE),  # Weights depend on data
      RE   = FALSE
    ),
  family = "Weibull",
  data   = coalgov
)

# Fixed coefficients (offsets)
m5 <- bml(
  Surv(govdur, earlyterm) ~ 1 + majority +
    mm(
      id   = id(pid, gid),
      vars = vars(fix(fdep, 1.0) + rile),  # Fix fdep coefficient to 1.0
      fn   = fn(w ~ 1/n, c = TRUE),
      RE   = FALSE
    ),
  family = "Weibull",
  data   = coalgov
)

# Custom priors
m6 <- bml(
  Surv(govdur, earlyterm) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE),
  family = "Weibull",
  priors = list(
    "b[1] ~ dnorm(0, 0.01)",       # Intercept prior
    "b.mm.1 ~ dnorm(0, 0.1)",      # MM coefficient prior
    "tau.mm ~ dgamma(2, 0.5)"      # MM precision prior
  ),
  data   = coalgov
)

# Cross-classified model (multiple hm blocks)
m7 <- bml(
  Y ~ 1 + x1 +
    hm(id = id(region), type = "RE") +
    hm(id = id(year), type = "RE"),
  family = "Gaussian",
  data   = mydata
)
} # }
```
