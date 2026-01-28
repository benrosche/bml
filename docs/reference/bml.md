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
  chains = 3,
  seed = NULL,
  run = TRUE,
  parallel = FALSE,
  monitor = TRUE,
  modelfile = FALSE,
  data = NULL
)
```

## Arguments

- formula:

  A symbolic model formula.

- family:

  A character string specifying the model family.

- priors:

  A named list or character vector mapping parameter groups to JAGS
  priors.

- inits:

  A list of initial values used for all chains.

- n.iter:

  Total number of MCMC iterations.

- n.burnin:

  Number of burn-in iterations to discard.

- n.thin:

  Thinning rate for MCMC sampling.

- chains:

  Number of MCMC chains to run.

- seed:

  Random seed for reproducibility.

- run:

  Logical; if `TRUE` (default), JAGS is executed.

- parallel:

  Logical; if `TRUE`, run chains in parallel.

- monitor:

  Logical; if `TRUE`, store additional outputs.

- modelfile:

  Character or logical for saving/loading JAGS model.

- data:

  A data frame where each row represents a level-1 observation.

## Value

A list of class "bml" containing model outputs.

## Details

In addition to hierarchical and cross-classified multilevel models, the
**bml** package allows users to fit Bayesian multiple-membership models.
Unlike tools such as `brms` or MLwiN
(<https://www.bristol.ac.uk/cmm/software/mlwin/>), **bml** lets users
specify and estimate models in which membership weights are
parameterized through flexible formula syntax. This enables a more
nuanced examination of how effects from lower-level units aggregate to
higher levels (the micro-macro link).

The package automatically generates JAGS code to fit the model and
processes the output to facilitate interpretation of model parameters
and diagnostics.

For accessible introductions to multiple-membership models, see Fielding
and Goldstein (2006) and Beretvas (2010). Advanced treatments include
Goldstein (2011, Ch. 13), Rasbash and Browne (2001, 2008), Browne et al.
(2001), and Leckie (2013). The package and modeling framework are
introduced in: Rosche, B. (2025). *A Multilevel Model for Coalition
Governments: Uncovering Party-Level Dependencies Within and Between
Governments*. *Political Analysis*.

## Formula Components

- **Outcome (Y):** The dependent variable. For survival models, use
  `Surv(time, event)`.

- **Intercept (1):** Includes an intercept term; use `0` to omit it.

- **Level-2 predictors (X.L2):** Variables defined at level 2, separated
  by `+`.

- **Level-3 predictors (X.L3):** Variables defined at level 3, separated
  by `+`.

- **Multiple membership object
  ([`mm()`](https://benrosche.github.io/bml/reference/mm.md)):** Defines
  how lower-level (level-1) units are associated with higher-level
  (level-2) constructs using a user-specified weighting function.
  Multiple [`mm()`](https://benrosche.github.io/bml/reference/mm.md)
  objects can be specified with different weight functions.

- **Hierarchical membership
  ([`hm()`](https://benrosche.github.io/bml/reference/hm.md)):**
  Specifies nesting of level-2 units within level-3 entities.
  Cross-classified structures can be modeled by including multiple
  [`hm()`](https://benrosche.github.io/bml/reference/hm.md) objects.

**Important:** The formula parser does not support
[`I()`](https://rdrr.io/r/base/AsIs.html) inside `bml()`. Create
transformations and interactions in your data object before modeling.

## Multiple Membership Object [`mm()`](https://benrosche.github.io/bml/reference/mm.md)


    mm(
      id   = id(l1id, l2id),
      vars = vars(X.L1),
      fn   = fn(w ~ 1/n, c = TRUE, ar = FALSE),
      RE   = TRUE
    )

**Components:**

- `id(l1id, l2id)`: Specifies identifiers linking each level-1 unit
  (`l1id`) to its corresponding level-2 entities (`l2id`).

- `vars(X.L1)`: Specifies level-1 covariates aggregated across
  memberships. Use `+` to include multiple variables. Set to `NULL` for
  RE-only blocks.

- `fn(w ~ ..., c, ar)`: Defines the weight function (micro-macro link).

- `RE`: Logical; if `TRUE`, include random effects for this block.
  Automatically `TRUE` if `vars = NULL`.

**Multiple mm() blocks:** You can specify multiple
[`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks with
different weight functions:


    mm(id = id(pid, gid), vars = vars(X1), fn = fn(w ~ 1/n), RE = FALSE) +
    mm(id = id(pid, gid), vars = vars(X2), fn = fn(w ~ max(n)), RE = FALSE) +
    mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE)

## Hierarchical Membership Object [`hm()`](https://benrosche.github.io/bml/reference/hm.md)


    hm(id = id(l3id), vars = vars(X.L3), name = l3name, type = "RE", showFE = FALSE)

**Components:**

- `id = id(l3id)`: Variable identifying level-3 groups.

- `vars = vars(X.L3)`: Level-3 variables, or `NULL`.

- `name = l3name`: Optional labels for level-3 units.

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
      "b.l1.1 ~ dnorm(0, 0.01)",
      "b.w.1 ~ dnorm(0, 0.1)",
      "tau.l1 ~ dscaled.gamma(25, 1)"
    )

## Author

Benjamin Rosche \<benrosche@nyu.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
data(coalgov)

# Single mm() block
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
  data    = coalgov
)

# Multiple mm() blocks with different weight functions
m2 <- bml(
  Y ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE) +
    mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE) +
    hm(id = id(cid), type = "RE"),
  family  = "Gaussian",
  data    = coalgov
)
} # }
```
