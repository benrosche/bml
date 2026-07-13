# Getting Started with bml

## Introduction

The `bml` package estimates *micro-to-macro regressions*: Bayesian
multilevel models in which the aggregation from member-level records to
group-level outcomes is an explicit, estimable object. Instead of
aggregating lower-level data to higher levels with an assumed rule
(risking aggregation bias) or disaggregating outcomes to lower levels
(artificially inflating sample size), `bml` models how lower-level units
jointly shape higher-level outcomes — including *which members count*
(the weights) and *how their contributions combine* (the aggregation
function).

## Installation

Install the stable version from CRAN:

``` r

install.packages("bml")
```

You’ll also need JAGS installed on your system. See the [installation
vignette](https://benrosche.github.io/bml/articles/installation.md) for
details.

## Basic Example

Let’s start with a simple example using coalition government data. Each
coalition (higher-level unit) is composed of multiple political parties
(lower-level units), and we want to model how party characteristics
influence coalition outcomes.

``` r

library(bml)
data(coalgov)

# Examine the data structure
head(coalgov[, c("gid", "pid", "pname", "n", "finance", "dur_wkb", "event_wkb")])
```

### Understanding the Data Structure

The `coalgov` dataset is in **long format** where each row represents a
party’s participation in a coalition:

- `gid`: Government (coalition) identifier
- `pid`: Party identifier
- `n`: Number of parties in the coalition
- `finance`: Party’s financial dependence on member contributions
  (standardized)
- `dur_wkb`: Coalition duration in days
- `event_wkb`: Early termination indicator (1 = terminated early, 0 =
  censored)

### Model 1: Basic Multiple-Membership Model

Our first model examines coalition duration as a function of party-level
characteristics, aggregated using equal weights:

``` r

mod1 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(
      id   = id(pid, gid),
      vars = vars(finance),
      w    = w(~ 1/n, scale = TRUE),
      fn   = fn("sum"),
      RE   = TRUE
    ),
  data   = coalgov,
  family = weibull(),
  seed   = 1
)

summary(mod1)
```

**Breaking down the formula:**

- `Surv(dur_wkb, event_wkb)`: Survival outcome (duration and event
  indicator)
- `majority`: Government-level covariate (binary indicator)
- [`mm()`](https://benrosche.github.io/bml/reference/mm.md):
  Multiple-membership block containing:
  - `id = id(pid, gid)`: Party ID (member) and government ID (group)
  - `vars = vars(finance)`: Party-level attribute
  - `w = w(~ 1/n, scale = TRUE)`: The weights — *who counts* (equal
    weights, normalized to sum to 1)
  - `fn = fn("sum")`: The aggregation function — *how contributions
    combine* (the additive weighted mean; this is the default)
  - `RE = TRUE`: Party-specific random intercepts (shorthand for
    `RE = re(1)`)

The feature’s coefficient appears in the output as `A_finance` — an
ordinary class-`"b"` coefficient of the main model.

### Model 2: Parameterizing the Weight Function

A key feature of `bml` is the ability to estimate the aggregation
weights from member characteristics. Instead of assuming all parties
have equal influence, we can test whether seat share affects a party’s
weight in the coalition:

``` r

mod2 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(
      id   = id(pid, gid),
      vars = vars(finance),
      w    = w(~ 1 / (1 + (n - 1) * exp(-(b0 + b1 * pseat))), scale = TRUE),
      fn   = fn("sum"),
      RE   = TRUE
    ),
  data   = coalgov,
  family = weibull(),
  prior  = prior(normal(0, 1), class = "w"),  # weakly informative prior on the weight parameters
  seed   = 1
)

summary(mod2)
```

**New elements:**

- The weight formula includes `pseat` (party seat share). By the *one
  parameter rule*, the symbols `b0` and `b1` are not data columns and
  therefore become free parameters — the build messages
  `Estimating parameters: b0, b1 (w.1)`.
- `prior = prior(normal(0, 1), class = "w")`: A structured prior on the
  weight parameters, written on the standard-deviation scale and
  translated to JAGS automatically. Run
  `get_prior(formula, data, family)` to list everything that is
  settable.

**Interpretation:**

- If \\b_1 \approx 0\\: Seat share doesn’t affect weights (equal
  influence)
- If \\b_1 \> 0\\: Larger parties have more weight
- If \\b_1 \< 0\\: Smaller parties have more weight

### Model 3: Beyond the sum — emergent features

The aggregation function
[`fn()`](https://benrosche.github.io/bml/reference/fn.md) is not limited
to the weighted mean. Emergent features are properties of the *whole
member set*:

``` r

mod3 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(finance), w = w(~ 1/n), fn = fn("sum"), RE = TRUE) +      # mean:     A_finance
    mm(id = id(pid, gid), vars = vars(finance), w = w(~ 1/n), fn = fn("var")) + # spread:   V_finance
    mm(id = id(pid, gid), w = w(~ pseat, scale = TRUE), fn = fn("hhi"),         # dominance: C_w
       RE = FALSE),
  data   = coalgov,
  family = weibull(),
  seed   = 1
)
```

Other built-ins include `fn("threshold", c = est())` (critical mass with
an estimated cutpoint), `fn("smax", kappa = est())` (is the aggregation
mean-like or max/min-like?), `fn("gmean", p = est())` (CES/power mean),
and `fn("cov")`. You can also write your own reduction with the
expression DSL: `fn(~ E((x - E(x))^2))`.

### Visualizing Results

``` r

# Coefficients (all class-"b" terms, labeled by feature name)
fixef(mod2)

# Diagnostic plot for a weight parameter
monetPlot(mod2, parameter = "b.w.1[1]", label = "Seat share effect")

# MCMC diagnostics
mcmcDiag(mod2, parameters = "b.w.1[1]")

# The posterior toolchain
draws <- as_draws_df(mod2)
```

## Next Steps

- See the [Model
  documentation](https://benrosche.github.io/bml/articles/model.md) for
  mathematical details
- Explore
  [Examples](https://benrosche.github.io/bml/articles/examples.md) for
  applications
- Check the [FAQ](https://benrosche.github.io/bml/articles/faq.md) for
  common questions
- Read the accompanying paper: [Rosche
  (2026)](https://osf.io/preprints/socarxiv/4bafr)

## Key Concepts

- **Multiple-membership structure:** Higher-level units (coalitions)
  contain multiple lower-level units (parties), and lower-level units
  can appear in multiple higher-level units. `id(member, group)` uses
  long format: one row per membership, any number of members per group.
- **Weights [`w()`](https://benrosche.github.io/bml/reference/w.md) —
  who counts:** a within-group measure over members. Fixed (`w(~ 1/n)`),
  observed (`w(~ importance, scale = TRUE)`), or estimated
  (`w(~ ilogit(b0 + b1 * q))`). `scale = TRUE` normalizes weights to sum
  to 1 within each group.
- **Aggregation function
  [`fn()`](https://benrosche.github.io/bml/reference/fn.md) — how
  contributions combine:** `fn("sum")` is the additive case;
  `fn("var")`, `fn("hhi")`, `fn("threshold")`, `fn("smax")`,
  `fn("gmean")`, `fn("cov")` are emergent features; `fn(~ E(...))` is
  the user-written DSL. Parameters that act on external member
  characteristics belong in
  [`w()`](https://benrosche.github.io/bml/reference/w.md); parameters
  that act on the effect attribute belong in
  [`fn()`](https://benrosche.github.io/bml/reference/fn.md).
- **Effects:** member/unit random effects with `RE = re(1 + x)` (partial
  pooling; `cor = TRUE` opts into correlated intercept-slope), or fixed
  effects with `FE = fe(1 + x)` (no pooling).
  [`hm()`](https://benrosche.github.io/bml/reference/hm.md) adds
  hierarchical nesting; unit-level covariates go in the main formula
  (`Y ~ gdp + hm(id(cid), RE = re(1 + gdp))`).
- **Named blocks:** give a block a `name =` and reference its feature in
  interactions (`Ax:education`, `Ax:Vx`).
