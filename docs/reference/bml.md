# Bayesian Multilevel Models for Micro-Macro Analysis (Multiple Membership) Using JAGS

The **bml** package estimates micro-to-macro regressions: Bayesian
multiple-membership multilevel models in which the aggregation from
member-level records to group-level outcomes is an explicit, estimable
object — weights
([`w`](https://benrosche.github.io/bml/reference/w.md)), aggregation
functions ([`fn`](https://benrosche.github.io/bml/reference/fn.md):
mean, variance, concentration, thresholds, smooth max/min, CES means,
covariance, or a user-written reduction), member/unit random and fixed
effects ([`re`](https://benrosche.github.io/bml/reference/re.md),
[`fe`](https://benrosche.github.io/bml/reference/fe.md)), and
cross-level or feature-by-feature interactions via named blocks.

JAGS must be installed separately:
<https://sourceforge.net/projects/mcmc-jags/>.

## Usage

``` r
bml(
  formula,
  data,
  family = stats::gaussian(),
  prior = NULL,
  inits = NULL,
  iter = 1000,
  warmup = 500,
  thin = max(1, floor((iter - warmup)/1000)),
  chains = 3,
  cores = 1,
  seed = NULL,
  run = TRUE,
  monitor = "summary",
  modelfile = FALSE,
  ...
)
```

## Arguments

- formula:

  A symbolic model formula; see Details.

- data:

  Data frame in member-level (long) format: one row per membership. Must
  contain all variables referenced in the formula, including the
  identifiers in
  [`id()`](https://benrosche.github.io/bml/reference/id.md).

- family:

  Model family: a family function
  ([`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`bernoulli()`](https://benrosche.github.io/bml/reference/bml-families.md),
  [`weibull()`](https://benrosche.github.io/bml/reference/bml-families.md),
  [`cox()`](https://benrosche.github.io/bml/reference/bml-families.md))
  or a string (`"gaussian"`, `"bernoulli"`, `"weibull"`, `"cox"`). Cox
  baseline-hazard intervals travel with the family:
  `cox(intervals = 10)`.

- prior:

  Prior specifications built with
  [`prior`](https://benrosche.github.io/bml/reference/prior.md) and
  combined with [`c()`](https://rdrr.io/r/base/c.html); see
  [`get_prior`](https://benrosche.github.io/bml/reference/get_prior.md)
  for the settable parameters. Raw JAGS strings (e.g.
  `"b.fn.1[1] ~ dnorm(0, 0.1)"`) pass through untranslated as an escape
  hatch.

- inits:

  List of initial values for MCMC chains (applied to all chains). Weight
  and `fn` shape parameters get stable defaults automatically; user
  values override.

- iter:

  Total number of MCMC iterations per chain. Default: 1000.

- warmup:

  Number of warmup (burn-in) iterations discarded from the start of each
  chain. Default: 500.

- thin:

  Thinning rate. Default targets ~1000 retained draws.

- chains:

  Number of MCMC chains. Default: 3.

- cores:

  Number of cores; `cores > 1` runs chains in parallel. Default: 1.

- seed:

  Integer random seed for reproducibility.

- run:

  Logical; if `FALSE`, returns the generated JAGS model string, data,
  and monitors without fitting (see also
  [`make_jagscode`](https://benrosche.github.io/bml/reference/make_jagscode.md)).

- monitor:

  Character vector selecting retained posterior capabilities. The
  default, `"summary"`, retains posterior summaries and convergence
  diagnostics but no draws. Add `"parameters"`, `"random_effects"`,
  `"weights"`, `"fitted"`, `"predictive"`, or `"log_lik"` as needed.
  `"full"` stores all family-supported draws once in a compact
  [`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html);
  `"jags"` explicitly retains the raw R2jags object.

- modelfile:

  Logical or character path: `TRUE` saves the generated JAGS code to
  `modelstring.txt`; a path reads JAGS code from that file instead of
  generating it.

- ...:

  Not used; catches removed legacy arguments (`n.iter`, `n.burnin`,
  `n.thin`, `n.chains`, `parallel`, `priors`, `cox_intervals`) with a
  migration message.

## Value

A list of class `"bml"` with elements `reg.table` (posterior summaries;
labeled terms), `w` (weight matrices per block), `re.mm`/`re.hm` (random
effects), `pred` (posterior predictive means), `fitted` (posterior means
of the linear predictor), `input` (model metadata), `diagnostics`,
`storage`, `draws` (requested compact posterior draws), and `jags.out`
(non-`NULL` only when `monitor = "jags"`).

## Details

The general formula structure is


    outcome ~ 1 + predictors + mm(...) + hm(...)

where [`mm`](https://benrosche.github.io/bml/reference/mm.md) defines a
multiple-membership block (the micro-macro link) and
[`hm`](https://benrosche.github.io/bml/reference/hm.md) a hierarchical
nesting level. Multiple
[`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks stack
features (e.g. mean + variance + concentration); multiple
[`hm()`](https://benrosche.github.io/bml/reference/hm.md) blocks give
cross-classification. For survival families use `Surv(time, event)` on
the left-hand side.

Named [`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks
can be referenced in main-formula interactions:


    Y ~ education + Ax:education +
      mm(name = Ax, id = id(task, occ), vars = vars(x), w = w(~ importance), fn = fn("sum"))

Both cross-level (feature x macro variable) and block x block (feature x
feature) interactions are supported; the macro variable must also appear
as a main effect, and a bare (non-interacted) block reference is an
error.

The package is introduced in Rosche (2026), *Political Analysis*.

## References

Rosche, B. (2026). A Multilevel Model for Coalition Governments:
Uncovering Party-Level Dependencies Within and Between Governments.
*Political Analysis*.

Browne, W. J., Goldstein, H., & Rasbash, J. (2001). Multiple membership
multiple classification (MMMC) models. *Statistical Modelling*, 1(2),
103-124.

## See also

[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`hm`](https://benrosche.github.io/bml/reference/hm.md),
[`w`](https://benrosche.github.io/bml/reference/w.md),
[`fn`](https://benrosche.github.io/bml/reference/fn.md),
[`re`](https://benrosche.github.io/bml/reference/re.md),
[`fe`](https://benrosche.github.io/bml/reference/fe.md),
[`prior`](https://benrosche.github.io/bml/reference/prior.md),
[`get_prior`](https://benrosche.github.io/bml/reference/get_prior.md),
[`summary.bml`](https://benrosche.github.io/bml/reference/summary.bml.md),
[`fixef.bml`](https://benrosche.github.io/bml/reference/fixef.bml.md),
[`ranef.bml`](https://benrosche.github.io/bml/reference/ranef.bml.md),
[`loo.bml`](https://benrosche.github.io/bml/reference/loo.bml.md),
[`varDecomp`](https://benrosche.github.io/bml/reference/varDecomp.md),
[`bmlCompare`](https://benrosche.github.io/bml/reference/bmlCompare.md),
[`make_jagscode`](https://benrosche.github.io/bml/reference/make_jagscode.md)

## Author

Benjamin Rosche \<benrosche@nyu.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
data(coalgov)

# Additive aggregation with member random effects, country nesting
m1 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("sum"), RE = TRUE) +
    hm(id = id(cid)),
  data = coalgov,
  family = weibull()
)
summary(m1)

# Estimated weights (b0, b1 are free parameters by the one-parameter rule)
m2 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion),
       w = w(~ ilogit(b0 + b1 * pseat), scale = TRUE),
       fn = fn("sum")),
  data = coalgov,
  family = weibull()
)

# Emergent features: variance + concentration, stacked
m3 <- bml(
  sim.y ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("var")) +
    mm(id = id(pid, gid), w = w(~ pseat, scale = TRUE), fn = fn("hhi")),
  data = coalgov,
  family = gaussian()
)

# Structured priors
m4 <- bml(
  sim.y ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("sum"), RE = TRUE),
  data = coalgov,
  family = gaussian(),
  prior = c(
    prior(normal(0, 5), class = "b"),
    prior(exponential(1), class = "sd")
  )
)
} # }
```
