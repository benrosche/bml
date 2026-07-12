# Model families for bml

Family functions specify the outcome distribution and link function of a
[`bml`](https://benrosche.github.io/bml/reference/bml.md) model,
following the brms/[`glm`](https://rdrr.io/r/stats/glm.html) convention
of passing a family *function* rather than a string:

- [`gaussian()`](https://rdrr.io/r/stats/family.html): normal outcome,
  identity link

- `bernoulli()` (alias
  [`binomial()`](https://rdrr.io/r/stats/family.html)): binary outcome,
  logit link

- `weibull()`: Weibull survival model; outcome `Surv(time, event)`

- `cox(intervals = NULL)`: Cox proportional hazards model; outcome
  `Surv(time, event)`. `intervals` controls the baseline hazard: `NULL`
  uses all unique event times (non-parametric); an integer k uses a
  piecewise-constant baseline hazard with k intervals (faster).

Family *strings* are accepted as first-class input too (as in brms):
`family = "gaussian"`, `"bernoulli"`, `"weibull"`, `"cox"`.

## Usage

``` r
bernoulli(link = "logit")

weibull(link = "log")

cox(intervals = NULL, link = "log")
```

## Arguments

- link:

  Link function. Only the canonical links are supported by the JAGS
  backend (`"identity"` for gaussian, `"logit"` for bernoulli, `"log"`
  for weibull/cox); passing anything else errors.

- intervals:

  For `cox()`: `NULL` for a non-parametric baseline hazard (all unique
  event times), or a positive integer number of intervals for a
  piecewise-constant baseline hazard.

## Value

A `bml_family` object (a list with elements `family` and `link`, plus
`intervals` for `cox()`).

## Examples

``` r
if (FALSE) { # \dontrun{
bml(y ~ x + mm(...), data = dat, family = gaussian())
bml(y01 ~ x + mm(...), data = dat, family = bernoulli())
bml(Surv(t, ev) ~ x + mm(...), data = dat, family = weibull())
bml(Surv(t, ev) ~ x + mm(...), data = dat, family = cox(intervals = 10))
} # }
```
