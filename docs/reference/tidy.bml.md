# Tidy a fitted bml model

Turns a fitted [`bml`](https://benrosche.github.io/bml/reference/bml.md)
model into a tidy tibble following the broom conventions: one row per
parameter with columns `term`, `estimate`, `std.error`, `conf.low`, and
`conf.high`. Because the model is Bayesian, `estimate` is the posterior
mean, `std.error` the posterior standard deviation, and
`conf.low`/`conf.high` the bounds of an equal-tailed credible interval.

With these methods, ecosystem tools that build on broom work out of the
box, e.g. `modelsummary::modelsummary(list(m1, m2))` for side-by-side
regression tables or `dotwhisker::dwplot(m1)` for coefficient plots.

## Usage

``` r
# S3 method for class 'bml'
tidy(
  x,
  conf.int = TRUE,
  conf.level = 0.95,
  component = c("all", "fixed", "random", "weights", "fn"),
  ...
)
```

## Arguments

- x:

  A fitted model object of class `"bml"` returned by
  [`bml`](https://benrosche.github.io/bml/reference/bml.md).

- conf.int:

  Logical. Include credible interval columns `conf.low`/`conf.high`?
  Default: `TRUE`.

- conf.level:

  Width of the equal-tailed credible interval. Default: 0.95, which is
  read directly from the fitted model. Other levels are computed from
  the posterior draws and therefore require `monitor = "parameters"`.

- component:

  Which parameters to return:

  - `"all"` (default): all of the below

  - `"fixed"`: all class-`"b"` coefficients (main-formula terms, block
    features, interactions)

  - `"random"`: variance parameters (`sd(...)`, `cor(...)`, `sigma`,
    Weibull `shape`)

  - `"weights"`: weight-model parameters (`w[...]`)

  - `"fn"`: aggregation-function shape parameters (`fn[...]`)

- ...:

  Additional arguments (currently unused).

## Value

A `tibble` with columns `term`, `estimate`, `std.error`, `conf.low`,
`conf.high` (if `conf.int`), and `component`.

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`glance.bml`](https://benrosche.github.io/bml/reference/glance.bml.md),
[`bmlCompare`](https://benrosche.github.io/bml/reference/bmlCompare.md),
[`coefPlot`](https://benrosche.github.io/bml/reference/coefPlot.md),
[`summary.bml`](https://benrosche.github.io/bml/reference/summary.bml.md)

## Author

Benjamin Rosche <benrosche@nyu.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
data(coalgov)

m1 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("sum"), RE = TRUE),
  family = weibull(),
  data = coalgov
)

tidy(m1)                      # regression coefficients
tidy(m1, component = "all")   # including weight and variance parameters
tidy(m1, conf.level = 0.9)    # requires monitor = "parameters"

# Multi-model comparison table (see also bmlCompare()):
models <- list(base = m1, weighted = m2)
purrr::imap_dfr(models, \(m, name)
  tidy(m) |>
    dplyr::mutate(model = name, N = glance(m)$nobs, DIC = glance(m)$DIC)
)
} # }
```
