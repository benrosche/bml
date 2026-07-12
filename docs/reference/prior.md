# Structured prior specifications

Specify priors brms-style: a distribution call plus a parameter `class`,
optionally narrowed to a single coefficient (`coef`) or block (`group`).
Distributions are written on the natural (standard-deviation) scale and
translated to JAGS's precision parameterization automatically.


    prior(normal(0, 10),    class = "b")                      # all fixed coefficients
    prior(normal(0, 1),     class = "b", coef = "education")  # one coefficient
    prior(exponential(1),   class = "sd")                     # RE standard deviations
    prior(lkj(2),           class = "cor")                    # RE correlation (re(cor=TRUE))
    prior(student_t(3,0,5), class = "Intercept")
    prior(normal(0, 1),     class = "w")                      # weight-model parameters
    prior(uniform(0, 1),    class = "fn", coef = "c")         # fn shape parameters

Combine several with [`c()`](https://rdrr.io/r/base/c.html) and pass as
`bml(..., prior = c(...))`. Raw JAGS strings remain a power-user escape
hatch: character elements such as `"b.fn.1[1] ~ dnorm(0, 0.1)"` are
passed through untranslated.

## Usage

``` r
prior(dist, class = "b", coef = "", group = "")
```

## Arguments

- dist:

  A distribution call. Supported: `normal(mu, sigma)`,
  `student_t(nu, mu, sigma)`, `cauchy(mu, sigma)`, `exponential(rate)`,
  `gamma(shape, rate)`, `uniform(lower, upper)`, `lkj(eta)` (class
  `"cor"` only). Scale parameters are standard deviations, not
  precisions.

- class:

  Parameter class. One of:

  - `"Intercept"`: the main intercept

  - `"b"`: fixed coefficients (main-formula terms and block features);
    narrow with `coef` (term label, e.g. `"education"` or a feature name
    like `"A_x"`)

  - `"sd"`: random-effect standard deviations; narrow with `group`
    (`"mm"`/`"hm"` or a block id like `"mm.1"`) and `coef`
    (`"Intercept"` or a slope variable)

  - `"cor"`: the RE intercept-slope correlation (only exists when the
    user opts into `re(..., cor = TRUE)`) — use `lkj(eta)`

  - `"sigma"`: the Gaussian residual SD

  - `"shape"`: the Weibull shape parameter

  - `"w"`: weight-model parameters (free symbols in
    [`w()`](https://benrosche.github.io/bml/reference/w.md))

  - `"fn"`: [`fn()`](https://benrosche.github.io/bml/reference/fn.md)
    shape parameters (`c`, `kappa`, `p`, DSL free symbols); narrow with
    `coef`

- coef:

  Optional coefficient/parameter name to narrow the prior to.

- group:

  Optional block identifier (e.g. `"mm.1"`, `"hm.2"`, or just
  `"mm"`/`"hm"`).

## Value

A `bml_prior` data frame (one row per prior); combine with
[`c()`](https://rdrr.io/r/base/c.html).

## See also

[`get_prior`](https://benrosche.github.io/bml/reference/get_prior.md),
[`bml`](https://benrosche.github.io/bml/reference/bml.md)
