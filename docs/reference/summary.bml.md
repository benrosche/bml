# Summarize a fitted bml model

S3 method for summarizing `bml` model objects. Returns a formatted table
of parameter estimates with posterior means, standard deviations, and
credible intervals, along with model information and convergence
statistics.

## Usage

``` r
# S3 method for class 'bml'
summary(object, r = 3, diagnostics = TRUE, ...)
```

## Arguments

- object:

  A fitted model object of class `"bml"` returned by
  [`bml`](https://benrosche.github.io/bml/reference/bml.md).

- r:

  Number of decimal places for rounding numeric output. Default: 3.

- diagnostics:

  Logical; include per-parameter R-hat and effective sample sizes in the
  printed table? Default: `TRUE`.

- ...:

  Additional arguments (currently unused).

## Value

A data frame of class `"bml_summary"` containing rounded parameter
estimates with the following columns:

- `Parameter`: Labeled parameter names

- `Estimate`: Posterior mean

- `Est.Error`: Posterior standard deviation

- `Q2.5`/`Q97.5`: 95% credible interval

- `Rhat`, `ESS_bulk`, and `ESS_tail`: convergence diagnostics when
  `diagnostics = TRUE`

The object includes metadata attributes printed above the table:

- Outcome family and link function

- Estimate type (posterior mean from MCMC)

- Credible interval specification (95% equal-tailed)

- Level specification (mm and hm block details)

- DIC (Deviance Information Criterion) for model comparison

## Details

The summary method rounds all numeric values for readability while
preserving the underlying structure and metadata from the fitted model.
All columns remain accessible via standard data frame indexing (e.g.,
`$Parameter`, `$Estimate`).

For Cox models with piecewise baseline hazards (when `cox_intervals` is
specified), the outcome description includes the number of intervals
used.

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`tidy.bml`](https://benrosche.github.io/bml/reference/tidy.bml.md),
[`coefPlot`](https://benrosche.github.io/bml/reference/coefPlot.md),
[`monetPlot`](https://benrosche.github.io/bml/reference/monetPlot.md),
[`mcmcDiag`](https://benrosche.github.io/bml/reference/mcmcDiag.md)

## Author

Benjamin Rosche \<benrosche@nyu.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
data(coalgov)

# Fit model
m1 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("sum"), RE = TRUE) +
    hm(id = id(cid)),
  family = weibull(),
  data = coalgov
)

# View summary
summary(m1)

# Summary with more decimal places
summary(m1, r = 4)

# Access specific columns
s <- summary(m1)
s$Parameter  # Parameter names
s$Estimate   # Posterior means
s$Q2.5       # Lower credible bounds

# Custom posterior summaries (requires monitor = "parameters")
draws <- posterior::as_draws_df(m1)

# Select specific parameters and compute custom summaries
draws |>
  dplyr::select(dplyr::starts_with("b[")) |>
  tidyr::pivot_longer(everything(), names_to = "param") |>
  dplyr::group_by(param) |>
  dplyr::summarise(
    median = median(value),
    mad    = mad(value),
    q05    = quantile(value, 0.05),
    q95    = quantile(value, 0.95)
  )
} # }
```
