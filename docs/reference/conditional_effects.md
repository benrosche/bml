# Conditional effects of main-formula predictors

Marginal-effect curves for main-formula covariates: the expected outcome
on the linear-predictor scale across a grid of the chosen predictor,
holding the other main-formula covariates at their means and the block
features, random effects, and interactions at their sample-average
contribution.

## Usage

``` r
conditional_effects(x, ...)

# S3 method for class 'bml'
conditional_effects(x, effects = NULL, resolution = 100, ...)
```

## Arguments

- x:

  A fitted `bml` model with `monitor = c("parameters", "fitted")`.

- ...:

  Unused.

- effects:

  Character vector of main-formula term names. Default: all numeric
  main-formula terms.

- resolution:

  Grid size. Default: 100.

## Value

A named list of data frames (one per effect) with columns `value`,
`estimate`, `lb`, `ub`, suitable for plotting.
