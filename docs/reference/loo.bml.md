# Efficient approximate leave-one-out cross-validation for bml models

PSIS-LOO via the loo package, computed from the retained pointwise
log-likelihood. Available for gaussian and bernoulli families fitted
with `monitor = "log_lik"`. Compare models with
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).

## Usage

``` r
# S3 method for class 'bml'
loo(x, ...)
```

## Arguments

- x:

  A fitted `bml` model.

- ...:

  Passed to
  [`loo::loo.matrix()`](https://mc-stan.org/loo/reference/loo.html).

## Value

A `loo` object.
