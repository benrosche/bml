# Pointwise log-likelihood of a bml model

Pointwise log-likelihood of a bml model

## Usage

``` r
log_lik(object, ...)

# S3 method for class 'bml'
log_lik(object, ...)
```

## Arguments

- object:

  A fitted `bml` model. Available for gaussian and bernoulli families
  fitted with `monitor = TRUE` (the JAGS model monitors a pointwise
  `log_lik` node).

- ...:

  Unused.

## Value

A draws x observations matrix of pointwise log-likelihood values.
