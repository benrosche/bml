# Posterior expectations of the linear predictor

Posterior expectations of the linear predictor

## Usage

``` r
# S3 method for class 'bml'
fitted(object, summary = TRUE, ...)
```

## Arguments

- object:

  A fitted `bml` model with `monitor = "fitted"`.

- summary:

  Logical; if `TRUE` (default), return a summary matrix (one row per
  observation), else the draws matrix (draws x observations).

- ...:

  Unused.

## Value

A matrix of posterior summaries of `mu`, or a draws matrix.
