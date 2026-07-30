# Predictions from a bml model

Posterior predictive summaries (one row per observation). For the raw
draws use
[`posterior_predict.bml`](https://benrosche.github.io/bml/reference/posterior_predict.md);
for the linear predictor use
[`fitted.bml`](https://benrosche.github.io/bml/reference/fitted.bml.md).

## Usage

``` r
# S3 method for class 'bml'
predict(object, summary = TRUE, ...)
```

## Arguments

- object:

  A fitted `bml` model with `monitor = "predictive"`.

- summary:

  Logical; if `TRUE` (default), return a summary matrix, else the draws
  matrix.

- ...:

  Unused (`newdata` is not supported yet).

## Value

A matrix of posterior predictive summaries, or a draws matrix.
