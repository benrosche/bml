# Widely applicable information criterion for bml models

Widely applicable information criterion for bml models

## Usage

``` r
# S3 method for class 'bml'
waic(x, ...)
```

## Arguments

- x:

  A fitted `bml` model (see
  [`loo.bml`](https://benrosche.github.io/bml/reference/loo.bml.md) for
  availability).

- ...:

  Passed to
  [`loo::waic.matrix()`](https://mc-stan.org/loo/reference/waic.html).

## Value

A `waic` object.
