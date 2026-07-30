# Posterior predictive checks for bml models

Graphical posterior predictive checks via bayesplot (must be installed).
Compares the observed outcome with replicated outcomes from
[`posterior_predict.bml`](https://benrosche.github.io/bml/reference/posterior_predict.md).

## Usage

``` r
pp_check(object, ...)

# S3 method for class 'bml'
pp_check(object, type = NULL, ndraws = 30, ...)
```

## Arguments

- object:

  A fitted `bml` model (gaussian or bernoulli, fitted with
  `monitor = "predictive"`).

- ...:

  Passed to the bayesplot function.

- type:

  A bayesplot PPC type (without the `"ppc_"` prefix), e.g.
  `"dens_overlay"` (default for gaussian), `"bars"` (default for
  bernoulli), `"hist"`, `"stat"`.

- ndraws:

  Number of replicated datasets to plot. Default: 30.

## Value

A ggplot object.
