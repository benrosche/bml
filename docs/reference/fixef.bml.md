# Extract fixed coefficients from a bml model

Returns all class-`"b"` coefficients — main-formula terms, block
features, and named-feature interactions — labeled by term, one row per
coefficient. `fn` shape parameters are deliberately excluded (they are
not coefficients; see the "Function parameters" section of
[`summary.bml`](https://benrosche.github.io/bml/reference/summary.bml.md)).

## Usage

``` r
# S3 method for class 'bml'
fixef(object, ...)
```

## Arguments

- object:

  A fitted `bml` model.

- ...:

  Unused.

## Value

A matrix with columns `Estimate`, `Est.Error`, `Q2.5`, `Q97.5` and term
rownames.
