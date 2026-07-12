# Extract combined coefficients from a bml model

Following brms semantics loosely: the fixed coefficients plus, per
grouping structure, the unit-level effects (fixed part + random
deviation).

## Usage

``` r
# S3 method for class 'bml'
coef(object, ...)
```

## Arguments

- object:

  A fitted `bml` model.

- ...:

  Unused.

## Value

A list with elements `fixef` (matrix) and `ranef` (list).
