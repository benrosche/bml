# Fixed effects: no pooling

The no-pooling peer of
[`re`](https://benrosche.github.io/bml/reference/re.md): identical
formula grammar, but unit effects get flat independent priors instead of
a hierarchical prior (Gelman-Hill: "FE = RE with infinite variance").
`fe(1)` gives classic unit dummies; `fe(1 + x)` additionally gives
unit-specific slopes. The first unit is the reference category (its
effects are fixed to 0), so `fe()` remains identified alongside a global
intercept and global slopes.

## Usage

``` r
fe(x = 1, showFE = FALSE)
```

## Arguments

- x:

  The effects expression: `1`, `1 + var`, `0 + var`.

- showFE:

  Logical; if `TRUE`, per-unit estimates are monitored and reported.
  Default `FALSE`.

## Value

A `bml_fe` object.

## See also

[`re`](https://benrosche.github.io/bml/reference/re.md),
[`hm`](https://benrosche.github.io/bml/reference/hm.md),
[`mm`](https://benrosche.github.io/bml/reference/mm.md)
