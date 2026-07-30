# Extract random effects from a bml model

Returns the posterior means of the unit-level random effects: member
effects per [`mm()`](https://benrosche.github.io/bml/reference/mm.md)
member-id group (intercepts, and slopes when `re(1 + x)` was used) and
unit effects per
[`hm()`](https://benrosche.github.io/bml/reference/hm.md) block.

## Usage

``` r
# S3 method for class 'bml'
ranef(object, ...)
```

## Arguments

- object:

  A fitted `bml` model with `monitor = "random_effects"`.

- ...:

  Unused.

## Value

A list with elements `mm` and `hm`.
