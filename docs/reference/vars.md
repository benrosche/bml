# Specify member-level attributes for mm() blocks

Helper function used within
[`mm`](https://benrosche.github.io/bml/reference/mm.md) to specify the
member-level attributes (the effect covariates \\x\_{kt}\\) that the
block's aggregation function
[`fn`](https://benrosche.github.io/bml/reference/fn.md) reads. Supports
both free variables (with coefficients estimated by the main model) and
fixed variables (coefficients held constant via
[`fix`](https://benrosche.github.io/bml/reference/fix.md)).

## Usage

``` r
vars(...)
```

## Arguments

- ...:

  Unquoted variable names, combined with `+` (formula-style) or as
  separate arguments. Supports:

  - Simple variables: `vars(x + y)` or `vars(x, y)`

  - Interactions: `vars(x * y)` (main effects + interaction) or
    `vars(x:y)` (interaction only; the member-paired interaction
    \\H\_{xz} = \sum_k w_k x_k z_k\\)

  - Transformations: `vars(I(x^2))`

  - Fixed coefficients: `vars(fix(x, 1.0) + y)`

  Note: Numeric literals like `1`, `0`, or `-1` are ignored (no
  intercept support in mm blocks).

## Value

A `bml_vars` object (list with `formula`, `free`, `fixed`), or `NULL` if
empty.

## See also

[`fix`](https://benrosche.github.io/bml/reference/fix.md),
[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`fn`](https://benrosche.github.io/bml/reference/fn.md)

## Examples

``` r
if (FALSE) { # \dontrun{
vars(income)
vars(income + education)
vars(income:education)      # member-paired interaction only
vars(x, z)                  # two attributes (e.g. for fn("cov"))
vars(fix(exposure, 1.0) + income)
} # }
```
