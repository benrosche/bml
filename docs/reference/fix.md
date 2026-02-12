# Fix a coefficient to a known value

Specify a covariate whose coefficient should be held constant at a fixed
value rather than estimated from the data. This is useful for offset
variables or when you want to impose theoretical constraints. Fixed
coefficients are handled efficiently by pre-computing their contribution
in R before passing data to JAGS.

## Usage

``` r
fix(var, value)
```

## Arguments

- var:

  Unquoted variable name from your data

- value:

  Numeric value for the coefficient (e.g., `1.0` for a standard offset)

## Value

A `bml_fix` object that can be used within
[`vars`](https://benrosche.github.io/bml/reference/vars.md).

## See also

[`vars`](https://benrosche.github.io/bml/reference/vars.md),
[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`hm`](https://benrosche.github.io/bml/reference/hm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fix a coefficient to 1.0 (standard offset)
fix(exposure, 1.0)

# Use within vars() for multiple-membership models
vars(fix(population, 0.5) + income + education)
} # }
```
