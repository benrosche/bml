# Glance at a fitted bml model

Returns a one-row tibble with model-level information following the
broom conventions, suitable for model comparison.

## Usage

``` r
# S3 method for class 'bml'
glance(x, ...)
```

## Arguments

- x:

  A fitted model object of class `"bml"` returned by
  [`bml`](https://benrosche.github.io/bml/reference/bml.md).

- ...:

  Additional arguments (currently unused).

## Value

A one-row `tibble` with columns:

- `nobs`: Number of main-level (outcome) units

- `n.members`: Number of member-level units across mm blocks

- `DIC`: Deviance Information Criterion

- `family`: Outcome family

- `iter`, `chains`: MCMC settings

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`tidy.bml`](https://benrosche.github.io/bml/reference/tidy.bml.md),
[`bmlCompare`](https://benrosche.github.io/bml/reference/bmlCompare.md)

## Author

Benjamin Rosche <benrosche@nyu.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
glance(m1)
} # }
```
