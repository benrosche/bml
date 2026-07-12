# List every settable prior of a bml model

Returns the table of all parameters whose priors can be set via
[`prior`](https://benrosche.github.io/bml/reference/prior.md), together
with their default priors — the bml analogue of `brms::get_prior()`. Use
it to discover valid `class`/`coef`/ `group` combinations before
overriding.

## Usage

``` r
get_prior(formula, data, family = stats::gaussian())
```

## Arguments

- formula:

  The model formula (same as for
  [`bml`](https://benrosche.github.io/bml/reference/bml.md)).

- data:

  The data (member-level long format).

- family:

  Model family (function, object, or string). Default
  [`gaussian()`](https://rdrr.io/r/stats/family.html).

## Value

A data frame with columns `class`, `coef`, `group`, `node` (the internal
JAGS node), and `default`.

## See also

[`prior`](https://benrosche.github.io/bml/reference/prior.md),
[`bml`](https://benrosche.github.io/bml/reference/bml.md)

## Examples

``` r
if (FALSE) { # \dontrun{
get_prior(
  y ~ x + mm(id = id(pid, gid), vars = vars(z), w = w(~ 1/n)),
  data = dat, family = gaussian()
)
} # }
```
