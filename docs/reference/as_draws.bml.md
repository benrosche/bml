# Extract posterior draws from a bml model

Converts the retained MCMC chains into the posterior package's draws
formats, unlocking the whole posterior/bayesplot toolchain
(`summarise_draws()`, `rhat()`, `ess_bulk()`, ...). Variables keep their
internal JAGS node names (`b[1]`, `b.fn.1[1]`, `sigma.mm.1`, ...); the
mapping to term labels is in `model$reg.table` (rownames = node,
`Parameter` = label).

## Usage

``` r
# S3 method for class 'bml'
as_draws(x, ...)

# S3 method for class 'bml'
as_draws_df(x, ...)

# S3 method for class 'bml'
as_draws_matrix(x, ...)

# S3 method for class 'bml'
as_draws_array(x, ...)
```

## Arguments

- x:

  A fitted `bml` model with retained posterior draws.

- ...:

  Passed on to the posterior converters.

## Value

A `draws_array` (or the format of the called variant) containing the
capabilities requested when the model was fitted.
