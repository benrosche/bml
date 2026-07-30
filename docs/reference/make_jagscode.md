# Generate the JAGS model code / data without fitting

The bml analogues of `brms::make_stancode()` / `make_standata()`: build
the model exactly as
[`bml`](https://benrosche.github.io/bml/reference/bml.md) would, but
return the generated JAGS model string (`make_jagscode`) or the JAGS
data list (`make_jagsdata`) instead of fitting.

## Usage

``` r
make_jagscode(
  formula,
  data,
  family = stats::gaussian(),
  prior = NULL,
  monitor = "summary"
)

make_jagsdata(
  formula,
  data,
  family = stats::gaussian(),
  prior = NULL,
  monitor = "summary"
)
```

## Arguments

- formula, data, family, prior:

  As in [`bml`](https://benrosche.github.io/bml/reference/bml.md).

- monitor:

  As in [`bml`](https://benrosche.github.io/bml/reference/bml.md)
  (affects which nodes the model defines).

## Value

`make_jagscode()`: a character scalar of class `"bml_jagscode"` (prints
nicely). `make_jagsdata()`: a named list.

## Examples

``` r
if (FALSE) { # \dontrun{
code <- make_jagscode(
  y ~ x + mm(id = id(pid, gid), vars = vars(z), w = w(~ 1/n), fn = fn("sum")),
  data = dat, family = gaussian()
)
code
} # }
```
