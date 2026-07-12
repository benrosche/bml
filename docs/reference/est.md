# Mark a shape parameter as estimated

Used inside [`fn`](https://benrosche.github.io/bml/reference/fn.md)'s
named-shortcut arguments to say "estimate this parameter" rather than
fixing it to a number: `fn("threshold", c = est(), kappa = 10)`. In the
formula DSL a free symbol alone declares a parameter, so `est()` is only
needed where a value slot must distinguish estimate-vs-fix.

## Usage

``` r
est(lower = NULL, upper = NULL)
```

## Arguments

- lower, upper:

  Optional bounds; if both given, the default prior is
  `dunif(lower, upper)`. Otherwise the parameter gets a
  parameter-specific default (e.g. the observed range of `x` for a
  threshold cutpoint). Override via `prior(..., class = "fn")`.

## Value

A `bml_est` marker object.

## See also

[`fn`](https://benrosche.github.io/bml/reference/fn.md),
[`prior`](https://benrosche.github.io/bml/reference/prior.md)
