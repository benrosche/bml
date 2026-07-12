# Random effects: partial pooling

Specifies the random-effect structure of an
[`mm`](https://benrosche.github.io/bml/reference/mm.md) or
[`hm`](https://benrosche.github.io/bml/reference/hm.md) block, mirroring
the left-hand side of an lme4/brms `( ... | group )` term. The grouping
is implied by the block's
[`id()`](https://benrosche.github.io/bml/reference/id.md) and therefore
not repeated.


    RE = TRUE            # shorthand for re(1)
    re(1)                # random intercept
    re(1 + x)            # intercept + slope on x (INDEPENDENT, the default)
    re(1 + x, cor = TRUE)# intercept + slope, correlated (lme4/brms `|`) - opt-in
    re(0 + x)            # slope only (also re(x - 1))
    re(1, ar = TRUE)     # random walk over participation order
    re(1, ar = year)     # random walk over calendar time (gap-scaled variance)

The default `cor = FALSE` deliberately diverges from brms's correlated
default: the correlated draw is a real JAGS cost, and in
[`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks the
member random effects are only read through the weighted aggregate
`sum_k w_k u_k`, so the correlation is rarely of interest and often
weakly identified. Opt in with `cor = TRUE` (supported for a single
slope, i.e. the 2x2 case).

## Usage

``` r
re(x = 1, cor = FALSE, ar = FALSE)
```

## Arguments

- x:

  The effects expression: `1`, `1 + var`, `0 + var`, or `var - 1`. For
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks,
  slope variables are member-level columns; for
  [`hm()`](https://benrosche.github.io/bml/reference/hm.md) blocks they
  are main-level columns varying within the nesting units.

- cor:

  Logical; if `TRUE`, intercept and slope get a joint bivariate-normal
  prior (correlated). Default `FALSE` (independent).

- ar:

  Autoregressive random walk for the intercept across a unit's repeated
  participations/observations. `FALSE` (default): independent effects.
  `TRUE`: random walk over participation order, \\u_t \sim N(u\_{t-1},
  \sigma^2)\\. An unquoted column name (e.g. `ar = year`): random walk
  over that (numeric) time variable, with the step variance scaled by
  the normalized time gaps, \\u_t \sim N(u\_{t-1}, \sigma^2 \\ \Delta_t
  / \bar\Delta)\\ — equal spacing reproduces the `ar = TRUE` model, and
  \\\sigma\\ keeps its per-average-step interpretation. Requires an
  intercept-only `re(1)`; time values must be unique within each unit.

## Value

A `bml_re` object.

## See also

[`fe`](https://benrosche.github.io/bml/reference/fe.md),
[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`hm`](https://benrosche.github.io/bml/reference/hm.md)
