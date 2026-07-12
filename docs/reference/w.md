# Specify the membership weights of an mm() block

Defines the weights \\w\_{ikt}\\ — how much each member counts in the
block's aggregation. `w()` defines a within-group *measure* over
members; the aggregation function
[`fn`](https://benrosche.github.io/bml/reference/fn.md) chooses which
function of that measure (and of the
[`vars`](https://benrosche.github.io/bml/reference/vars.md) attributes)
to read off. Weights can be fixed, taken from observed data, or
estimated from member characteristics:


    w(~ 1/n)                          # equal weights
    w(~ importance, scale = TRUE)     # observed, normalized
    w(~ ilogit(b0 + b1 * q))          # estimated: b0, b1 are free parameters

**One parameter rule (shared with
[`fn`](https://benrosche.github.io/bml/reference/fn.md)):** any symbol
in the formula that is not a data column or reserved word (`n`, the
group size) is a *free parameter* with a default prior. The build always
messages the detected parameters (e.g.
`Estimating parameters: b0, b1 (w)`), so a misspelled column name shows
up immediately. Priors are set via `prior(..., class = "w")`.

**Who counts vs. how contributions combine:** parameters that act on
*external* member characteristics (importance, exposure time) belong in
`w()`; parameters that act on the *effect attribute* (whatever is in
[`vars()`](https://benrosche.github.io/bml/reference/vars.md)) belong in
[`fn`](https://benrosche.github.io/bml/reference/fn.md). Building
weights from the same variable that is being aggregated (e.g.
`w(~ exp(om * x))` with `vars(x)`) is near-equivalent to a smooth-max in
[`fn()`](https://benrosche.github.io/bml/reference/fn.md) and is flagged
with a warning.

## Usage

``` r
w(formula = ~1/n, scale = TRUE)
```

## Arguments

- formula:

  One-sided formula defining the raw weights. Group-level aggregation
  shortcuts (`min(x)`, `max(x)`, `mean(x)`, `sum(x)`, `sd(x)`, `var(x)`,
  `median(x)`, `mode(x)`, `range(x)`, `first(x)`, `last(x)`,
  `quantile(x, p)`) are pre-computed in R. JAGS math functions (`exp`,
  `log`, `ilogit`, ...) pass through.

- scale:

  Logical; if `TRUE` (default), weights are normalized to sum to 1
  within each group, making them a probability measure over the group's
  members. Emergent
  [`fn()`](https://benrosche.github.io/bml/reference/fn.md) types
  require `scale = TRUE`.

## Value

A `bml_w` object.

## Details

Weight functions with estimated parameters must produce bounded,
positive values across plausible parameter values; unbounded weight
functions can crash the sampler. Use bounded forms such as
`w(~ ilogit(b0 + b1 * q), scale = TRUE)` or the generalized logistic
`w(~ 1 / (1 + (n - 1) * exp(-(b0 + b1 * q))), scale = TRUE)`, keep
`scale = TRUE`, standardize covariates, and use informative priors
(`prior(normal(0, 1), class = "w")`). Weight parameters are initialized
at 0 for stable starting values.

## References

Rosche, B. (2026). A Multilevel Model for Coalition Governments:
Uncovering Party-Level Dependencies Within and Between Governments.
*Political Analysis*.

## See also

[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`fn`](https://benrosche.github.io/bml/reference/fn.md),
[`bml`](https://benrosche.github.io/bml/reference/bml.md)

## Examples

``` r
if (FALSE) { # \dontrun{
w(~ 1/n)                                # equal weights
w(~ importance, scale = TRUE)           # observed weights, normalized
w(~ ilogit(b0 + b1 * seniority))        # estimated weights
w(~ duration, scale = FALSE)            # unnormalized (aggregating sums)
} # }
```
