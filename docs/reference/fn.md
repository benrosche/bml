# Specify a weight function for multiple-membership models

Defines how member-level contributions are weighted when aggregating to
the group level (the "micro-macro link"). The weight function can be a
simple formula (e.g., `1/n` for equal weights) or can include parameters
to be estimated from the data.

## Usage

``` r
fn(w = w ~ 1/n, c = TRUE)
```

## Arguments

- w:

  A two-sided formula specifying the weight function. The left-hand side
  must be `w`; the right-hand side defines the weighting scheme:

  - Simple: `w ~ 1/n` (equal weights based on group size)

  - Parameterized: `w ~ b0 + b1 * tenure` (weights depend on member
    characteristics and estimated parameters)

  - With group aggregates: `w ~ b1 * min(x) + (1-b1) * mean(x)` (weights
    based on group-level summaries; see Details)

  Parameters must be named `b0`, `b1`, `b2`, etc.

- c:

  Logical; if `TRUE` (default), weights are normalized to sum to 1
  within each group. Set to `FALSE` for unnormalized weights.

## Value

A `bml_fn` object containing the parsed weight function specification.

## Details

**Weight Function Components:**

- **Variables** (e.g., `n`, `tenure`): Data from your dataset

- **Parameters** (e.g., `b0`, `b1`): Estimated from the data

- **Operations**: Standard R arithmetic (`+`, `-`, `*`, `/`, `^`, etc.)

**Common Weight Functions:**

- Equal weights: `w ~ 1/n`

- Duration-based: `w ~ duration`

- Flexible parameterized: `w ~ b0 + b1 * seniority`

- Group aggregates: `w ~ b1 * min(x) + (1-b1) * mean(x)`

When `c = TRUE`, the weights are constrained: \\\sum\_{k \in group} w_k
= 1\\.

**Group-Level Aggregation Functions:**

The weight function supports aggregation functions that compute
summaries within each group (mainid). These are pre-computed in R before
passing to JAGS. Supported functions:

- `min(var)`, `max(var)`: Minimum/maximum value within the group

- `mean(var)`, `sum(var)`: Mean/sum of values within the group

- `median(var)`, `mode(var)`: Median/mode (most frequent) value within
  the group

- `sd(var)`, `var(var)`, `range(var)`: Standard deviation/variance/range
  (max-min) within the group

- `first(var)`, `last(var)`: First/last value (based on data order)

- `quantile(var, prob)`: Quantile at probability `prob` (0 to 1). For
  example, `quantile(x, 0.25)` computes the 25th percentile.

Example: `fn(w ~ b1 * min(tenure) + (1-b1) * max(tenure))` creates
weights that blend the minimum and maximum tenure within each group,
with the blend controlled by the estimated parameter `b1`.

Example with quantile: `fn(w ~ quantile(tenure, 0.75) / max(tenure))`
uses the 75th percentile relative to the maximum within each group.

Note: Nested aggregation functions (e.g., `min(max(x))`) are not
supported.

**JAGS Mathematical Functions:**

The following mathematical functions are passed directly to JAGS and can
be used in weight formulas:

- `exp`, `log`, `log10`, `sqrt`, `abs`, `pow`

- `sin`, `cos`, `tan`, `asin`, `acos`, `atan`

- `sinh`, `cosh`, `tanh`

- `round`, `trunc`, `floor`, `ceiling`

Example: `fn(w ~ 1 / (1 + (n - 1) * exp(-(b1 * x))))` uses an
exponential decay function where weights depend on member
characteristics. See Rosche (2026) for more details on parameterized
weight functions.

## References

Browne, W. J., Goldstein, H., & Rasbash, J. (2001). Multiple membership
multiple classification (MMMC) models. *Statistical Modelling*, 1(2),
103-124.

## See also

[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`bml`](https://benrosche.github.io/bml/reference/bml.md)

## Examples

``` r
# Equal weights (standard multiple-membership)
fn(w ~ 1/n, c = TRUE)
#> $formula
#> w ~ 1/n
#> <environment: 0x0000022d9249c9a8>
#> 
#> $string
#> [1] "1/n"
#> 
#> $vars
#> [1] "n"
#> 
#> $vars_p
#> character(0)
#> 
#> $params
#> character(0)
#> 
#> $constraint
#> [1] TRUE
#> 
#> $agg_funcs
#> NULL
#> 
#> $agg_vars
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_fn"

# Tenure-based weights (proportional to time served)
fn(w ~ tenure, c = TRUE)
#> $formula
#> w ~ tenure
#> <environment: 0x0000022d9249c9a8>
#> 
#> $string
#> [1] "tenure"
#> 
#> $vars
#> [1] "tenure"
#> 
#> $vars_p
#> character(0)
#> 
#> $params
#> character(0)
#> 
#> $constraint
#> [1] TRUE
#> 
#> $agg_funcs
#> NULL
#> 
#> $agg_vars
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_fn"

# Flexible parameterized weights
fn(w ~ b0 + b1 * seniority, c = TRUE)
#> $formula
#> w ~ b0 + b1 * seniority
#> <environment: 0x0000022d9249c9a8>
#> 
#> $string
#> [1] "b0 * X0 + b1 * seniority"
#> 
#> $vars
#> [1] "X0"        "seniority"
#> 
#> $vars_p
#> [1] "X0"        "seniority"
#> 
#> $params
#> [1] "b0" "b1"
#> 
#> $constraint
#> [1] TRUE
#> 
#> $agg_funcs
#> NULL
#> 
#> $agg_vars
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_fn"

# Unconstrained weights
fn(w ~ importance, c = FALSE)
#> $formula
#> w ~ importance
#> <environment: 0x0000022d9249c9a8>
#> 
#> $string
#> [1] "importance"
#> 
#> $vars
#> [1] "importance"
#> 
#> $vars_p
#> character(0)
#> 
#> $params
#> character(0)
#> 
#> $constraint
#> [1] FALSE
#> 
#> $agg_funcs
#> NULL
#> 
#> $agg_vars
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_fn"

# Weights based on group aggregates
fn(w ~ b1 * min(tenure) + (1 - b1) * mean(tenure), c = TRUE)
#> $formula
#> w ~ b1 * min(tenure) + (1 - b1) * mean(tenure)
#> <environment: 0x0000022d9249c9a8>
#> 
#> $string
#> [1] "b1 * tenure_min + (1 - b1) * tenure_mean"
#> 
#> $vars
#> [1] "tenure_min"  "tenure_mean"
#> 
#> $vars_p
#> [1] "tenure_min"
#> 
#> $params
#> [1] "b1"
#> 
#> $constraint
#> [1] TRUE
#> 
#> $agg_funcs
#> $agg_funcs[[1]]
#> $agg_funcs[[1]]$original
#> [1] "min(tenure)"
#> 
#> $agg_funcs[[1]]$func
#> [1] "min"
#> 
#> $agg_funcs[[1]]$var
#> [1] "tenure"
#> 
#> $agg_funcs[[1]]$prob
#> NULL
#> 
#> $agg_funcs[[1]]$col_name
#> [1] "tenure_min"
#> 
#> 
#> $agg_funcs[[2]]
#> $agg_funcs[[2]]$original
#> [1] "mean(tenure)"
#> 
#> $agg_funcs[[2]]$func
#> [1] "mean"
#> 
#> $agg_funcs[[2]]$var
#> [1] "tenure"
#> 
#> $agg_funcs[[2]]$prob
#> NULL
#> 
#> $agg_funcs[[2]]$col_name
#> [1] "tenure_mean"
#> 
#> 
#> 
#> $agg_vars
#> [1] "tenure"
#> 
#> attr(,"class")
#> [1] "bml_fn"

# Combining individual and aggregate measures
fn(w ~ b0 + b1 * (tenure / max(tenure)), c = TRUE)
#> $formula
#> w ~ b0 + b1 * (tenure/max(tenure))
#> <environment: 0x0000022d9249c9a8>
#> 
#> $string
#> [1] "b0 * X0 + b1 * (tenure/tenure_max)"
#> 
#> $vars
#> [1] "X0"         "tenure"     "tenure_max"
#> 
#> $vars_p
#> [1] "X0"
#> 
#> $params
#> [1] "b0" "b1"
#> 
#> $constraint
#> [1] TRUE
#> 
#> $agg_funcs
#> $agg_funcs[[1]]
#> $agg_funcs[[1]]$original
#> [1] "max(tenure)"
#> 
#> $agg_funcs[[1]]$func
#> [1] "max"
#> 
#> $agg_funcs[[1]]$var
#> [1] "tenure"
#> 
#> $agg_funcs[[1]]$prob
#> NULL
#> 
#> $agg_funcs[[1]]$col_name
#> [1] "tenure_max"
#> 
#> 
#> 
#> $agg_vars
#> [1] "tenure"
#> 
#> attr(,"class")
#> [1] "bml_fn"

# Using median for robust central tendency
fn(w ~ tenure / median(tenure), c = TRUE)
#> $formula
#> w ~ tenure/median(tenure)
#> <environment: 0x0000022d9249c9a8>
#> 
#> $string
#> [1] "tenure/tenure_median"
#> 
#> $vars
#> [1] "tenure"        "tenure_median"
#> 
#> $vars_p
#> character(0)
#> 
#> $params
#> character(0)
#> 
#> $constraint
#> [1] TRUE
#> 
#> $agg_funcs
#> $agg_funcs[[1]]
#> $agg_funcs[[1]]$original
#> [1] "median(tenure)"
#> 
#> $agg_funcs[[1]]$func
#> [1] "median"
#> 
#> $agg_funcs[[1]]$var
#> [1] "tenure"
#> 
#> $agg_funcs[[1]]$prob
#> NULL
#> 
#> $agg_funcs[[1]]$col_name
#> [1] "tenure_median"
#> 
#> 
#> 
#> $agg_vars
#> [1] "tenure"
#> 
#> attr(,"class")
#> [1] "bml_fn"

# Using quantiles for percentile-based weights
fn(w ~ quantile(tenure, 0.75) - quantile(tenure, 0.25), c = TRUE)
#> $formula
#> w ~ quantile(tenure, 0.75) - quantile(tenure, 0.25)
#> <environment: 0x0000022d9249c9a8>
#> 
#> $string
#> [1] "tenure_q75 - tenure_q25"
#> 
#> $vars
#> [1] "tenure_q75" "tenure_q25"
#> 
#> $vars_p
#> character(0)
#> 
#> $params
#> character(0)
#> 
#> $constraint
#> [1] TRUE
#> 
#> $agg_funcs
#> $agg_funcs[[1]]
#> $agg_funcs[[1]]$original
#> [1] "quantile(tenure, 0.75)"
#> 
#> $agg_funcs[[1]]$func
#> [1] "quantile"
#> 
#> $agg_funcs[[1]]$var
#> [1] "tenure"
#> 
#> $agg_funcs[[1]]$prob
#> [1] 0.75
#> 
#> $agg_funcs[[1]]$col_name
#> [1] "tenure_q75"
#> 
#> 
#> $agg_funcs[[2]]
#> $agg_funcs[[2]]$original
#> [1] "quantile(tenure, 0.25)"
#> 
#> $agg_funcs[[2]]$func
#> [1] "quantile"
#> 
#> $agg_funcs[[2]]$var
#> [1] "tenure"
#> 
#> $agg_funcs[[2]]$prob
#> [1] 0.25
#> 
#> $agg_funcs[[2]]$col_name
#> [1] "tenure_q25"
#> 
#> 
#> 
#> $agg_vars
#> [1] "tenure"
#> 
#> attr(,"class")
#> [1] "bml_fn"
```
