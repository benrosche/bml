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

- `logit`, `ilogit`, `probit`, `iprobit`, `cloglog`, `icloglog`

- `round`, `trunc`, `floor`, `ceiling`

Example: `fn(w ~ 1 / (1 + (n - 1) * exp(-(b1 * x))))` uses an
exponential decay function where weights depend on member
characteristics.

**Ensuring Numerical Stability:**

Weight functions with estimated parameters (`b0`, `b1`, ...) must
produce bounded, positive values across all plausible parameter values.
Unbounded weight functions can cause the MCMC sampler to crash (e.g.,
`"Error in node w.1[...]: Invalid parent values"`). During sampling,
weight parameters can take on extreme values, and if the weight function
is not bounded, this will destabilize the likelihood.

Recommendations:

- **Use bounded weight functions.** Two options:

  - `ilogit()`: Bounds weights between 0 and 1 with a zero-point at 0.5:
    `fn(w ~ ilogit(b0 + b1 * x), c = TRUE)`

  - **Generalized logistic** (Rosche, 2026): Bounds weights between 0
    and 1 with a zero-point at \\1/n\\ (equal weights), so deviations
    from equal weighting are estimated as a function of covariates:
    `fn(w ~ 1 / (1 + (n - 1) * exp(-(b0 + b1 * x))), c = TRUE)`

- **Use `c = TRUE`** (weight normalization) to prevent weights from
  growing without bound

- **Standardize covariates** in the weight function. Variables with
  large ranges (e.g., income in thousands) can cause `b * x` to overflow

- **Use informative priors** for weight parameters via the `priors`
  argument in [`bml`](https://benrosche.github.io/bml/reference/bml.md)
  (e.g., `priors = list("b.w.1[1] ~ dnorm(0, 1)")`)

- **Avoid unbounded functions** like `exp(b * x)` without normalization
  (`c = TRUE`) or wrapping (e.g., inside `ilogit()`)

Weight parameters are initialized at 0 by default to ensure numerically
stable starting values. See
[`vignette("faq")`](https://benrosche.github.io/bml/articles/faq.md)
(Question 7) for detailed troubleshooting of numerical issues.

## References

Rosche, B. (2026). A Multilevel Model for Theorizing and Estimating the
Micro-Macro Link. *Political Analysis*.

Browne, W. J., Goldstein, H., & Rasbash, J. (2001). Multiple membership
multiple classification (MMMC) models. *Statistical Modelling*, 1(2),
103-124.

## See also

[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`vignette("model")`](https://benrosche.github.io/bml/articles/model.md)
for the model structure,
[`vignette("faq")`](https://benrosche.github.io/bml/articles/faq.md) for
troubleshooting

## Examples

``` r
if (FALSE) { # \dontrun{
# Equal weights (standard multiple-membership)
fn(w ~ 1/n, c = TRUE)

# Tenure-based weights (proportional to time served)
fn(w ~ tenure, c = TRUE)

# Flexible parameterized weights
fn(w ~ b0 + b1 * seniority, c = TRUE)

# Unconstrained weights
fn(w ~ importance, c = FALSE)

# Weights based on group aggregates
fn(w ~ b1 * min(tenure) + (1 - b1) * mean(tenure), c = TRUE)

# Combining individual and aggregate measures
fn(w ~ b0 + b1 * (tenure / max(tenure)), c = TRUE)

# Using median for robust central tendency
fn(w ~ tenure / median(tenure), c = TRUE)

# Using quantiles for percentile-based weights
fn(w ~ quantile(tenure, 0.75) - quantile(tenure, 0.25), c = TRUE)
} # }
```
