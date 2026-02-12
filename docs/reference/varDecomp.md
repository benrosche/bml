# Variance decomposition for fitted bml models

Computes a posterior variance decomposition and intraclass correlation
coefficients (ICCs) from a fitted `bml` model. The function
automatically discovers all variance components (sigma parameters) in
the model, applies weight adjustments for multiple-membership levels,
and returns posterior summaries.

## Usage

``` r
varDecomp(model, uncertainty = "sd", r = 2)
```

## Arguments

- model:

  A fitted model object of class `"bml"` returned by
  [`bml`](https://benrosche.github.io/bml/reference/bml.md). Must have
  been fitted with `monitor = TRUE` (the default).

- uncertainty:

  Uncertainty measure to report. One of `"sd"` (posterior standard
  deviation, the default), `"mad"` (median absolute deviation), or
  `"ci"` (95% credible interval with lower/upper bounds).

- r:

  Number of decimal places for rounding numeric output. Default: 2.

## Value

A data frame of class `"bml_varDecomp"` with one row per variance
component. Always includes `Component`, `sigma`, and `ICC` columns.
Additional columns depend on `uncertainty`:

- `"sd"`: `sigma_sd` and `ICC_sd`

- `"mad"`: `sigma_mad` and `ICC_mad`

- `"ci"`: `sigma_lb`, `sigma_ub`, `ICC_lb`, `ICC_ub`

## Details

**Variance decomposition.** The total variance of the outcome is
partitioned into additive components, one for each level in the model:

\$\$\mathrm{Var}(y) = \sigma_1^2 + \sigma_2^2 + \ldots\$\$

Each component contributes variance \\\sigma^2\\, except for
multiple-membership (MM) levels, where the effective variance
contribution is scaled by the average of the summed squared weights
across groups.

\$\$\mathrm{Var}\_{\mathrm{mm}} = \sigma\_{\mathrm{mm}}^2 \cdot
\overline{w^2}, \quad \overline{w^2} = \frac{1}{N} \sum\_{i=1}^{N}
\sum\_{k} w\_{ik}^2\$\$

This weight adjustment accounts for the fact that the member-level
variance is distributed across multiple members with potentially unequal
influence. With equal weights (\\w\_{ik} = 1/n_i\\), the effective
variance shrinks as group size increases.

**Intraclass Correlation Coefficient (ICC).** The ICC for a given level
is the proportion of total variance attributable to that level:

\$\$\rho_l = \frac{\sigma_l^2}{\sum\_{l'} \sigma\_{l'}^2}\$\$

Intuitively, the ICC answers: "What fraction of the total variation in
the outcome is due to differences between units at this level?" An ICC
of 0.30 for the country level, for example, means that 30% of the
outcome variation can be attributed to between-country differences.

ICCs are computed per posterior draw and then summarized, properly
propagating uncertainty from the MCMC samples.

**Family-specific handling:**

- **Gaussian / Weibull**: The residual `sigma` from the model is used
  directly.

- **Binomial**: There is no residual sigma. The latent logistic residual
  variance \\\pi^2/3 \approx 3.29\\ is used instead.

- **Cox**: There is no residual variance. ICCs are computed among the
  non-residual components only.

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`summary.bml`](https://benrosche.github.io/bml/reference/summary.bml.md)

## Author

Benjamin Rosche \<benrosche@nyu.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
data(coalgov)

m1 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE) +
    hm(id = id(cid), type = "RE"),
  family = "Weibull",
  data = coalgov
)

varDecomp(m1)
varDecomp(m1, uncertainty = "ci")
varDecomp(m1, uncertainty = "mad")
} # }
```
