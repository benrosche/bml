# Visualize posterior distributions with density and trace plots

Creates a combined diagnostic plot showing both the posterior density
and MCMC trace plot for a specified parameter. Helps assess convergence
and visualize posterior uncertainty. The plot displays the median and
90% highest posterior density (HPD) interval.

## Usage

``` r
monetPlot(bml, parameter, label = NULL, r = 2, yaxis = T)
```

## Arguments

- bml:

  A fitted model object of class `"bml"` returned by
  [`bml`](https://benrosche.github.io/bml/reference/bml.md). Must be
  fitted with `monitor = TRUE` to store MCMC chains.

- parameter:

  Character string specifying the parameter to plot. Must use the
  internal parameter name (i.e., row names from `bml$reg.table`).
  Examples: `"b[1]"` (intercept), `"b[2]"` (first covariate), `"b.mm.1"`
  (first mm block coefficient), `"sigma.mm"` (mm random effect SD).

- label:

  Optional character string for the parameter label displayed on the
  plot. If `NULL` (default), uses the internal parameter name.

- r:

  Number of decimal places for displayed quantiles and statistics.
  Default: 2.

- yaxis:

  Logical; if `TRUE` (default), display axis titles ("Density" and
  "Scans"). If `FALSE`, omit axis titles for cleaner appearance when
  combining multiple plots.

## Value

A `ggplot` object (using `patchwork`) combining two panels:

- **Top panel**: Posterior density with shaded 90% HPD interval. Solid
  vertical line at zero, dashed line at posterior median.

- **Bottom panel**: Trace plot showing MCMC iterations across chains.
  Same reference lines as top panel. Helps diagnose convergence and
  mixing.

## Details

**Interpreting the Plot:**

- **Density panel**: Shows the posterior distribution. The dashed line
  marks the median (central estimate). Shading indicates the 90%
  credible region.

- **Trace panel**: Shows parameter values across MCMC iterations for
  each chain. Good mixing looks like "fuzzy caterpillars" with chains
  overlapping. Poor mixing shows trends, stickiness, or separation
  between chains.

**Convergence Checks:**

- Chains should overlap and explore the same space

- No sustained trends or drift

- Rapid mixing (no long autocorrelation)

Use [`mcmcDiag`](https://benrosche.github.io/bml/reference/mcmcDiag.md)
for formal convergence statistics (Gelman-Rubin, Geweke, etc.).

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`mcmcDiag`](https://benrosche.github.io/bml/reference/mcmcDiag.md),
[`summary.bml`](https://benrosche.github.io/bml/reference/summary.bml.md)

## Author

Benjamin Rosche \<benrosche@nyu.edu\>

## Examples

``` r
if (FALSE) { # \dontrun{
data(coalgov)

# Fit model with monitoring enabled
m1 <- bml(
  Surv(govdur, earlyterm) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE) +
    hm(id = id(cid), type = "RE"),
  family = "Weibull",
  monitor = TRUE,  # Required for monetPlot
  data = coalgov
)

# Plot intercept
monetPlot(m1, parameter = "b[1]", label = "Intercept")

# Plot majority coefficient with custom label
monetPlot(m1, parameter = "b[2]", label = "Majority Government Effect")

# Plot mm coefficient
monetPlot(m1, parameter = "b.mm.1", label = "Party Fragmentation")

# Plot random effect SD
monetPlot(m1, parameter = "sigma.mm")

# List available parameters
rownames(m1$reg.table)
} # }
```
