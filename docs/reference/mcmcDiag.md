# Summarize MCMC convergence diagnostics

Computes common convergence diagnostics for selected parameters from a
JAGS/BUGS fit and returns a compact, report-ready table. The diagnostics
include Gelman–Rubin \\\hat{R}\\, Geweke z-scores, Heidelberger-Welch
stationarity p-values, and autocorrelation at a user-specified lag.

## Usage

``` r
mcmcDiag(bml.out, parameters, lag = 50)
```

## Arguments

- bml.out:

  A model fit object containing JAGS output, typically as returned by
  [`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html), with
  component `$jags.out$BUGSoutput`.

- parameters:

  Character vector of parameter names (or patterns) to extract. These
  may be exact names or patterns (e.g., a prefix like `"b"` that matches
  `"b[1]"`, `"b[2]"`, …).

- lag:

  Integer specifying the lag at which to compute autocorrelation.
  Default: 50. Lower values (e.g., 10) capture short-range dependence;
  higher values assess whether the chain has mixed well over longer
  intervals.

## Value

A `data.frame` with one row per diagnostic and one column per parameter;
cell entries are the average diagnostic values across chains. Row names
include: `"Gelman/Rubin convergence statistic"`, `"Geweke z-score"`,
`"Heidelberger/Welch p-value"`, `"Autocorrelation (lag <lag>)"`.

## Details

Internally, the function converts the BUGS/JAGS output to a
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html), then
computes per-chain diagnostics and averages them across chains for each
parameter:

- **Gelman–Rubin** (\\\hat{R}\\):
  [`coda::gelman.diag()`](https://rdrr.io/pkg/coda/man/gelman.diag.html).
  Values close to 1 indicate convergence; a common heuristic is
  \\\hat{R} \le 1.1\\.

- **Geweke** z-score:
  [`coda::geweke.diag()`](https://rdrr.io/pkg/coda/man/geweke.diag.html).
  Large absolute values (e.g., \\\|z\|\>2\\) suggest lack of
  convergence.

- **Heidelberger–Welch** p-value:
  [`coda::heidel.diag()`](https://rdrr.io/pkg/coda/man/heidel.diag.html)
  tests the null of stationarity in the chain segment.

- **Autocorrelation**:
  [`coda::autocorr()`](https://rdrr.io/pkg/coda/man/autocorr.html) at
  the specified `lag`, averaged across chains. Values near zero indicate
  good mixing; persistent autocorrelation suggests the chain needs
  thinning or reparameterization.

All statistics are rounded to three decimals. The returned table is
transposed so that *rows are diagnostics* and *columns are parameters*.

## References

Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
using multiple sequences. *Statistical Science*, 7(4), 457-472.

Brooks, S. P., & Gelman, A. (1998). General methods for monitoring
convergence of iterative simulations. *Journal of Computational and
Graphical Statistics*, 7(4), 434-455.

## See also

[`gelman.diag`](https://rdrr.io/pkg/coda/man/gelman.diag.html),
[`geweke.diag`](https://rdrr.io/pkg/coda/man/geweke.diag.html),
[`heidel.diag`](https://rdrr.io/pkg/coda/man/heidel.diag.html),
[`autocorr`](https://rdrr.io/pkg/coda/man/autocorr.html)

## Author

Benjamin Rosche <benrosche@nyu.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
data(coalgov)

# Fit model
m1 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE),
  family = "Weibull",
  monitor = TRUE,
  data = coalgov
)

# Check convergence for main parameters
mcmcDiag(m1, parameters = "b")  # All b coefficients

# Check specific parameters
mcmcDiag(m1, parameters = c("b[1]", "b[2]", "shape"))

# Check mm block parameters
mcmcDiag(m1, parameters = c("b.mm.1", "sigma.mm.1"))

# Custom autocorrelation lag
mcmcDiag(m1, parameters = "b", lag = 100)

# Interpreting results:
# - Gelman-Rubin < 1.1: Good convergence
# - |Geweke z| < 2: No evidence against convergence
# - Heidelberger p > 0.05: Chain appears stationary
# - Low autocorrelation: Good mixing
} # }
```
