# Summarize MCMC convergence diagnostics

Computes common convergence diagnostics for selected parameters from a
JAGS/BUGS fit and returns a compact, report-ready table. The diagnostics
include Gelman–Rubin \\\hat{R}\\, Geweke z-scores, Heidelberger-Welch
stationarity p-values, and autocorrelation at lag 50.

## Usage

``` r
mcmcDiag(bml.out, parameters)
```

## Arguments

- bml.out:

  A model fit object containing JAGS output, typically as returned by
  [`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html), with
  component `$jags.out$BUGSoutput`.

- parameters:

  Character vector of parameter names (or patterns) to extract, passed
  to `crMCMC`. Depending on your `crMCMC()` implementation, these may be
  exact names or patterns (e.g., a prefix like `"b.l2"` that matches
  `"b.l2[1]"`, `"b.l2[2]"`, …).

## Value

A `data.frame` with one row per diagnostic and one column per parameter;
cell entries are the average diagnostic values across chains. Row names
include: `"Gelman/Rubin convergence statistic"`, `"Geweke z-score"`,
`"Heidelberger/Welch p-value"`, `"Autocorrelation (lag 50)"`.

## Details

Internally, the function converts the BUGS/JAGS output to a
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) via
`crMCMC`, then computes per-chain diagnostics and averages them across
chains for each parameter:

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

- **Autocorrelation (lag 50)**:
  [`coda::autocorr()`](https://rdrr.io/pkg/coda/man/autocorr.html) at
  lag 50, averaged across chains.

All statistics are rounded to three decimals. The returned table is
transposed so that *rows are diagnostics* and *columns are parameters*.

## See also

`crMCMC`,
[`gelman.diag`](https://rdrr.io/pkg/coda/man/gelman.diag.html),
[`geweke.diag`](https://rdrr.io/pkg/coda/man/geweke.diag.html),
[`heidel.diag`](https://rdrr.io/pkg/coda/man/heidel.diag.html),
[`autocorr`](https://rdrr.io/pkg/coda/man/autocorr.html)

## Author

Benjamin Rosche <benrosche@nyu.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
mcmcdiag(fit, parameters = c("beta", "sigma"))
} # }
```
