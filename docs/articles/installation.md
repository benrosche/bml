# Installation

`bml` fits Bayesian multiple-membership multilevel models via **JAGS**,
so installation has two parts:

1.  Install JAGS (system software)
2.  Install the `bml` R package (from CRAN or GitHub)

## 1. Install JAGS

Install the latest JAGS version for your operating system from:

- <http://mcmc-jags.sourceforge.net/>

After installing JAGS, restart R so `rjags` can find it.

## 2. Install `bml` R package

### Option A: Install from CRAN (Recommended)

The stable release version is available on CRAN:

``` r

install.packages("bml")
```

### Option B: Install development version from GitHub

For the latest development features, install from GitHub using
[remotes](https://remotes.r-lib.org):

``` r

install.packages("remotes")
remotes::install_github("benrosche/bml")
```

If you want to build vignettes locally during installation:

``` r

remotes::install_github("benrosche/bml", build_vignettes = TRUE)
```

## Troubleshooting

- `Error: (converted from warning) ...` while installing.
  [remotes](https://remotes.r-lib.org) treats an installation warning as
  an error on your machine. You can opt out:

  ``` r

  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
  remotes::install_github("benrosche/bml")
  ```

- JAGS not found / `rjags` fails to load: Make sure JAGS is installed
  (not just the R packages), then restart R. On Windows, confirm you
  installed the 64-bit JAGS build if you are using 64-bit R.
