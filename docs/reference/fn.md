# The aggregation function of an mm() block

Selects the function `f` that reduces the weighted member records to a
group-level feature — the `f` in the framework's
\\\theta^{micro,f}(M\_{it})\\. `fn("sum")` (the default) is the additive
weighted mean; the other named types are "emergent" features of the
whole member set. `fn()` also accepts a one-sided formula written in a
restricted expression DSL, so users can define their own reductions.

**Named types:**


    fn("sum")                          # A_x  = sum_k w_k x_k          (default)
    fn("var")                          # V_x  = sum_k w_k (x_k - A_x)^2
    fn("var", moment = 3)              # higher central moments
    fn("hhi")                          # C    = sum_k w_k^2  (weights only)
    fn("effn")                         # 1 / C
    fn("entropy")                      # -sum_k w_k log w_k
    fn("threshold", c = 0.7, kappa = 10)  # T(c) = sum_k w_k ilogit(kappa (x_k - c))
    fn("threshold", c = est())         # estimate the cutpoint
    fn("smax", kappa = 5)              # (1/kappa) log sum_k w_k exp(kappa x_k)
    fn("smax", kappa = est())          # estimate the aggregation function itself:
                                       #   kappa<0 -> min, kappa->0 -> mean, kappa>0 -> max
    fn("gmean", p = est())             # power/CES mean (sum_k w_k x_k^p)^(1/p); x > 0
    fn("cov")                          # C_xz = sum_k w_k (x_k - A_x)(z_k - A_z); two attributes

**Expression DSL:**
[`w()`](https://benrosche.github.io/bml/reference/w.md) normalizes
weights to a probability measure over the group's members; `E(e)` is the
expectation \\\sum_k w_k e_k\\ under that measure. Member-level
quantities are the
[`vars()`](https://benrosche.github.io/bml/reference/vars.md) attributes
and the reserved symbol `w`; anything wrapped in `E()` is a group
scalar. Whitelisted operations: `+ - * / ^`, `exp`, `log`, `ilogit`,
`pow`. Any symbol that is not a data column or reserved word is a *free
parameter* with a default prior (one rule, shared with
[`w`](https://benrosche.github.io/bml/reference/w.md)); the build
messages the detected parameters.


    fn(~ E((x - E(x))^2))                       # variance, written out
    fn(~ E(ilogit(kappa * (x - c))))            # threshold with free c, kappa
    fn(~ E((x - E(x)) * (z - E(z))))            # covariance

**Identification:** internal parameters reach the outcome only through
`b * feature`; if the feature's coefficient is ~0 they are unidentified.
Dispersion functions need real within-group spread; threshold/tail
functions need mass of `x` in the region they read. A free parameter
that multiplies an attribute inside a nonlinear `E()` is confounded with
the feature's coefficient — treat inner parameters as shape parameters
(cutpoints, sharpness), not slopes.

## Usage

``` r
fn(type = "sum", ...)
```

## Arguments

- type:

  A type string (one of `"sum"`, `"var"`, `"hhi"`, `"effn"`,
  `"entropy"`, `"threshold"`, `"smax"`, `"gmean"`, `"cov"`) or a
  one-sided formula (`~ E(...)`).

- ...:

  Shape parameters for the named types: `moment` (var), `c` and `kappa`
  (threshold), `kappa` (smax), `p` (gmean). Each is a number (fixed) or
  [`est()`](https://benrosche.github.io/bml/reference/est.md)
  (estimated).

## Value

A `bml_fn` object.

## See also

[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`w`](https://benrosche.github.io/bml/reference/w.md),
[`est`](https://benrosche.github.io/bml/reference/est.md),
[`prior`](https://benrosche.github.io/bml/reference/prior.md)
