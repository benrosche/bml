# Define a hierarchical nesting structure

Specifies a hierarchical (nesting) level: each group-level unit belongs
to exactly one nesting-level unit (e.g., governments within countries).
`hm()` carries only the effect *structure* — random effects
(`RE = re(...)`, partial pooling) or fixed effects (`FE = fe(...)`, no
pooling). Fixed effects of unit-level *covariates* belong in the main
formula, following the lme4/brms convention:


    Y ~ gdp + hm(id = id(cid), RE = re(1 + gdp))   # like y ~ gdp + (1 + gdp | cid)

## Usage

``` r
hm(id, RE = NULL, FE = NULL, labels = NULL, ...)
```

## Arguments

- id:

  An [`id`](https://benrosche.github.io/bml/reference/id.md) object:
  `id(hmid)`.

- RE:

  Random effects: `TRUE` (shorthand for `re(1)`), an
  [`re`](https://benrosche.github.io/bml/reference/re.md) object, or
  `NULL`. Slope variables are main-level columns varying within the
  nesting units. Default when neither `RE` nor `FE` is given: `re(1)`.

- FE:

  Fixed effects: a
  [`fe`](https://benrosche.github.io/bml/reference/fe.md) object (unit
  dummies via `fe(1)`; unit-specific slopes via `fe(1 + x)`; report
  per-unit estimates with `fe(1, showFE = TRUE)`). Mutually exclusive
  with `RE`.

- labels:

  Unquoted variable name holding display labels for the nesting units
  (used when reporting per-unit fixed effects).

- ...:

  Not used; catches the removed `type =`, `vars =`, `name =`,
  `showFE =`, and `ar =` arguments (the last moved into the effects
  grammar: `RE = re(1, ar = TRUE)`) with a migration message.

## Value

A `bml_hm` object.

## Details

Cross-classified structures are modeled by including multiple `hm()`
blocks. The old `type = "RE"/"FE"` and `vars =` arguments were removed:
use `RE = re(...)` / `FE = fe(...)`, and move unit-level covariates into
the main formula (their coefficients are ordinary main-model
coefficients).

## References

Goldstein, H. (2011). *Multilevel Statistical Models* (4th ed.). Wiley.

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`id`](https://benrosche.github.io/bml/reference/id.md),
[`re`](https://benrosche.github.io/bml/reference/re.md),
[`fe`](https://benrosche.github.io/bml/reference/fe.md)

## Examples

``` r
if (FALSE) { # \dontrun{
hm(id = id(cid))                          # random intercepts (default)
hm(id = id(cid), RE = re(1 + x))          # random intercept + slope on x
hm(id = id(cid), FE = fe(1))              # country dummies (no pooling)
hm(id = id(cid), FE = fe(1, showFE = TRUE), labels = cname)
} # }
```
