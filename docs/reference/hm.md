# Define a hierarchical nesting structure

Specifies a hierarchical (nesting) level in the model where group-level
units are nested within higher-level entities. Unlike
multiple-membership structures, each group belongs to exactly one
nesting-level unit. Can model either random effects or fixed effects at
the nesting level.

## Usage

``` r
hm(id, vars = NULL, name = NULL, type = "RE", showFE = FALSE, ar = FALSE)
```

## Arguments

- id:

  An [`id`](https://benrosche.github.io/bml/reference/id.md) object
  specifying the nesting-level identifier: `id(hmid)` where `hmid`
  identifies the higher-level units (e.g., countries, regions).

- vars:

  A [`vars`](https://benrosche.github.io/bml/reference/vars.md) object
  specifying nesting-level covariates, or `NULL` for intercept-only
  effects. Supports interactions (`*`, `:`) and transformations
  ([`I()`](https://rdrr.io/r/base/AsIs.html)).

- name:

  Unquoted variable name for nesting-level labels (optional). If
  provided, these labels will be displayed in model output for fixed
  effects.

- type:

  Character; either `"RE"` for random effects (default) or `"FE"` for
  fixed effects at the nesting level.

- showFE:

  Logical; if `TRUE` and `type = "FE"`, fixed effect estimates for each
  nesting-level unit are included in output. Default: `FALSE`.

- ar:

  Logical; if `TRUE`, random effects evolve autoregressively across
  participations at the nesting level. Requires sequential participation
  indicators in the data. Default: `FALSE`.

## Value

A `bml_hm` object containing the hierarchical specification.

## Details

**Hierarchical vs. Multiple-Membership:**

Hierarchical structures (`hm`) model strict nesting: each group belongs
to exactly one higher-level unit. Use
[`mm`](https://benrosche.github.io/bml/reference/mm.md) when groups can
have memberships in multiple units.

**Random vs. Fixed Effects:**

- **Random effects** (`type = "RE"`): Nesting-level units are treated as
  a random sample from a population. Best when you have many units and
  want to generalize.

- **Fixed effects** (`type = "FE"`): Each unit gets its own parameter.
  Best when you have few units or want to estimate unit-specific
  effects.

**Cross-Classification:**

Multiple `hm()` blocks create cross-classified models where groups are
simultaneously nested within multiple non-nested hierarchies (e.g.,
schools within both neighborhoods and districts).

## References

Goldstein, H. (2011). *Multilevel Statistical Models* (4th ed.). Wiley.

Rabe-Hesketh, S., & Skrondal, A. (2012). *Multilevel and Longitudinal
Modeling Using Stata* (3rd ed.). Stata Press.

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`id`](https://benrosche.github.io/bml/reference/id.md),
[`vars`](https://benrosche.github.io/bml/reference/vars.md)

## Examples

``` r
# Random effects with covariates
hm(
  id = id(cid),
  vars = vars(gdp + democracy),
  name = cname,
  type = "RE"
)
#> $id
#> [1] "cid"
#> 
#> $vars
#> $formula
#> ~0 + gdp + democracy
#> <environment: 0x0000022d90960c08>
#> 
#> $free
#> [1] "gdp"       "democracy"
#> 
#> $fixed
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_vars"
#> 
#> $name
#> [1] "cname"
#> 
#> $type
#> [1] "RE"
#> 
#> $showFE
#> [1] FALSE
#> 
#> $ar
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "bml_hm"

# Random intercepts only
hm(
  id = id(cid),
  vars = NULL,
  type = "RE"
)
#> $id
#> [1] "cid"
#> 
#> $vars
#> NULL
#> 
#> $name
#> NULL
#> 
#> $type
#> [1] "RE"
#> 
#> $showFE
#> [1] FALSE
#> 
#> $ar
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "bml_hm"

# Fixed effects
hm(
  id = id(cid),
  vars = NULL,
  name = cname,
  type = "FE",
  showFE = TRUE  # Show estimates for each country
)
#> $id
#> [1] "cid"
#> 
#> $vars
#> NULL
#> 
#> $name
#> [1] "cname"
#> 
#> $type
#> [1] "FE"
#> 
#> $showFE
#> [1] TRUE
#> 
#> $ar
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "bml_hm"

# Autoregressive random effects
hm(
  id = id(cid),
  vars = NULL,
  type = "RE",
  ar = TRUE  # Effects evolve over time
)
#> $id
#> [1] "cid"
#> 
#> $vars
#> NULL
#> 
#> $name
#> NULL
#> 
#> $type
#> [1] "RE"
#> 
#> $showFE
#> [1] FALSE
#> 
#> $ar
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "bml_hm"
```
