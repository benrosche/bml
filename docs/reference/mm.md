# Define a multiple-membership structure

Specifies a multiple-membership level in the model where group-level
units (e.g., governments) are composed of multiple member-level units
(e.g., political parties). Unlike pure hierarchical nesting, members can
belong to multiple groups, and their contributions are aggregated using
a user-specified weight function.

## Usage

``` r
mm(id, vars = NULL, fn = NULL, RE = NULL, ar = FALSE)
```

## Arguments

- id:

  An [`id`](https://benrosche.github.io/bml/reference/id.md) object
  specifying the member-level and group-level identifiers:
  `id(mmid, mainid)` where `mmid` identifies members and `mainid`
  identifies groups.

- vars:

  A [`vars`](https://benrosche.github.io/bml/reference/vars.md) object
  specifying member-level covariates to aggregate, or `NULL` for random
  effects only. Supports interactions (`*`, `:`) and transformations
  ([`I()`](https://rdrr.io/r/base/AsIs.html)). Variables are weighted
  according to the function specified in `fn`.

- fn:

  A [`fn`](https://benrosche.github.io/bml/reference/fn.md) object
  specifying the weight function (default: `fn(w ~ 1/n, c = TRUE)` for
  equal weights). Note: Weight functions do NOT support interactions or
  [`I()`](https://rdrr.io/r/base/AsIs.html) - pre-create any needed
  transformed variables in your data. See
  [`fn`](https://benrosche.github.io/bml/reference/fn.md) for details.

- RE:

  Logical; if `TRUE`, include random effects for member-level units.
  Automatically set to `TRUE` if `vars = NULL` (random effects only).

- ar:

  Logical; if `TRUE`, random effects evolve autoregressively across
  participations. Requires members to have sequential participation
  indicators in the data. Default: `FALSE`.

## Value

A `bml_mm` object containing the multiple-membership specification.

## Details

**Multiple-Membership Models:**

In standard hierarchical models, each observation belongs to exactly one
group. Multiple-membership models relax this assumption, allowing groups
to be composed of multiple members, with flexible weighting of member
contributions.

**Model Structure:**

The contribution from mm block \\k\\ to group \\j\\ is:

\$\$\text{mm}\_{kj} = \sum\_{i \in group_j} w\_{ki}
(\mathbf{x}\_{ki}'\boldsymbol{\beta}\_k + \alpha\_{ki})\$\$

where:

- \\w\_{ki}\\: Weight for member \\i\\ in group \\j\\ (from `fn`)

- \\\mathbf{x}\_{ki}\\: Member-level covariates (from `vars`)

- \\\boldsymbol{\beta}\_k\\: Regression coefficients (estimated)

- \\\alpha\_{ki}\\: Member-level random effect (if `RE = TRUE`)

**Multiple mm() Blocks:**

You can specify multiple `mm()` blocks with different weight functions,
variables, or random effect specifications. This allows modeling
different aggregation mechanisms simultaneously.

## References

Browne, W. J., Goldstein, H., & Rasbash, J. (2001). Multiple membership
multiple classification (MMMC) models. *Statistical Modelling*, 1(2),
103-124.

Fielding, A., & Goldstein, H. (2006). *Cross-classified and multiple
membership structures in multilevel models: An introduction and review*.
Research Report RR791, Department for Education and Skills.

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`id`](https://benrosche.github.io/bml/reference/id.md),
[`vars`](https://benrosche.github.io/bml/reference/vars.md),
[`fn`](https://benrosche.github.io/bml/reference/fn.md),
[`hm`](https://benrosche.github.io/bml/reference/hm.md)

## Examples

``` r
# Equal weights with variables
mm(
  id = id(pid, gid),
  vars = vars(rile + ipd),
  fn = fn(w ~ 1/n, c = TRUE),
  RE = FALSE
)
#> Error in mm(id = id(pid, gid), vars = vars(rile + ipd), fn = fn(w ~ 1/n,     c = TRUE), RE = FALSE): could not find function "mm"

# Random effects only (no variables)
mm(
  id = id(pid, gid),
  vars = NULL,
  fn = fn(w ~ 1/n, c = TRUE),
  RE = TRUE  # Automatically TRUE when vars = NULL
)
#> Error in mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n, c = TRUE),     RE = TRUE): could not find function "mm"

# Flexible weights with parameter
mm(
  id = id(pid, gid),
  vars = vars(party_strength),
  fn = fn(w ~ b0 + b1 * tenure, c = TRUE),
  RE = TRUE
)
#> Error in mm(id = id(pid, gid), vars = vars(party_strength), fn = fn(w ~     b0 + b1 * tenure, c = TRUE), RE = TRUE): could not find function "mm"

# Autoregressive random effects
mm(
  id = id(pid, gid),
  vars = NULL,
  fn = fn(w ~ 1/n, c = TRUE),
  RE = TRUE,
  ar = TRUE  # Random effects evolve over participations
)
#> Error in mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n, c = TRUE),     RE = TRUE, ar = TRUE): could not find function "mm"

# Interactions and transformations in vars
mm(
  id = id(pid, gid),
  vars = vars(rile * ipd),  # Main effects plus interaction
  fn = fn(w ~ 1/n, c = TRUE),
  RE = FALSE
)
#> Error in mm(id = id(pid, gid), vars = vars(rile * ipd), fn = fn(w ~ 1/n,     c = TRUE), RE = FALSE): could not find function "mm"

mm(
  id = id(pid, gid),
  vars = vars(rile + I(rile^2)),  # Quadratic term
  fn = fn(w ~ 1/n, c = TRUE),
  RE = FALSE
)
#> Error in mm(id = id(pid, gid), vars = vars(rile + I(rile^2)), fn = fn(w ~     1/n, c = TRUE), RE = FALSE): could not find function "mm"
```
