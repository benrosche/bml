# Define a multiple-membership block (the micro-to-macro aggregation)

Specifies a multiple-membership level where group-level units (e.g.,
occupations) are composed of member-level units (e.g., tasks). Each
block aggregates the weighted member records into one or more
*group-level features* whose coefficients are estimated by the main
model. The components map onto the framework's notation
\\\theta^{micro,f}(M\_{it})\\:

- `id = id(mmid, mainid)`: the membership structure \\S\_{it}\\

- `vars = vars(x)`: the member attributes \\x\_{kt}\\

- `w = w(~ ...)`: the weights \\w\_{ikt}\\ (who counts)

- `fn = fn("...")`: the aggregation function \\f\\ (how contributions
  combine); `fn("sum")` is the additive case, other types are emergent
  features (variance, concentration, thresholds, ...)

- `RE = re(...)`: member random effects, aggregated by the weights
  (\\\sum_k w\_{ikt} u\_{0k}\\, plus random slopes via `re(1 + x)`)

- `name`: block name for referencing the feature in main-formula
  interactions (auto-generated when omitted)

## Usage

``` r
mm(
  id,
  vars = NULL,
  w = NULL,
  fn = NULL,
  RE = NULL,
  FE = NULL,
  name = NULL,
  ...
)
```

## Arguments

- id:

  An [`id`](https://benrosche.github.io/bml/reference/id.md) object:
  `id(mmid, mainid)`.

- vars:

  A [`vars`](https://benrosche.github.io/bml/reference/vars.md) object
  with member attributes, or `NULL` for weights-only blocks (RE-only, or
  `fn("hhi")`-family).

- w:

  A [`w`](https://benrosche.github.io/bml/reference/w.md) object
  specifying the weights. Default: `w(~ 1/n, scale = TRUE)` (equal
  weights).

- fn:

  A [`fn`](https://benrosche.github.io/bml/reference/fn.md) object
  selecting the aggregation function. Default: `fn("sum")` (the additive
  weighted mean).

- RE:

  Random effects: `TRUE` (shorthand for `re(1)`), a
  [`re`](https://benrosche.github.io/bml/reference/re.md) object, or
  `NULL`/`FALSE` for none. Only available with `fn("sum")`: member
  random effects do not compose with dispersion-type features.

- FE:

  Fixed effects: a
  [`fe`](https://benrosche.github.io/bml/reference/fe.md) object, or
  `NULL`. Mutually exclusive with `RE`. With many members this is weakly
  identified — prefer `RE` (partial pooling) at the member level.

- name:

  Optional block name (unquoted or string) used to reference the block's
  feature in main-formula interactions (e.g. `Ax:education`).
  Auto-generated from the feature when omitted (e.g. `A_x`, `V_x`). Must
  not collide with a data column.

- ...:

  Not used; catches removed arguments (`c =`, and `ar =`, which moved
  into the effects grammar: `RE = re(1, ar = TRUE)`) with a migration
  message.

## Value

A `bml_mm` object.

## Details

**One mechanism.** Every block emits group-level feature(s) whose
coefficients are main-model coefficients (class `"b"`, labeled by the
feature name, e.g. `A_x`). `fn("sum")` is the additive case \\\theta =
\beta A_x + \sum_k w_k u\_{0k}\\: the feature \\A_x\\ plus the optional
`RE`. There is no separate "sum mode".

**Multiple blocks** can be combined with `+` to stack features (mean +
variance + concentration, ...). `RE` can be specified for one block per
member-id group.

## References

Rosche, B. (2026). A Multilevel Model for Coalition Governments:
Uncovering Party-Level Dependencies Within and Between Governments.
*Political Analysis*.

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`id`](https://benrosche.github.io/bml/reference/id.md),
[`vars`](https://benrosche.github.io/bml/reference/vars.md),
[`w`](https://benrosche.github.io/bml/reference/w.md),
[`fn`](https://benrosche.github.io/bml/reference/fn.md),
[`re`](https://benrosche.github.io/bml/reference/re.md),
[`fe`](https://benrosche.github.io/bml/reference/fe.md),
[`hm`](https://benrosche.github.io/bml/reference/hm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Additive aggregation (weighted mean effect)
mm(id = id(task, occ), vars = vars(x), w = w(~ importance, scale = TRUE))

# With member random intercepts
mm(id = id(task, occ), vars = vars(x), w = w(~ 1/n), RE = TRUE)

# Random intercept + slope (residual effect heterogeneity)
mm(id = id(task, occ), vars = vars(x), w = w(~ 1/n), RE = re(1 + x))

# Emergent features
mm(id = id(task, occ), vars = vars(x), w = w(~ 1/n), fn = fn("var"))
mm(id = id(task, occ), w = w(~ importance), fn = fn("hhi"))
mm(id = id(task, occ), vars = vars(x), fn = fn("threshold", c = est(), kappa = 10))
mm(id = id(task, occ), vars = vars(x), fn = fn("smax", kappa = est()))

# Named block, referenced in a cross-level interaction
bml(Y ~ education + Ax:education +
      mm(name = Ax, id = id(task, occ), vars = vars(x), w = w(~ 1/n)),
    data = dat)
} # }
```
