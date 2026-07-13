# Compare fitted bml models in a publication-style table

Builds a side-by-side regression table from one or more fitted
[`bml`](https://benrosche.github.io/bml/reference/bml.md) models: one
row per term, one column per model, with cells formatted as
`estimate (conf.low, conf.high)`. Goodness-of-fit rows (e.g. `N` and
`DIC`) are appended at the bottom. The result is a `data.frame` of
character cells, ready to pass to
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html),
[`gt::gt()`](https://gt.rstudio.com/reference/gt.html), etc.

For the underlying numeric values (e.g. to build custom tables or plots)
use the broom methods
[`tidy.bml`](https://benrosche.github.io/bml/reference/tidy.bml.md) and
[`glance.bml`](https://benrosche.github.io/bml/reference/glance.bml.md)
directly; they also make `modelsummary::modelsummary(list(m1, m2))` work
on bml models.

## Usage

``` r
bmlCompare(
  ...,
  component = "all",
  terms = NULL,
  labels = NULL,
  stats = c("N", "DIC"),
  digits = 3,
  conf.level = 0.95
)
```

## Arguments

- ...:

  Fitted `bml` model objects, optionally named (e.g.,
  `bmlCompare(base = m1, weighted = m2)`). The names become the column
  headers. Unnamed models are labeled with the expressions they were
  passed as (e.g., `"m1"`).

- component:

  Which parameters to include; passed to
  [`tidy.bml`](https://benrosche.github.io/bml/reference/tidy.bml.md).
  One of `"all"` (default), `"fixed"`, `"random"`, or `"weights"`.

- terms:

  Optional character vector of term names selecting which rows to show
  and in what order (matched against the `term` labels from
  [`tidy.bml`](https://benrosche.github.io/bml/reference/tidy.bml.md),
  e.g. `"smoking_alter (mm.1)"`). Terms not found are skipped with a
  warning. Default `NULL`: all terms in `component`, in their natural
  order.

- labels:

  Optional character vector renaming the rows for display. Must be the
  same length as the selected `terms` (so supply `terms` when using
  `labels`). Default `NULL`: keep the original term labels.

- stats:

  Goodness-of-fit rows to append, in order. Any of `"N"` (main-level
  units), `"n.members"` (member-level units), and `"DIC"`. Default:
  `c("N", "DIC")`. Use `character(0)` for none.

- digits:

  Number of decimal places for the estimate and interval bounds.
  Default: 3.

- conf.level:

  Width of the equal-tailed credible interval. Default: 0.95. Other
  levels require the models to be fitted with `monitor = TRUE` (see
  [`tidy.bml`](https://benrosche.github.io/bml/reference/tidy.bml.md)).

## Value

A `data.frame` whose first column `Term` holds the term labels (and the
requested `stats` labels), followed by one character column per model
containing `estimate (conf.low, conf.high)` cells.

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`tidy.bml`](https://benrosche.github.io/bml/reference/tidy.bml.md),
[`glance.bml`](https://benrosche.github.io/bml/reference/glance.bml.md),
[`coefPlot`](https://benrosche.github.io/bml/reference/coefPlot.md)

## Author

Benjamin Rosche <benrosche@nyu.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
data(coalgov)

m1 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("sum"), RE = TRUE),
  family = weibull(),
  data = coalgov
)
m2 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ pseat), fn = fn("sum"), RE = TRUE),
  family = weibull(),
  data = coalgov
)

bmlCompare("Equal" = m1, "Seat share" = m2)

# Select and relabel specific rows
bmlCompare(
  "Equal" = m1, "Seat share" = m2,
  terms  = c("majority", "cohesion (mm.1)"),
  labels = c("Majority", "Cohesion")
) |> knitr::kable()
} # }
```
