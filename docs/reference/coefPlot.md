# Coefficient plot for fitted bml models

Creates a coefficient ("dot-and-whisker") plot from one or more fitted
[`bml`](https://benrosche.github.io/bml/reference/bml.md) models:
posterior means as points and credible intervals as horizontal ranges,
with a dashed reference line at zero. When several models are supplied,
their coefficients are placed next to each other for each term so the
models can be compared directly.

## Usage

``` r
coefPlot(
  ...,
  parameters = NULL,
  intercept = TRUE,
  conf.level = 0.95,
  component = "fixed"
)
```

## Arguments

- ...:

  One or more fitted `bml` model objects, optionally named (e.g.,
  `coefPlot(base = m1, weighted = m2)`). Names are used in the legend;
  unnamed models are labeled with the expressions they were passed as
  (e.g., `"m1"`).

- parameters:

  Optional character vector of term labels (as shown in
  [`summary()`](https://rdrr.io/r/base/summary.html)) to include, in the
  order they should appear (top to bottom). Default: all terms of the
  selected `component`.

- intercept:

  Logical. Include the intercept? Default: `TRUE`.

- conf.level:

  Width of the equal-tailed credible interval. Default: 0.95. Other
  levels require the models to be fitted with `monitor = TRUE` (see
  [`tidy.bml`](https://benrosche.github.io/bml/reference/tidy.bml.md)).

- component:

  Which parameters to plot; passed to
  [`tidy.bml`](https://benrosche.github.io/bml/reference/tidy.bml.md).
  Default: `"fixed"` (regression coefficients).

## Value

A `ggplot` object that can be further customized with ggplot2 layers or
saved with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Details

The plot uses a colorblind-friendly palette (Okabe-Ito) to distinguish
models and intentionally substantial line widths and point sizes for
readability in presentations and publications. All defaults can be
overridden by adding ggplot2 layers to the returned object, e.g.
`coefPlot(m1) + ggplot2::scale_colour_brewer(palette = "Dark2")`.

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md),
[`tidy.bml`](https://benrosche.github.io/bml/reference/tidy.bml.md),
[`bmlCompare`](https://benrosche.github.io/bml/reference/bmlCompare.md),
[`monetPlot`](https://benrosche.github.io/bml/reference/monetPlot.md)

## Author

Benjamin Rosche <benrosche@nyu.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
data(coalgov)

m1 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE),
  family = weibull(),
  data = coalgov
)

# Single model
coefPlot(m1)

# Without intercept
coefPlot(m1, intercept = FALSE)

# Compare two models side by side
coefPlot(equal = m1, seatshare = m2)

# Save the plot
p <- coefPlot(m1, m2)
ggplot2::ggsave("coefplot.pdf", p, width = 6, height = 4)
} # }
```
