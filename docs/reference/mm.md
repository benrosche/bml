# Multiple membership specification

Specifies a multiple membership structure including variables, weight
function, and random effects.

## Usage

``` r
mm(id, vars = NULL, fn = fn(), RE = NULL)
```

## Arguments

- id:

  An id() object specifying level-1 and level-2 identifiers: id(l1id,
  l2id)

- vars:

  A vars() object specifying level-1 variables to aggregate, or NULL for
  RE-only

- fn:

  A fn() object specifying the weight function (default: fn(w ~ 1/n, c =
  TRUE, ar = FALSE))

- RE:

  Logical; if TRUE, include random effects. Automatically TRUE if vars
  is NULL.

## Value

A bml_mm object containing the multiple membership specification
