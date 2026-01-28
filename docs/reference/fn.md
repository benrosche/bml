# Specify weight function

Helper function to specify the weight function in mm() objects.

## Usage

``` r
fn(w = w ~ 1/n, c = TRUE, ar = FALSE)
```

## Arguments

- w:

  A formula specifying the weight function (default: w ~ 1/n)

- c:

  Logical; if TRUE (default), weights are constrained to sum to 1

- ar:

  Logical; if TRUE, autoregressive random effects are used (default:
  FALSE)

## Value

A bml_fn object containing the weight function specification
