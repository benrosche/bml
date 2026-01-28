# Specify variables

Helper function to specify variables in mm() and hm() objects.

## Usage

``` r
vars(...)
```

## Arguments

- ...:

  Unquoted variable names, can use + syntax (e.g., vars(x1 + x2) or
  vars(x1, x2))

## Value

A bml_vars object containing the variable names, or NULL if empty
