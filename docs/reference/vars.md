# Specify covariates for multiple-membership or hierarchical models

Helper function used within
[`mm`](https://benrosche.github.io/bml/reference/mm.md) and
[`hm`](https://benrosche.github.io/bml/reference/hm.md) to specify which
variables should be included at each level of the model. Supports both
free variables (with coefficients to be estimated) and fixed variables
(with coefficients held constant using
[`fix`](https://benrosche.github.io/bml/reference/fix.md)).

## Usage

``` r
vars(...)
```

## Arguments

- ...:

  Unquoted variable names from your data, combined using `+`
  (formula-style). Supports:

  - Simple variables: `vars(x + y)`

  - Interactions: `vars(x * y)` or `vars(x:y)`

  - Transformations: `vars(I(x^2))` or `vars(I(x + y))`

  - Fixed coefficients: `vars(fix(x, 1.0) + y)`

  Note: Numeric literals like `1`, `0`, or `-1` are ignored (no
  intercept support in mm/hm blocks).

## Value

A `bml_vars` object containing:

- `formula`: Formula object for use with
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)

- `free`: Character vector of base variable names

- `fixed`: List of variables with fixed coefficients (if any)

Returns `NULL` if no variables are specified.

## See also

[`fix`](https://benrosche.github.io/bml/reference/fix.md),
[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`hm`](https://benrosche.github.io/bml/reference/hm.md)

## Examples

``` r
# Simple variable specification (formula-style with +)
vars(income + education)
#> $formula
#> ~0 + income + education
#> <environment: 0x0000022d917c07b0>
#> 
#> $free
#> [1] "income"    "education"
#> 
#> $fixed
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_vars"

# Single variable
vars(income)
#> $formula
#> ~0 + income
#> <environment: 0x0000022d918677b0>
#> 
#> $free
#> [1] "income"
#> 
#> $fixed
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_vars"

# Interactions
vars(income * education)  # expands to income + education + income:education
#> $formula
#> ~0 + income * education
#> <environment: 0x0000022d918b8b30>
#> 
#> $free
#> [1] "income"    "education"
#> 
#> $fixed
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_vars"
vars(income:education)    # interaction only
#> $formula
#> ~0 + income:education
#> <environment: 0x0000022d91982000>
#> 
#> $free
#> [1] "income"    "education"
#> 
#> $fixed
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_vars"

# Transformations
vars(I(income^2))         # squared term
#> $formula
#> ~0 + I(income^2)
#> <environment: 0x0000022d91b66118>
#> 
#> $free
#> [1] "income"
#> 
#> $fixed
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_vars"
vars(income + I(income^2)) # linear and squared
#> $formula
#> ~0 + income + I(income^2)
#> <environment: 0x0000022d91cd7bd8>
#> 
#> $free
#> [1] "income"
#> 
#> $fixed
#> NULL
#> 
#> attr(,"class")
#> [1] "bml_vars"

# Mix free and fixed variables
vars(fix(exposure, 1.0) + income + education)
#> $formula
#> ~0 + income + education
#> <environment: 0x0000022d91d8e238>
#> 
#> $free
#> [1] "income"    "education"
#> 
#> $fixed
#> $fixed[[1]]
#> $fixed[[1]]$var
#> [1] "exposure"
#> 
#> $fixed[[1]]$value
#> [1] 1
#> 
#> 
#> 
#> attr(,"class")
#> [1] "bml_vars"

# Use in mm() specification
if (FALSE) { # \dontrun{
mm(
  id = id(pid, gid),
  vars = vars(rile + ipd),
  fn = fn(w ~ 1/n),
  RE = FALSE
)
} # }
```
