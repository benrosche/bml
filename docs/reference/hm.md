# Hierarchical membership specification

Specifies a hierarchical (level-3) structure including variables and
effect type.

## Usage

``` r
hm(id, vars = NULL, name = NULL, type = "RE", showFE = FALSE)
```

## Arguments

- id:

  An id() object specifying the level-3 identifier: id(l3id)

- vars:

  A vars() object specifying level-3 variables, or NULL

- name:

  Unquoted variable name for level-3 labels (optional)

- type:

  Character; "RE" for random effects (default) or "FE" for fixed effects

- showFE:

  Logical; if TRUE and type = "FE", report the fixed effects (default:
  FALSE)

## Value

A bml_hm object containing the hierarchical membership specification
