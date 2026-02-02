# Specify identifier variables for multiple-membership and hierarchical structures

Helper function used within
[`mm`](https://benrosche.github.io/bml/reference/mm.md) and
[`hm`](https://benrosche.github.io/bml/reference/hm.md) to specify the
identifier variables that define memberships and nesting structures. In
multiple-membership models, `id()` links member-level units (e.g., party
IDs) to group-level units (e.g., government IDs). In hierarchical
models, `id()` specifies the nesting-level identifier (e.g., country
ID).

## Usage

``` r
id(...)
```

## Arguments

- ...:

  Unquoted variable names from your data:

  - For [`mm()`](https://benrosche.github.io/bml/reference/mm.md): Two
    identifiers `id(mmid, mainid)` where `mmid` identifies member-level
    units and `mainid` identifies group-level units

  - For [`hm()`](https://benrosche.github.io/bml/reference/hm.md): One
    identifier `id(hmid)` where `hmid` identifies nesting-level units

## Value

A `bml_id` object containing the variable names as character strings.

## See also

[`mm`](https://benrosche.github.io/bml/reference/mm.md),
[`hm`](https://benrosche.github.io/bml/reference/hm.md),
[`bml`](https://benrosche.github.io/bml/reference/bml.md)

## Examples

``` r
# Multiple-membership: parties (pid) within governments (gid)
id(pid, gid)
#> [1] "pid" "gid"
#> attr(,"class")
#> [1] "bml_id"

# Hierarchical: governments within countries
id(cid)
#> [1] "cid"
#> attr(,"class")
#> [1] "bml_id"
```
