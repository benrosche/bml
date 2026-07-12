# Reallocation of membership weights toward high-attribute members

The headline temporal micro-macro term: did the group shift weight mass
toward members with high \\x\\? Computes \\\sum_k \Delta w_k x_k\\ over
the union roster (using each member's observed attribute) plus the
reallocation *intensity* \\\sum_k \bar w_k (\Delta w_k)^2\\-style spread
measures.

## Usage

``` r
bml_realloc(data, group, member, time, weight = NULL, x, t1, t2)
```

## Arguments

- data:

  Member-level long data covering both snapshots.

- group, member, time:

  Column names (strings): group id, member id, and the time variable.

- weight:

  Column name (string) of the raw weight, or `NULL` for equal weights.

- x:

  Column name (string) of the member attribute.

- t1, t2:

  The two time values to compare.

## Value

A data frame with one row per group: `realloc_x` (\\\sum_k \Delta w_k
x_k\\; reallocation toward high-`x` members) and `realloc_intensity`
(\\\sum_k (\Delta w_k)^2\\; how much reshuffling happened, regardless of
direction).

## See also

[`bml_delta`](https://benrosche.github.io/bml/reference/bml_delta.md),
[`bml_turnover`](https://benrosche.github.io/bml/reference/bml_turnover.md)
