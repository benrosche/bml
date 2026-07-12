# Membership turnover between two snapshots

Set-level composition dynamics: how much of the group's membership mass
entered or exited between the two snapshots.

## Usage

``` r
bml_turnover(data, group, member, time, weight = NULL, t1, t2, x = NULL)
```

## Arguments

- data:

  Member-level long data covering both snapshots.

- group, member, time:

  Column names (strings): group id, member id, and the time variable.

- weight:

  Column name (string) of the raw weight, or `NULL` for equal weights.

- t1, t2:

  The two time values to compare.

- x:

  Not used (turnover is a property of the roster and weights alone).

## Value

A data frame with one row per group: `entry_mass` (weight mass of
entrants at `t2`), `exit_mass` (weight mass of exiters at `t1`),
`jaccard_dist` (1 - \|intersection\| / \|union\| of the member sets),
and `realloc_l1` (\\\frac{1}{2}\sum_k \|\Delta w_k\|\\, the total
variation distance between the two weight measures).

## See also

[`bml_delta`](https://benrosche.github.io/bml/reference/bml_delta.md),
[`bml_realloc`](https://benrosche.github.io/bml/reference/bml_realloc.md)
