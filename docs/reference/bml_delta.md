# Temporal composition change: matched decomposition of a feature difference

Computes the change in a group's weighted mean attribute between two
snapshots, \\\Delta A_x = A\_{x,t_2} - A\_{x,t_1}\\, together with its
matched (shift-share / Foster-Haltiwanger-Krizan style) decomposition
over membership spells: \$\$\Delta A_x = \underbrace{\sum_k \Delta w_k
\bar x_k}\_{reallocation} + \underbrace{\sum_k \bar w_k \Delta
x_k}\_{attribute\\ change} + entry/exit\\ terms\$\$ The components sum
to \\\Delta A_x\\ exactly (adding-up identity). The returned group-level
covariates can be entered into a
[`bml`](https://benrosche.github.io/bml/reference/bml.md) main formula.

## Usage

``` r
bml_delta(data, group, member, time, weight = NULL, x, t1, t2)
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

A data frame with one row per group: `dA` (total change), `realloc`
(reallocation among stayers), `attr_change` (attribute change among
stayers), `entry_exit` (net entry/exit contribution).

## Details

Weights are normalized within each group and snapshot, so \\\Delta w\\
sums to zero over the union roster. For stayers, \\\bar x\\ and \\\bar
w\\ are the across-snapshot means; entering/exiting members contribute
their full \\w x\\ term to `entry_exit`.

## See also

[`bml_realloc`](https://benrosche.github.io/bml/reference/bml_realloc.md),
[`bml_turnover`](https://benrosche.github.io/bml/reference/bml_turnover.md)
