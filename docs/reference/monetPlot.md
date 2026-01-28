# Create monetPlots for your bml results.

The function **monetPlot** creates a density plot of the posterior
distribution of your model parameters and the traceplot that led to this
density.

## Usage

``` r
monetPlot(bml, parameter, label = NULL, r = 2, yaxis = T)
```

## Arguments

- bml:

  A bml object. bml has to be run with monitor=T

- parameter:

  A string with the parameter name. The internal name has to be used,
  which are the rownames in the bml reg.table output.

- label:

  String to describe the parameter on the graph's x-axis. Optional. If
  not specified, the internal parameter name is used.

- r:

  Specify number of decimal places. Default equals 3.

- yaxis:

  Logical. If FALSE, the y-axis title is omitted.

## Value

Returns a plot. The solid vertical is at 0 and the dashed vertical line
is the mode of the posterior distributions.

## Author

Benjamin Rosche \<benrosche@nyu.edu\>

## Examples

``` r
data(coalgov)
m1 <- bml(Surv(govdur, earlyterm) ~ 1 + majority + mm(id(pid, gid), mmc(fdep), mmw(w ~ 1/n, constraint=T)) + hm(id=cid, name=cname, type=RE, showFE=F),
          family="Weibull", monitor=T, data=coalgov)
#> Error in c(ids, vars, l1, l3) %<-% dissectFormula(formula, family, data): could not find function "%<-%"
monetPlot(m1, parameter="b.l1")
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: tidyr
#> Error: object 'm1' not found
```
