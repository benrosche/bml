# summary() method for an bml object

summary() method for an bml object

## Usage

``` r
# S3 method for class 'bml'
summary(bml, r = 3)
```

## Value

Returns a table with regression results.

## Author

Benjamin Rosche \<benrosche@nyu.edu\>

## Examples

``` r
data(coalgov)
m1 <- bml(Surv(govdur, earlyterm, govmaxdur) ~ 1 + mm(id(pid, gid), mmc(fdep), mmw(w ~ 1/n, constraint=1)) + majority + hm(id=cid, name=cname, type=RE, showFE=F),
          family="Weibull", monitor=T, data=coalgov)
#> Error in c(ids, vars, l1, l3) %<-% dissectFormula(formula, family, data): could not find function "%<-%"
summary(m1)
#> Error: object 'm1' not found
```
