# Draws from the posterior predictive distribution

Returns draws of replicated outcomes `y_rep` (the monitored `pred`
node). `newdata` is not supported yet: multiple-membership prediction
for new units requires the new units' membership map, which is future
work.

## Usage

``` r
posterior_predict(object, ...)

# S3 method for class 'bml'
posterior_predict(object, ndraws = NULL, ...)
```

## Arguments

- object:

  A fitted `bml` model with `monitor = "predictive"`.

- ...:

  Unused.

- ndraws:

  Optional number of draws to return (subsampled without replacement).

## Value

A draws x observations matrix.
