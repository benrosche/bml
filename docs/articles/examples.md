# Examples

I present three examples to demonstrate how the MMMM can be applied to
analyze spatial, network, and aggregation problems.

- [Example 1: The effect of air quality on home values](#e1) (spatial
  regression)
- [Example 2: All friends or just your best friend?](#e2) (network
  regression)
- [Example 3: The effect of political parties’ financial dependency on
  the survival coalition governments](#e3) (aggregation regression)

### The effect of air quality on home values

The R file of this example can be found
[here](https://benrosche.github.io/bml/articles/examples_files/bml-examples.Rmd).

This **spatial analysis** example is based on [Harrison & Rubenfield
1978](https://www.sciencedirect.com/science/article/abs/pii/0095069678900062),
who study the effect of air quality on home values, and a re-analysis by
[Bivand 2017](https://openjournals.wu.ac.at/region/paper_107/107.html).

The study employs census tract data from the Boston Standard
Metropolitan Statistical Area in 1970. With tracts containing no housing
units or comprised entirely of institutions excluded, the sample
contains 506 census tracts. Air quality is measured by the concentration
of nitric oxides in the air, which is obtained from a meteorological
model (Transportation and Air Shed Simulation Model). I refer to their
paper for more information on data and operationalization.

Let us load in the data and plot the home values across town:

``` r


library(spData)     # spatial datasets
library(sf)         # read spatial datasets
library(spdep)      # create spatial weights
library(spatialreg) # spatial regression models
library(dplyr)
library(ggplot2)
library(bml)

# Load spatial data from Boston
boston <-
  read_sf(system.file("shapes/boston_tracts.gpkg", package = "spData")) %>%
  select(CMEDV, NOX, CRIM, RM, DIS, AGE, LSTAT, geom) %>%
  st_transform(crs = 5070) %>% # use Albers equal-area conic projection
  mutate(tid = row_number(), lnCMEDV=log(CMEDV), across(c(NOX, CRIM, RM, DIS, AGE, LSTAT), scale)) %>%
  relocate(tid, CMEDV, lnCMEDV)

# Dependent variable:
# lnCMEDV = ln(median home value in $1000)

# Plot median home values across Boston
ggplot(boston, aes(fill = CMEDV)) +
  geom_sf(color = NA) +
  labs(fill = "Median home value") +
  scale_fill_viridis_b() +
  theme(legend.position = "bottom")

# Explanatory variables (standardized):
# NOX     = nitric oxides concentration
# CRIM    = per capita crime
# RM      = avg. number of rooms per dwelling
# DIS     = weighted distance to five Boston employment centers
# AGE     = proportion of units built prior 1940
# LSTAT   = percentage working-class population
```

First, we examine three canonical spatial regression models:

1.  **Residual spatial effect**: $`Y = XB + \lambda Wu + \epsilon`$.
    This model is called the spatial error model because it incorporates
    the residuals of other spatial units into the regression equation of
    the focal unit.
2.  **Exogenous spatial effects**: $`Y = XB + WXB + \epsilon`$. This
    model is called the spatial lag-x model because it incorporates the
    covariates of other spatial units into the regression equation of
    the focal unit.
3.  **Endogenous spatial effect**: $`Y = \rho WY + XB + \epsilon`$. This
    model is called the spatial autoregressive model because it
    incorporates the outcomes of other spatial units into the regression
    equation of the focal unit.

More information on those models and combinations of them can be found
in [Gibbons, Overman & Patacchini
2015](https://www.sciencedirect.com/science/article/pii/B9780444595171000039).

To specify any of these models, we must construct a **spatial weight
matrix** $`W`$, which defines the neighborhood structure for each
location. $`W`$ is an $`N \times N`$ matrix, where $`N`$ represents the
number of neighborhoods. Each element $`w_{ij}`$ quantifies the
relationship between locations $`i`$ and $`j`$, with the convention that
$`w_{ii} = 0`$ along the diagonal.

Neighborhood relationships can be defined in two ways:  
- **Binary contiguity**: $`w_{ij} = 1`$ if $`i`$ and $`j`$ are
neighbors, and 0 otherwise.  
- **Distance-based weighting**: $`w_{ij}`$ is a continuous function of
distance.

Continuous weights are often row-standardized so that the sum of all
weights for a given location $`i`$ equals 1.

The R package `spatialreg` takes the weight matrix as list:

``` r

# Create row-standardized weight matrix
boston_nb <- poly2nb(as_Spatial(boston), row.names = boston$tid) # from polygon list to neighbor list

boston_wmat  <- nb2mat(boston_nb, zero.policy = TRUE) %>% as.matrix() # weight matrix
boston_wlist <- nb2listw(boston_nb, style = "W") # weight matrix as list
```

Now let us estimate the three models:

``` r


# Residual spatial effect (spatial error model):
mod1 <-
  errorsarlm(
    lnCMEDV ~ NOX + CRIM + RM + DIS + AGE,
    data = boston,
    listw = boston_wlist
  )
mod1 %>% summary()

# Exogenous spatial effect (spatial lag-x model):
mod2 <-
  lmSLX(
    lnCMEDV ~ NOX + CRIM + RM + DIS + AGE,
    data = boston,
    listw = boston_wlist,
    Durbin =  ~ NOX + CRIM + RM + DIS + AGE
  )
mod2 %>% summary()

# Endogenous spatial effect (spatial autoregressive model):
mod3 <-
  lagsarlm(
    lnCMEDV ~ NOX + CRIM + RM + DIS + AGE,
    data = boston,
    listw = boston_wlist
  )
mod3 %>% summary()
```

**Results**:

``` r

# 1. The spatial error model estimates that after conditioning on those five covariates, the residual is still spatially correlated, as $\lambda$=`r round(mod1$lambda,2)`.

# 2.  The spatial lag-x model estimates an exogenous spatial effect for each of the considered covariates. The effect of air quality of neighboring locations, for instance, is estimated to be $\beta$=`r round(mod2$coefficients['lag.NOX'],2)`. That is, the home values of a given neighborhood increases if the air quality of surrounding neighborhoods decreases.

# 3.  The spatial autoregessive (SAR) model estimates a endogenous spatial effect of $\rho$=`r round(mod3$rho,2)`. That is, the SAR summarizes the spatial dependency in one coefficient.
```

#### The MMMM for spatial analysis

Let $`y_{i}`$ denote the outcome for location $`i`$.

Using the MMMM, we can model this outcome based on:  
(i) the effects of the location’s own features,
$`x_{i}^{\intercal} \beta + \epsilon_{i}`$  
(ii) the effects of its neighbors’ features,
$`\sum_{j \in n(i)} w_{j} (z_{j}^{\intercal} \gamma + u_{j})`$, where
$`j`$ indexes the neighbors, $`n(i)`$ is the set of neighbors of
location $`i`$, $`z_{j}`$ represents the observed features of neighbor
$`j`$, and $`u_{j}`$ captures the combined influence of unobserved
features.

This model closely resembles a combination of the spatial lag and
spatial error models. The key distinction is that the error term for
each location is decomposed into two components: a random effect for its
role as a focal location and a separate random effect for its role as a
neighbor.

Conceptually, the combination of exogenous spatial effects and spatial
error is intuitive, as a location is influenced by its neighbors’ entire
right-hand side of the regression equation. In contrast, the endogenous
spatial effect model is more challenging to interpret and faces
identification issues when both endogenous and exogenous effects are
included (spatial Durbin model).

To estimate a spatial MMMM, the neighbors of each location must be
included in the dataframe as individual rows:

``` r

# Neighbor list to data.frame
nb2df <- function(nb) {
  return(
    data.frame(
      tid = rep(1:length(nb), sapply(nb, length)),
      tid_nb = unlist(nb)
    )
  )
}

boston_df <-
  nb2df(boston_nb) %>%
  group_by(tid) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  inner_join(boston, by=c("tid")) %>% # own features
  inner_join(                         # neighbor features
    as.data.frame(boston) %>%
      select(-CMEDV,-lnCMEDV, -geom) %>%
      rename_with(~paste0(.,"_nb")),
    by=c("tid_nb")
  )

head(boston_df %>% select(tid, tid_nb, NOX, CRIM, NOX_nb, CRIM_nb))
```

For each tract `tid`, we have one row for each of its neighbors
`tid_nb`. These rows contain the covariates of `tid`, which remain
constant across its neighbors, as well as the covariates of the
neighbors themselves, such as `NOX_nb`, `CRIM_nb`, ….

With this setup, we are now ready to estimate the MMMM:

``` r


# Spatial random effect:
mod.bml1 <-
  bml(
    lnCMEDV ~
      NOX + CRIM + RM + DIS + AGE +
      mm(
        id(tid_nb, tid),
        mmc(),
        mmw(w ~ 1/n, c=T)
      ),
    n.iter = 1000, n.burnin = 100, seed=1, monitor = T,
    data = boston_df
  )

names(mod.bml1)

mod.bml1 %>% summary()

# Spatial fixed effects + spatial random effect:
mod.bml2 <-
  bml(
    lnCMEDV ~
      NOX + CRIM + RM + DIS + AGE +
      mm(
        id(tid_nb, tid),
        mmc(NOX_nb + CRIM_nb + RM_nb + DIS_nb + AGE_nb),
        mmw(w ~ 1/n, c=T)
      ),
    n.iter = 1000, n.burnin = 100, seed=1, monitor = T,
    data = boston_df
  )

mod.bml2 %>% summary()

# Calculate spatial correlation in the residual
getLambda <- function(x) {
  s.l1 <- x$reg.table["sigma.l1", "coefficients"]
  s.l2 <- x$reg.table["sigma.l2", "coefficients"]
  return(s.l1^2/(s.l1^2+s.l2^2))
}

mod.bml1 %>% getLambda()
mod.bml2 %>% getLambda()
```

Let’s take a closer look at
the[`bml()`](https://benrosche.github.io/bml/reference/bml.md) function
by typing [`?bml`](https://benrosche.github.io/bml/reference/bml.md).

**The formula object**

The most important component is the formula object, which in our case
looks like this:
`lnCMEDV ~ NOX + CRIM + RM + DIS + AGE + mm(id(tid_nb, tid), mmc(NOX_nb + CRIM_nb + RM_nb + DIS_nb + AGE_nb), mmw(w ~ 1/n, constraint=1))`

The only difference compared to a
[`lm()`](https://rdrr.io/r/stats/lm.html) formula is the inclusion of
the [`mm()`](https://benrosche.github.io/bml/reference/mm.md) container.
Within this container, we define three sub-containers:  
- `ids()` for the identifiers  
- `mmc()` for the covariates being considered,  
- `mmw()` to endogenize the weight function. Here, `w ~ 1/n` specifies
the row-standardized weight.

The combination of the spatial lag and spatial error models can also be
estimated within the spatial regression framework. Let’s estimate it and
compare the results:

``` r


# Exogenous + residual spatial effect (combination of spatial lag model and spatial error model):
mod2 <-
  errorsarlm(
    lnCMEDV ~ NOX + CRIM + RM + DIS + AGE,
    data = boston,
    listw = boston_wlist,
    Durbin =  ~ NOX + CRIM + RM + DIS + AGE
  )
mod2 %>% summary()

# Plot coefficients next to each other
coefs <-
  mod.bml2$reg.table[,c(1,2,4,5)] %>%
  mutate(model="MMMM") %>%
  add_row(
    data.frame(model="SLM", coefficients=coef(mod2)[-1], confint(mod2, level=0.95)[-1,]) %>%
      rename(lb=X2.5.., ub=X97.5..) %>%
      tibble::rownames_to_column("variable") %>%
      mutate(
        variable=
          case_when(
            startsWith(variable, "lag.") ~ paste0(sub("lag.", "", variable), "_nb"),
            variable=="(Intercept)" ~ "X0",
            TRUE ~ variable))
  ) %>%
  filter(!variable %in% c("X0", "sigma.l1", "sigma.l2", "DIC"))

ggplot(coefs, aes(x=variable, y=coefficients, color=model)) +
  geom_hline(yintercept = 0, color="red") +
  geom_point(position=position_dodge(width=0.3)) +
  geom_pointrange(aes(ymin = lb, ymax = ub), position=position_dodge(width=0.3))+
  labs(title = "Coefficients", x = "Variables", y="", color="Model") +
  theme(legend.position = "bottom") +
  coord_flip()
```

The estimates are similar but not identical due to differences in the
estimation algorithms (maximum likelihood vs. Bayesian MCMC) and slight
differences in model specifications.

The advantage of using `errorsarlm` is its faster estimation and ability
to account for endogenous errors. If this model aligns with your needs,
it is generally best to estimate it using the `spatialreg` package.

**Weight function regression**

The MMMM, however, allows the weights to be modeled as a function of
covariates. Instead of assuming a fixed weighting scheme, we can
estimate whether a neighbor’s influence on a location varies based on
specific covariates.

In this example, I hypothesize that a neighbor’s weight in the overall
neighborhood effect depends on the similarity between the neighbor and
the focal location—an idea inspired by social network theory.

I define *similarity* as the inverse of the average absolute difference
between a neighbor and the focal location across the six considered
covariates. This measure captures dissimilarity, so I label the variable
`DIFF`:

``` r

# Difference between a focal location and its neighbors
boston_df2 <-
  boston_df %>%
  mutate(
    DIFF=1/6*(abs(NOX-NOX_nb)+abs(CRIM-CRIM_nb)+abs(RM-RM_nb)+abs(DIS-DIS_nb)+abs(AGE-AGE_nb)+abs(LSTAT-LSTAT_nb))
  )
```

The weight function is specified within the `mmw()` container, allowing
for any function that produces bounded weights. Logical operators such
as *==, \>, \<* can be used to define conditions that enable aggregation
functions, such as `min` or `max`.

Here, I use the following functional form:

``` math
w = \frac{1}{n^{\exp(-X\beta)}}
```

where $`n`$ is the number of neighbors of location $`i`$, and $`X\beta`$
represents the covariates used to determine the weights.

This formulation has two key advantages:  
1. It ensures weights remain bounded between 0 and 1.  
2. When $`X\beta = \boldsymbol{0}`$, the weights simplify to
$`w = \frac{1}{n}`$, which corresponds to the standard row-standardized
weights.

**The issue of scaling**:

With row-standardized weights, the total sum of weights is
$`\sum_{i=1}^{506} \sum_{j} w_{ij} = \sum_{i=1}^{506} 1 = 506`$.

To ensure this overall sum remains unchanged, we need to carefully
specify constraints. Recall that neighbor effects are aggregated as
$`\sum_{j \in n(i)} w_{j} (z_{j}^{\intercal} \gamma + u_{j})`$.

If the total sum of weights changes, the regression coefficients
$`\gamma`$ will be rescaled. To prevent this, we can apply one of two
constraints:

- **`constraint=1`**: Ensures that the weights of a location’s neighbors
  sum to 1 for each focal location, while allowing them to vary within
  each location.  
- **`constraint=2`**: Ensures that all weights sum to the total number
  of neighborhoods, allowing them to vary both within and across
  locations.

Both constraints identify the model but have different substantive
interpretations. Here, we will use weights that sum to 1 for each focal
location while allowing variation within locations (constraint=1):

``` r

# Spatial fixed effects + spatial random effect with "endogenized" weights:
mod.bml3 <-
  bml(
    lnCMEDV ~
      NOX + CRIM + RM + DIS + AGE +
      mm(
        id(tid_nb, tid),
        mmc(NOX_nb + CRIM_nb + RM_nb + DIS_nb + AGE_nb),
        mmw(w ~ 1/n^exp(-(b1*DIFF)), c=T)
      ),
    priors=c("b.w~dnorm(0,1)"), n.iter = 1000, n.burnin = 100, seed=1, monitor = T,
    data = boston_df2
  )

mod.bml3 %>% summary()

monetPlot(mod.bml3, parameter="b.w")
```

We find that the degree of dissimilarity between a location and its
neighbors significantly affects the aggregation weights, with more
similar neighbors exerting greater influence than less similar ones.

In other words, features of neighbors that differ significantly from the
focal location in terms of air quality, crime levels, and other factors
have less impact on home values in the focal neighborhood. This insight
would be impossible to uncover using conventional spatial regression
models!

If we set `monitor=T`, we can take a look at the weight estimates:

``` r


data.frame(mod.bml3$w, sum=rowSums(mod.bml3$w, na.rm = T)) %>% head()

rowSums(mod.bml3$w, na.rm = T) %>% sum()
```

The weights sum to 1 for each focal location but vary across neighbors
of each location.

Ben 2do:

- Introduce predict() function  
- Compare predictive performance

### All your friends or just your best friend?

A widely used peer effect model in economics is the linear-in-means
model, which assumes that all peers in a group exert an equal influence
on an individual. In contrast, sociology often emphasizes the
best-friend model, which assumes that an individual’s primary influence
comes from their closest friend. In this example, I use the MMMM to
empirically test which of these models better fits the data.

``` r


data(schoolnets)

SAS <- function(nodedat, vars=NULL, suffix) {

  # Helper function to select node features as _from or _to

  if(is.null(vars)) { # Return all variables if no variable-list specified
    return(
      nodedat %>%
        dplyr::rename_with(~paste0(., suffix), everything())
    )
  } else {
    return(
      nodedat %>%
        dplyr::select(!!vars) %>%
        dplyr::rename_with(~paste0(., suffix), everything())
    )
  }
}

# Create school network data

schoolnets <-

  nodedat %>%

  # add friend features

  left_join(
    edgedat %>%
      select(youthid_from, youthid_to, rank) %>%
      left_join(
        SAS(nodedat, vars=c("youthid", "parent_inc"), suffix="_to"),
        by=c("youthid_to")
      ) %>%
      group_by(youthid_from) %>%
      mutate(
        bestie = rank==1,
        parent_inc_mean = mean(parent_inc_to, na.rm=T),
      ),
    by=c("youthid"="youthid_from")
  ) %>%

  # count number of friends

  tidyr::drop_na() %>%
  add_count(youthid, name="n_friends") %>%
  filter(n_friends>=1) %>%

  mutate(across(everything(), ~as.numeric(.)))


# Naive linear-in-means model
mod4 <- lm(test_cogn ~ sex + etn + parent_edu + parent_inc + parent_inc_mean, data=schoolnets)

# MMMM comparing mean to max aggregation
mod.bml4 <-
  bml(
    test_cogn ~
      sex + etn + parent_edu + parent_inc +
      mm(
        id(youthid_to, youthid),
        mmc(parent_inc_to),
        mmw(w ~ b1*1/n_friends+(1-b1)*bestie, c=T)
      ),
    n.iter = 1000, n.burnin = 100, seed=1, monitor = T, priors = c("b.w~dunif(0,1)"),
    data = schoolnets
  )

mod.bml4 %>% summary()
monetPlot(mod.bml4, parameter="b.w") # posterior distribution of the weight coefficient

mod.bml4$w %>% head() # inspect estimated weight matrix
```

``` r


# The results suggest that the linear-in-means model is appropriate: the best friend does not appear to have a greater impact than the other nominated friends as b1=`r mod.bml4$reg.table %>% filter(variable=="b.w") %>% pull(coefficients) %>% round(.,2)`. Since the estimated aggregation is close to linear-in-means, the effect size of `parent_inc_to` does not differ between models:
```

``` r


coefs <-
  mod.bml4$reg.table[,c(1,2,4,5)] %>%
  filter(!variable %in% c("X0", "sigma.l1", "sigma.l2", "DIC")) %>%
  mutate(model="MMMM") %>%
  add_row(
    data.frame(model="SLM", coefficients=coef(mod4)[-1], confint(mod4, level=0.95)[-1,]) %>%
      rename(lb=X2.5.., ub=X97.5..) %>%
      tibble::rownames_to_column("variable") %>%
      mutate(variable=ifelse(variable=="parent_inc_mean", "parent_inc_to", variable))
  )
rownames(coefs) <- NULL

ggplot(coefs, aes(x=variable, y=coefficients, color=model)) +
  geom_hline(yintercept = 0, color="red") +
  geom_point(position=position_dodge(width=0.3)) +
  geom_pointrange(aes(ymin = lb, ymax = ub), position=position_dodge(width=0.3))+
  labs(title = "Coefficients", x = "Variables", y="", color="Model") +
  theme(legend.position = "bottom") +
  coord_flip()
```

### The effect of political parties’ financial dependency on the survival coalition governments

Here, I replicate the example from the paper, examining whether the
survival of coalition governments is influenced by parties’ financial
dependencies and whether this influence varies across coalition
partners. [Rosche (2025)](https://osf.io/preprints/socarxiv/4bafr_v1).

Tbd
