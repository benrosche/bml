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
        id = id(tid_nb, tid),
        vars = NULL,
        fn = fn(w ~ 1/n, c = TRUE),
        RE = TRUE
      ),
    n.iter = 1000, n.burnin = 100, seed = 1, monitor = TRUE,
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
        id = id(tid_nb, tid),
        vars = vars(NOX_nb + CRIM_nb + RM_nb + DIS_nb + AGE_nb),
        fn = fn(w ~ 1/n, c = TRUE),
        RE = TRUE
      ),
    n.iter = 1000, n.burnin = 100, seed = 1, monitor = TRUE,
    data = boston_df
  )

mod.bml2 %>% summary()

# Calculate spatial correlation in the residual
getLambda <- function(x) {
  s.mm <- x$reg.table["sigma.mm", "coefficients"]
  s    <- x$reg.table["sigma", "coefficients"]
  return(s.mm^2/(s.mm^2+s^2))
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
`lnCMEDV ~ NOX + CRIM + RM + DIS + AGE + mm(id = id(tid_nb, tid), vars = vars(NOX_nb + CRIM_nb + RM_nb + DIS_nb + AGE_nb), fn = fn(w ~ 1/n, c = TRUE), RE = TRUE)`

The only difference compared to a
[`lm()`](https://rdrr.io/r/stats/lm.html) formula is the inclusion of
the [`mm()`](https://benrosche.github.io/bml/reference/mm.md) container.
Within this container, we define three sub-containers:  
- `ids()` for the identifiers  
- `vars = NULL` for the covariates being considered,  
- `fn = fn()` to endogenize the weight function. Here, `w ~ 1/n`
specifies the row-standardized weight.

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
  filter(!variable %in% c("X0", "sigma.mm", "sigma", "DIC"))

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

The weight function is specified within the `fn = fn()` container,
allowing for any function that produces bounded weights. Logical
operators such as *==, \>, \<* can be used to define conditions that
enable aggregation functions, such as `min` or `max`.

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
constraints using the `c` parameter in
[`fn()`](https://benrosche.github.io/bml/reference/fn.md):

- **`c = TRUE`**: Ensures that the weights of a location’s neighbors sum
  to 1 for each focal location, while allowing them to vary within each
  location.  
- **`c = FALSE`**: No constraint applied; weights can vary freely both
  within and across locations.

Both options identify the model but have different substantive
interpretations. Here, we will use weights that sum to 1 for each focal
location while allowing variation within locations (`c = TRUE`):

``` r

# Spatial fixed effects + spatial random effect with "endogenized" weights:
mod.bml3 <-
  bml(
    lnCMEDV ~
      NOX + CRIM + RM + DIS + AGE +
      mm(
        id = id(tid_nb, tid),
        vars = vars(NOX_nb + CRIM_nb + RM_nb + DIS_nb + AGE_nb),
        fn = fn(w ~ 1/n^exp(-(b1*DIFF)), c = TRUE),
        RE = TRUE
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

### Peer Effects in High School Networks: Equal Influence or Reciprocity-Weighted?

A widely used peer effect model in economics is the **linear-in-means
model**, which assumes that all peers in a group exert equal influence
on an individual. An alternative specification allows peer influence to
vary based on relationship characteristics, such as reciprocity or
strength of ties. In this example, I use the MMMM to empirically test
which aggregation mechanism better fits friendship network data.

We use the `faux.highschool` dataset from the `ergm.count` package,
which contains friendship ties among high school students. For this
demonstration, we’ll create a simulated outcome variable representing
academic engagement and test whether peer characteristics have equal
effects or whether reciprocated friendships carry more weight.

``` r

library(bml)
library(dplyr)
library(tidyr)
library(ggplot2)
library(statnet)

data(faux.mesa.high)

# Extract node attributes
nodedat <- data.frame(
  student_id = network.vertex.names(faux.highschool),
  grade = faux.highschool %v% "grade",
  sex = faux.highschool %v% "sex",
  race = faux.highschool %v% "race"
) %>%
  mutate(student_id = as.numeric(student_id))

# Extract edge list (friendship ties)
edgedat <- as.data.frame(as.edgelist(faux.highschool)) %>%
  rename(from = V1, to = V2) %>%
  mutate(
    from = as.numeric(as.character(from)),
    to = as.numeric(as.character(to))
  )

# Identify reciprocated ties (mutual friendships)
edgedat <- edgedat %>%
  mutate(
    dyad = paste(pmin(from, to), pmax(from, to), sep = "-")
  ) %>%
  group_by(dyad) %>%
  mutate(
    reciprocated = n() == 2 # TRUE if both i->j and j->i exist
  ) %>%
  ungroup() %>%
  select(-dyad)

# Simulate an academic engagement outcome for demonstration
# In practice, this would be measured empirically
set.seed(123)
nodedat <- nodedat %>%
  mutate(
    ses = rnorm(n()), # Simulated socioeconomic status
    engagement = 50 + 5 * (grade - 10) + 3 * ses + rnorm(n(), sd = 10)
  )

# Create analysis dataset by merging node and edge data
network_data <- edgedat %>%
  left_join(
    nodedat %>% select(student_id, ses, engagement),
    by = c("from" = "student_id")
  ) %>%
  rename(ses_from = ses, engagement_from = engagement) %>%
  left_join(
    nodedat %>% select(student_id, grade, sex, race, ses, engagement),
    by = c("to" = "student_id")
  ) %>%
  rename(
    ses_to = ses,
    engagement_to = engagement,
    grade_to = grade,
    sex_to = sex,
    race_to = race
  ) %>%
  group_by(from) %>%
  mutate(
    n_friends = n(),
    ses_friends_mean = mean(ses_to, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  drop_na()


# Model 1: Naive linear-in-means model (OLS)
mod_lm <- lm(
  engagement_from ~ grade_to + sex_to + race_to + ses_from + ses_friends_mean,
  data = network_data
)
summary(mod_lm)

# Model 2: MMMM with equal weights (linear-in-means via MM framework)
mod_equal <- bml(
  engagement_from ~ 1 +
    mm(
      id = id(to, from),
      vars = vars(ses_to),
      fn = fn(w ~ 1 / n, c = TRUE),
      RE = FALSE
    ),
  family = "Gaussian",
  n.iter = 2000,
  n.burnin = 500,
  seed = 1,
  data = network_data
)
summary(mod_equal)

# Model 3: MMMM with reciprocity-weighted aggregation
# Test whether reciprocated friendships carry more weight
mod_reciprocity <- bml(
  engagement_from ~ 1 +
    mm(
      id = id(to, from),
      vars = vars(ses_to),
      fn = fn(w ~ 1 / n^exp(b1 * reciprocated), c = TRUE),
      RE = FALSE
    ),
  family = "Gaussian",
  priors = c("b.w ~ dnorm(0, 1)"),
  n.iter = 2000,
  n.burnin = 500,
  seed = 1,
  monitor = TRUE,
  data = network_data
)
summary(mod_reciprocity)

# Visualize the weight function coefficient
monetPlot(mod_reciprocity, parameter = "b.w")
```

**Interpreting the results:**

- If the weight coefficient (`b.w`) is close to 0, reciprocity doesn’t
  affect peer influence—supporting the linear-in-means model.
- If `b.w` is significantly positive, reciprocated friendships carry
  more weight, suggesting that mutual friends have stronger peer
  effects.
- The MMMM framework allows us to empirically test these competing
  mechanisms rather than imposing one a priori.

``` r

# Compare model fit using DIC
dic_equal <- mod_equal$DIC
dic_reciprocity <- mod_reciprocity$DIC

cat("DIC (equal weights):", dic_equal, "\n")
cat("DIC (reciprocity-weighted):", dic_reciprocity, "\n")
cat("Difference:", dic_equal - dic_reciprocity, "\n")

# Visualize coefficient comparison
coefs <- bind_rows(
  mod_equal$reg.table %>%
    filter(variable == "ses_to") %>%
    mutate(model = "Equal weights"),
  mod_reciprocity$reg.table %>%
    filter(variable == "ses_to") %>%
    mutate(model = "Reciprocity-weighted")
)

ggplot(coefs, aes(x = model, y = coefficients)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2) +
  labs(
    title = "Peer SES Effect on Academic Engagement",
    subtitle = "Comparing equal vs. reciprocity-weighted aggregation",
    x = "Model",
    y = "Coefficient estimate (95% CI)"
  ) +
  theme_minimal()
```

This example demonstrates how the MMMM allows researchers to test
substantive theories about aggregation mechanisms in network data,
moving beyond the assumption of equal peer influence.

\`\`\`

### The effect of political parties’ financial dependency on the survival coalition governments

Here, I replicate the example from the paper, examining whether the
survival of coalition governments is influenced by parties’ financial
dependencies and whether this influence varies across coalition
partners. [Rosche (2025)](https://osf.io/preprints/socarxiv/4bafr_v1).

Tbd
