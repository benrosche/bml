# Model Structure

## Overview

The `bml` package implements Bayesian **multiple-membership multilevel
models (MMMM)** with **parameterizable weight functions**. These models
are designed for data structures where:

1.  **Lower-level units (members)** can belong to **multiple
    higher-level units (groups)**
2.  **Groups** are composed of **multiple members**
3.  Group outcomes depend on the combined influence of their constituent
    members

This structure is common across many empirical settings: political
parties participating in multiple coalition governments, individuals
participating in multiple teams, spatial units influenced by multiple
neighbors, or actors embedded in multiple networks.

The key innovation of the extended MMMM is that it allows researchers to
**model rather than impose** the aggregation mechanism through which
member-level effects combine to influence group-level outcomes.

## Notation and Data Structure

Let:

- $`i = 1, \ldots, N`$ index **groups** (the primary unit of analysis)
- $`k = 1, \ldots, K`$ index **members** (lower-level units)
- $`j = 1, \ldots, J`$ index **nesting units** (higher-level contexts)
- $`p[i]`$ denote the set of members belonging to group $`i`$
- $`n_i = |p[i]|`$ denote the number of members in group $`i`$

**Example (Coalition Governments):** Groups are governments ($`i`$),
members are political parties ($`k`$), and nesting units are countries
($`j`$). The function $`p[i]`$ returns all parties in government $`i`$,
and $`n_i`$ is coalition size.

The outcome $`y_i^{(g)}`$ is measured at the group level, where
superscript $`(g)`$ denotes the level of measurement. In what follows,
we present the model for Gaussian (normal) outcomes, though the
framework extends to other distributions (see Extensions section).

## The Extended Multiple-Membership Multilevel Model

### Full Model Specification

The extended MMMM decomposes a group-level outcome into three components
corresponding to the levels at which explanatory factors operate:

``` math
\begin{equation}
y_i^{(g)} = \theta_i^{(g)} + \theta_{j[i]}^{(c)} + \theta_i^{(p)} \tag{1}
\end{equation}
```

where:

- $`\theta_i^{(g)}`$ is the **group-level effect**
- $`\theta_{j[i]}^{(c)}`$ is the **nesting-level effect** (for nesting
  unit $`j`$ containing group $`i`$)
- $`\theta_i^{(p)}`$ is the **aggregated member-level effect**

Each component includes both systematic effects (observed covariates)
and random effects (unobserved heterogeneity).

### Group-Level Effect

The group-level effect models the direct influence of group-level
covariates and group-specific unobserved factors:

``` math
\begin{equation}
\theta_i^{(g)} = \boldsymbol{\beta}^{(g)} \cdot \mathbf{x}_i^{(g)} + u_i^{(g)}, \quad u_i^{(g)} \overset{\text{iid}}{\sim} N(0, \sigma_{u^{(g)}}^2) \tag{2}
\end{equation}
```

where:

- $`\mathbf{x}_i^{(g)} = [1, x_{1i}^{(g)}, \ldots, x_{Gi}^{(g)}]`$ is a
  vector of group-level covariates
- $`\boldsymbol{\beta}^{(g)} = [\beta_0^{(g)}, \beta_1^{(g)}, \ldots, \beta_G^{(g)}]`$
  are the corresponding regression coefficients
- $`u_i^{(g)}`$ are independently and identically distributed (i.i.d.)
  random effects
- $`\sigma_{u^{(g)}}^2`$ captures unobserved heterogeneity at the group
  level

**Example:** For government survival, $`\mathbf{x}_i^{(g)}`$ might
include majority status, ideological heterogeneity, or economic
conditions.

### Nesting-Level Effect

When groups are nested within higher-level contexts, the nesting-level
effect accounts for hierarchical dependencies:

``` math
\begin{equation}
\theta_j^{(c)} = \boldsymbol{\beta}^{(c)} \cdot \mathbf{x}_j^{(c)} + u_j^{(c)}, \quad u_j^{(c)} \overset{\text{iid}}{\sim} N(0, \sigma_{u^{(c)}}^2) \tag{3}
\end{equation}
```

where:

- $`\mathbf{x}_j^{(c)} = [x_{1j}^{(c)}, \ldots, x_{Cj}^{(c)}]`$ are
  nesting-level covariates
- $`\boldsymbol{\beta}^{(c)} = [\beta_1^{(c)}, \ldots, \beta_C^{(c)}]`$
  are the corresponding coefficients
- $`u_j^{(c)}`$ are i.i.d. random effects
- $`\sigma_{u^{(c)}}^2`$ captures unobserved heterogeneity at the
  nesting level

Including $`u_j^{(c)}`$ allows groups within the same nesting unit to be
more similar than groups across nesting units, even after conditioning
on observed covariates.

**Example:** Governments are nested within countries. Country-level
variables might include electoral system or investiture requirements.
The random effect captures country-specific factors affecting all
governments in that country.

**Note:** Multiple
[`hm()`](https://benrosche.github.io/bml/reference/hm.md) blocks can be
specified to accommodate cross-classified structures (e.g., groups
nested within multiple non-hierarchical contexts).

### Aggregated Member-Level Effect

The core innovation of the extended MMMM is modeling how member-level
effects aggregate to influence group outcomes. Individual member effects
are first specified, then aggregated via a weighted sum:

``` math
\begin{equation}
\theta_i^{(p)} = \sum_{k \in p[i]} w_{ik} \left( \boldsymbol{\beta}^{(p)} \cdot \mathbf{x}_{ik}^{(p)} + u_k^{(p)} \right) \tag{4}
\end{equation}
```

This can be decomposed into systematic and random components:

``` math
\begin{equation}
\theta_i^{(p)} = \underbrace{\boldsymbol{\beta}^{(p)} \sum_{k \in p[i]} w_{ik} \mathbf{x}_{ik}^{(p)}}_{\text{Systematic component}} + \underbrace{\sum_{k \in p[i]} w_{ik} u_k^{(p)}}_{\text{Random component}}, \quad u_k^{(p)} \overset{\text{iid}}{\sim} N(0, \sigma_{u^{(p)}}^2) \tag{5}
\end{equation}
```

where:

- $`\mathbf{x}_{ik}^{(p)} = [x_{1ik}^{(p)}, \ldots, x_{Pik}^{(p)}]`$ are
  member-level covariates
- $`\boldsymbol{\beta}^{(p)} = [\beta_1^{(p)}, \ldots, \beta_P^{(p)}]`$
  are structural regression coefficients measuring the effect of member
  attributes on group outcomes
- $`w_{ik}`$ are **aggregation weights** determining how member $`k`$’s
  contribution enters group $`i`$’s outcome
- $`u_k^{(p)}`$ are i.i.d. member-specific random effects
- $`\sigma_{u^{(p)}}^2`$ captures unobserved heterogeneity at the member
  level

**Interpretation of Weights:** The weights $`w_{ik}`$ define the
**micro-macro link**—how individual member characteristics combine to
produce group-level outcomes. Since the structural effect
$`\boldsymbol{\beta}^{(p)}`$ is constant, the weights serve to aggregate
the member-level covariate $`\mathbf{x}_{ik}^{(p)}`$ into a group-level
construct $`\sum_{k \in p[i]} w_{ik} \mathbf{x}_{ik}^{(p)}`$.

**Covariance Structure:** Including member-specific random effects
$`u_k^{(p)}`$ induces covariance among groups sharing members. The
covariance between groups $`i`$ and $`i'`$ is:

``` math
\begin{equation}
\text{Cov}(y_i^{(g)}, y_{i'}^{(g)}) = \sigma_{u^{(p)}}^2 \sum_k w_{ik} w_{i'k}
\end{equation}
```

This increases with the degree of member overlap and accounts for
dependencies across groups involving recurring members.

**Example:** In coalition governments, $`\mathbf{x}_{ik}^{(p)}`$ might
include party cohesion or organizational structure. The random effect
$`u_k^{(p)}`$ captures party-specific unobserved factors that affect all
governments in which party $`k`$ participates.

## Parameterizable Weight Functions

### Conventional MMMM: Fixed Weights

The **conventional MMMM** requires weights to be specified *a priori*.
The most common choice is equal weighting (arithmetic mean):

``` math
\begin{equation}
w_{ik} = \frac{1}{n_i} \tag{6}
\end{equation}
```

This assumes all members contribute equally to group outcomes,
regardless of their attributes or context.

**Limitations:** Equal weighting is often theoretically implausible. In
coalition governments, for instance, larger parties, parties holding the
prime ministership, or ideologically central parties may exert
disproportionate influence. In spatial analysis, nearby neighbors may
matter more than distant ones. Imposing equal weights when they vary in
practice introduces **aggregation bias**.

### Extended MMMM: Endogenizing Weights

The **extended MMMM** allows weights to be modeled as functions of
covariates and estimated parameters:

``` math
\begin{equation}
w_{ik}(\mathbf{x}_{ik}^{(w)}, n_i) = \frac{1}{1 + (n_i - 1) \exp\left(-\boldsymbol{\beta}^{(w)} \cdot \mathbf{x}_{ik}^{(w)}\right)} \tag{7}
\end{equation}
```

subject to the normalization constraint $`\sum_{k \in p[i]} w_{ik} = 1`$
for all $`i`$.

where:

- $`\mathbf{x}_{ik}^{(w)} = [1, x_{1ik}^{(w)}, \ldots, x_{Wik}^{(w)}]`$
  are **weight covariates**
- $`\boldsymbol{\beta}^{(w)} = [\beta_0^{(w)}, \beta_1^{(w)}, \ldots, \beta_W^{(w)}]`$
  are **weight function parameters** to be estimated

**Functional Form:** Equation (7) generalizes the logistic function by
centering it at $`1/n_i`$ when
$`\boldsymbol{\beta}^{(w)} \cdot \mathbf{x}_{ik}^{(w)} = 0`$. This sets
equal weighting as the baseline: absent systematic differences, each
member receives weight $`1/n_i`$. When
$`\boldsymbol{\beta}^{(w)} \neq 0`$, weights deviate from this baseline,
reflecting heterogeneous member influence.

**Properties:**

- Bounded: $`w_{ik} \in [0, 1]`$
- Baseline-centered: When $`\boldsymbol{\beta}^{(w)} = 0`$, weights
  reduce to $`w_{ik} = 1/n_i`$
- Monotonic: Weights increase/decrease smoothly with the linear
  predictor
- Normalized: Constraint ensures $`\sum_k w_{ik} = 1`$ (unless
  normalization is disabled)

**Interpretation of Weight Parameters:** While
$`\boldsymbol{\beta}^{(p)}`$ measure *structural effects* of member
attributes on group outcomes, $`\boldsymbol{\beta}^{(w)}`$ measure
*aggregation effects*—how member attributes determine their relative
importance in the weighted sum. Together, they reveal both *what*
matters (structural effects) and *how much* different members’
characteristics matter (aggregation effects).

**Example:** Let $`x_{ik}^{(p)}`$ = party cohesion and $`x_{ik}^{(w)}`$
= indicator for holding the prime ministership in a government survival
model.

- If $`\beta^{(w)} = 0`$: Prime minister’s cohesion matters as much as
  coalition partners’ cohesion (equal weights)
- If $`\beta^{(w)} > 0`$: Prime minister’s cohesion carries more weight
- If $`\beta^{(w)} < 0`$: Coalition partners’ cohesion carries more
  weight

### Testing Alternative Aggregation Mechanisms

Weight function regression enables empirical tests of substantive
hypotheses about aggregation:

**1. Mean vs. Extremes:** Does the group outcome depend on the average
member attribute or on extreme values (minimum/maximum)?

Create indicator variables for $`\text{argmin}_k(x_{ik}^{(p)})`$ and
$`\text{argmax}_k(x_{ik}^{(p)})`$ and include them in
$`\mathbf{x}_{ik}^{(w)}`$. If aggregation is better captured by
extremes, the corresponding member’s weight will approach 1.

**2. Mean vs. Sum:** Are member effects *contextual* (mean aggregation:
$`\sum_k w_{ik} = 1`$) or *autonomous* (sum aggregation:
$`\sum_k w_{ik} = n_i`$)?

Disable normalization (`c = FALSE` in
[`fn()`](https://benrosche.github.io/bml/reference/fn.md)). The
intercept $`\beta_0^{(w)}`$ then reveals the aggregation form: -
$`\beta_0^{(w)} \approx 0`$: Mean aggregation - $`\beta_0^{(w)} \gg 0`$:
Sum aggregation

Under mean aggregation, a one-unit change in $`x_{ik}^{(p)}`$ affects
the outcome by $`\beta^{(p)}`$ only if *all* members change. Under sum
aggregation, a one-unit change by a *single* member already produces an
effect of $`\beta^{(p)}`$.

**3. Asymmetric Influence:** Do certain members (e.g., larger, more
central, more powerful) carry disproportionate weight?

Include member attributes hypothesized to determine influence (e.g.,
size, centrality, power) in $`\mathbf{x}_{ik}^{(w)}`$.

## Variance Decomposition and Intraclass Correlation

The multilevel structure partitions total variance into group-level,
nesting-level, and member-level components:

``` math
\begin{equation}
\text{Var}(y_i^{(g)}) = \sigma_{u^{(g)}}^2 + \sigma_{u^{(c)}}^2 + \sigma_{u^{(p)}}^2 \sum_{k \in p[i]} w_{ik}^2
\end{equation}
```

**Intraclass Correlation Coefficients (ICC):**

The proportion of variance attributable to each level quantifies the
importance of group-level, nesting-level, and member-level unobserved
factors:

``` math
\begin{align}
\rho^{(g)} &= \frac{\sigma_{u^{(g)}}^2}{\sigma_{u^{(g)}}^2 + \sigma_{u^{(c)}}^2 + \sigma_{u^{(p)}}^2 \mathbb{E}[\sum_k w_{ik}^2]} \\
\rho^{(c)} &= \frac{\sigma_{u^{(c)}}^2}{\sigma_{u^{(g)}}^2 + \sigma_{u^{(c)}}^2 + \sigma_{u^{(p)}}^2 \mathbb{E}[\sum_k w_{ik}^2]} \\
\rho^{(p)} &= \frac{\sigma_{u^{(p)}}^2 \mathbb{E}[\sum_k w_{ik}^2]}{\sigma_{u^{(g)}}^2 + \sigma_{u^{(c)}}^2 + \sigma_{u^{(p)}}^2 \mathbb{E}[\sum_k w_{ik}^2]}
\end{align}
```

When weights are equal ($`w_{ik} = 1/n_i`$), the member-level variance
component simplifies to $`\sigma_{u^{(p)}}^2 / \bar{n}`$ where
$`\bar{n}`$ is average group size.

## Extensions

### Multiple MM Blocks

The `bml` package allows multiple
[`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks with
distinct aggregation mechanisms. This is useful when different
member-level covariates aggregate through different processes:

``` math
\begin{equation}
\theta_i^{(p)} = \sum_{b=1}^B \sum_{k \in p[i]} w_{ik}^{(b)} \left( \boldsymbol{\beta}^{(p,b)} \cdot \mathbf{x}_{ik}^{(p,b)} + u_k^{(p,b)} \right)
\end{equation}
```

where $`b`$ indexes MM blocks, each with its own weight function
$`w_{ik}^{(b)} = f(\mathbf{x}_{ik}^{(w,b)}, n_i)`$.

**Example:** Model both observed party characteristics (with estimated
weights) and party random effects (with equal weights):

``` r

mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n^exp(b*primeminister)), RE = FALSE) +
mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE)
```

**Important:** Only one
[`mm()`](https://benrosche.github.io/bml/reference/mm.md) block can have
`RE = TRUE` in a given model.

### Autoregressive Random Effects

Member random effects may evolve over time rather than remain constant.
The `ar = TRUE` option specifies a random walk process:

``` math
\begin{equation}
u_{ik}^{(p)} \sim \begin{cases}
N(0, \sigma_{u^{(p)}}^2) & \text{if } g[i,k] = \emptyset \\
N(u_{g[i,k],k}^{(p)}, \sigma_{u^{(p)}}^2) & \text{if } g[i,k] \neq \emptyset
\end{cases}
\end{equation}
```

where $`g[i,k]`$ returns the most recent previous group in which member
$`k`$ participated. For $`g[i,k] = \emptyset`$ (first participation),
the effect is drawn from the prior.

**Interpretation:** Autoregressive effects capture dynamics where a
member’s unobserved heterogeneity changes across successive group
participations. The variance $`\sigma_{u^{(p)}}^2`$ now represents the
*innovation variance* governing temporal change.

**Example:** A party’s propensity to destabilize governments may evolve
as its leadership, organization, or electoral fortunes change.

**Usage:**

``` r

mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE, ar = TRUE)
```

### Opposition Members

When analyzing group outcomes influenced by both members and non-members
(e.g., coalition and opposition parties), separate
[`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks can
model their distinct contributions:

``` math
\begin{equation}
\theta_i^{(p)} = \sum_{k \in p[i]} w_{ik}^{(g)} \left( \boldsymbol{\beta}^{(g)} \cdot \mathbf{x}_{ik}^{(g)} + u_k^{(g)} \right) + \sum_{k \notin p[i]} w_{ik}^{(o)} \left( \boldsymbol{\beta}^{(o)} \cdot \mathbf{x}_{ik}^{(o)} + u_k^{(o)} \right)
\end{equation}
```

where superscripts $`(g)`$ and $`(o)`$ denote member and non-member
effects, respectively. This allows both groups to have different
structural effects and different aggregation mechanisms.

### Generalized Outcomes

The model structure extends to non-Gaussian outcomes via generalized
linear and survival models:

**Binomial (Logistic Regression):**

``` math
\begin{equation}
\text{logit}\left(P(y_i^{(g)} = 1)\right) = \theta_i^{(g)} + \theta_{j[i]}^{(c)} + \theta_i^{(p)}
\end{equation}
```

**Weibull Survival Model:**

``` math
\begin{equation}
\log(t_i) = \theta_i^{(g)} + \theta_{j[i]}^{(c)} + \theta_i^{(p)} + \sigma \epsilon_i
\end{equation}
```

where $`t_i`$ is survival time, $`\epsilon_i`$ follows an extreme value
distribution, and $`\sigma`$ is the scale parameter related to the
Weibull shape parameter.

**Cox Proportional Hazards Model:**

``` math
\begin{equation}
h(t_i) = h_0(t) \exp\left(\theta_i^{(g)} + \theta_{j[i]}^{(c)} + \theta_i^{(p)}\right)
\end{equation}
```

where $`h_0(t)`$ is the baseline hazard function, estimated
non-parametrically or as a piecewise constant function.

All outcome types support the full multilevel and multiple-membership
structure, including parameterizable weight functions.

## Comparison with Conventional MMMM

The **conventional MMMM** (as implemented in `brms` or `MLwiN`) is a
special case:

``` math
\begin{equation}
y_i^{(g)} = \boldsymbol{\beta}^{(g)} \cdot \mathbf{x}_i^{(g)} + \boldsymbol{\beta}^{(p)} \cdot \bar{\mathbf{x}}_{i}^{(p)} + \sum_{k \in p[i]} w_{ik} u_k^{(p)} + u_i^{(g)} \tag{Conventional MMMM}
\end{equation}
```

where
$`\bar{\mathbf{x}}_{i}^{(p)} = \sum_k w_{ik} \mathbf{x}_{ik}^{(p)}`$ is
pre-aggregated using *imposed* weights $`w_{ik}`$.

The **extended MMMM** generalizes this by:

1.  **Aggregating both systematic and random components** within the
    model
2.  **Endogenizing weights** as functions of covariates with estimable
    parameters
3.  **Enabling multiple aggregation mechanisms** via multiple
    [`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks

``` math
\begin{equation}
y_i^{(g)} = \boldsymbol{\beta}^{(g)} \cdot \mathbf{x}_i^{(g)} + \sum_{k \in p[i]} w_{ik} \left( \boldsymbol{\beta}^{(p)} \cdot \mathbf{x}_{ik}^{(p)} + u_k^{(p)} \right) + u_i^{(g)}, \quad w_{ik} = f(\mathbf{x}_{ik}^{(w)}, n_i) \tag{Extended MMMM}
\end{equation}
```

## Model Assumptions

The extended MMMM assumes:

1.  **Independence of random components** across levels:
    $`u_i^{(g)} \perp u_j^{(c)} \perp u_k^{(p)}`$. If violated, variance
    estimates may be biased. More complex covariance structures can be
    accommodated if needed.

2.  **Random effects are uncorrelated with covariates** (random effects
    assumption). Including within-cluster means of covariates or using
    within-between specifications isolates within-cluster effects,
    yielding identical point estimates to fixed effect models.

3.  **Correct functional form** for weight function. The logistic-type
    function in Equation (7) is flexible, but other forms can be
    specified using the
    [`fn()`](https://benrosche.github.io/bml/reference/fn.md) syntax.

4.  **Proper model specification:** Misspecifying a weight covariate as
    a structural covariate (or vice versa) yields null effects for the
    misspecified variable while correctly identifying properly specified
    ones, provided there is no omitted variable bias.

## Estimation

Model parameters are estimated using **Bayesian Markov chain Monte Carlo
(MCMC)** via JAGS because:

1.  The likelihood has no closed-form solution due to the parameterized
    weight function
2.  Bayesian estimation captures uncertainty in variance estimates,
    yielding credible intervals even when maximum likelihood methods
    report zero variance
3.  Weakly informative priors can be specified for all parameters to
    regularize estimates in small samples

The `bml` package provides a user-friendly interface to JAGS with
automated prior specification, convergence diagnostics, and
post-estimation tools.

## Summary

The extended MMMM provides a flexible framework for analyzing group
outcomes when:

- Groups are composed of multiple members
- Members participate in multiple groups
- The aggregation of member-level effects to group-level outcomes is
  theoretically uncertain

By modeling rather than imposing aggregation weights, researchers can
test substantive hypotheses about micro-macro linkages, avoid
aggregation bias, and account for complex dependence structures in
multilevel data.
