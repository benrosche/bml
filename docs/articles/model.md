# Model structure

Let $`y_{i}^{(2)}`$ be an outcome at level 2, where subscript $`i`$
indexes level 2 units and superscript $`(2)`$ indicates that the outcome
is located at the second level. Using the generalized multiple
membership multilevel model with endogenized weights (MMMM), we can
model this outcome in terms of an aggregated level 1 effect
$`\theta_{i}^{(1)}`$ and a level 2 effect $`\theta_{i}^{(2)}`$:

``` math
\begin{equation} y_{i}^{(2)}=\theta_{i}^{(1)}+\theta_{i}^{(2)} \tag{1} \end{equation}
```

The level 2 effect in equation (1) is determined by a systematic
component of observed explanatory variables at level 2 $`x_{i}^{(2)}`$
with regression coefficients $`\beta^{(2)}`$and a random component
$`u_{i}^{(2)}`$:

``` math
\begin{equation} \theta_{i}^{(2)}=x_{i}^{\intercal (2)}\beta^{(2)}+u_{i}^{(2)} \tag{2} \end{equation}
```

The random component at this level is a disturbance term, which is
assumed to be normally distributed with a mean of zero and constant
variance $`u_{i}^{(2)}\sim N(0,\sigma_{u^{(2)}}^{2})`$. This part of the
model represents the conventional single-level model structure.

The aggregated level 1 effect $`\theta_{i}^{(1)}`$ in equation (1)
models the aggregation of the effects of level 1 units to the second
level. It is determined by a weighted sum of the effect of each level 1
unit $`j`$ on the level 2 outcome in the set of level 1 units $`z(i)`$
embedded in level 2 unit $`i`$:

``` math
\begin{equation} \theta_{i}^{(1)}=\sum_{j \in z(i)}w_{ij}\zeta_{ij} \tag{3} \end{equation}
```

That is, subscript $`j=1,...,J`$ indexes level 1 units and the indexing
function $`z(i)`$ returns all level 1 units that are members of level 2
unit $`i`$. The multiple membership construct aggregates individual
level 1 effects $`\zeta_{ij}`$ by taking their weighted sum with weights
$`w_{ij}=w_{ij}^{*}`$ for all parties $`j \in z(i)`$ and $`w_{ij}=0`$
for all parties $`j \notin z(i)`$.

The individual-level 1 effects $`\zeta_{ij}`$ are determined by a
systematic component of observed explanatory variables at level 1
$`x_{ij}^{(1)}`$ with regression coefficients $`\beta^{(1)}`$ and a
random component $`u_{ij}^{(1)}`$, representing the joint impact of
unobserved variables:

``` math
\begin{equation} \zeta_{ij}=x_{ij}^{\intercal (1)}\beta^{(1)}+u_{ij}^{(1)} \tag{4} \end{equation}
```

The random component at this level is also assumed to be normally
distributed with a mean of zero and constant variance
$`u_{ij}^{(1)}\sim N(0,\sigma_{u^{(1)}}^{2})`$.

To examine how the effects of level 1 units propagates to the second
level, we endogenize the weights instead of assigning fixed values to
each unit:

``` math
\begin{equation} \begin{split} w_{ij}=\frac{1}{n_{i}^{exp(-(x_{ij}^{\intercal W}\beta^{W}))}} \\ s.t. \sum_{i} \sum_{j} w_{ij}=N \end{split} \tag{5} \end{equation}
```

where $`x_{ij}^{W}`$ are explanatory variables with regression
coefficients $`\beta^{W}`$, $`n_{i}`$ is the number of members level 2
unit $`i`$, and $`N`$ equals the total number of observations in the
dataset. In this form, the weights are bounded by $`[0,1]`$.

The weight regression coefficients estimate the impact of explanatory
variables on the specific weight of a level 1 unit in its effect on the
level 2 outcome. If the weight variables have no impact on the
aggregation process, i.e. $`\beta^{W}=0`$, the weights reduce to
$`w_{ij}=\frac{1}{n_{i}}`$ (mean aggregation). If $`\beta^{W} \neq 0`$,
the weights reveal a more complex interplay of level 1 units in their
effect on the level 2 outcome. That is, weights will deviate from
$`\frac{1}{n_{i}}`$ and are no longer constant within and between level
2 units. Instead, they depend on $`x_{ij}^{W}`$.
