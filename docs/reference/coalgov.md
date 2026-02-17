# Coalition Governments in Western Democracies (1944-2014)

A dataset containing information on coalition governments and their
member parties across 30 parliamentary democracies. The data are in long
format where the unit of analysis is parties in governments, making it
suitable for multiple-membership multilevel models where governments
(groups) are composed of multiple parties (members).

## Usage

``` r
coalgov
```

## Format

A tibble with 2,077 rows and 18 variables. Each row represents a party's
participation in a specific coalition government. The sample contains
628 governments formed by 312 unique parties across 29 countries.

**Identifiers:**

- gid:

  Government identifier (group-level unit in
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md)
  specification). Range: \[3, 1105\]

- pid:

  Party identifier (member-level unit in
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md)
  specification). Range: \[11110, 96955\]

- cid:

  Country identifier (nesting-level unit in
  [`hm()`](https://benrosche.github.io/bml/reference/hm.md)
  specification). Range: \[11, 96\]

- cname:

  Three-letter country code (ISO 3166-1 alpha-3)

- pname:

  Full party name

**Government-level variables:**

- election:

  Date of the preceding election that led to the government's formation.
  Range: \[1939-04-02, 2014-12-14\]

- n:

  Number of parties in the coalition (group size for weight functions).
  Range: \[2, 9\], mean: 3.31

- dur_wkb:

  Government duration in days, measured from investiture to termination
  (outcome variable for survival models). Range: \[7, 1840\], mean:
  554.5

- event_wkb:

  Early termination indicator: 1 = government terminated due to
  political conflict (voluntary resignation, dissension within
  government, lack of parliamentary support, or head of state
  intervention) more than one year before the official end of term; 0 =
  censored (regular elections, other reasons, or termination within one
  year of scheduled elections). Range: \[0, 1\], mean: 0.39

- majority:

  Majority government indicator: 1 = coalition controls majority of
  parliamentary seats, 0 = minority government. Range: \[0, 1\], mean:
  0.80

- mwc:

  Minimal winning coalition indicator: 1 = coalition would lose its
  majority if any party left, 0 = oversized coalition. Range: \[0, 1\],
  mean: 0.35

- rile_SD:

  Inter-party ideological heterogeneity. Standard deviation of coalition
  parties' left-right positions (from CMP) relative to the ideological
  distribution of all parties in parliament. Standardized and inverted
  so higher values indicate greater ideological cohesion. Range:
  \[-8.40, 2.12\], mean: 0.04

**Country-level variables:**

- investiture:

  Investiture vote requirement (time-constant country characteristic): 1
  = country requires formal parliamentary investiture vote, 0 = no
  formal requirement. Range: \[0, 1\], mean: 0.46

**Party-level variables:**

- pseat:

  Party's relative seat share within the coalition, computed as
  `pseat / sum(pseat)` within each government. Sums to 1 within each
  coalition. Range: \[0.00, 1.00\], mean: 0.33

- prime:

  Prime minister party indicator: `TRUE` = party holds prime
  ministership (n = 628), `FALSE` = junior coalition partner (n = 1,449)

- cohesion:

  Intra-party ideological cohesion, measured using an adaptation of the
  Cowles-Jones ratio. Computed as the ratio of continuous ideological
  shifts to reversals in a party's left-right position over time. Higher
  values indicate more consistent ideological trajectories (greater
  cohesion). Standardized. Range: \[-1.13, 3.85\], mean: 0.00

- rile:

  Party's left-right ideological position (from CMP). Measured on a
  continuous scale where higher values indicate more right-wing
  positions and lower values indicate more left-wing positions.
  Standardized. Range: \[-3.21, 3.68\], mean: 0.00

- finance:

  Party's economic dependence on member contributions (from PPDB).
  Measured as the share of party funding from member dues relative to
  total income. Standardized; higher values indicate greater dependence
  on member financing. Treated as time-constant due to data limitations.
  Range: \[-0.98, 4.40\], mean: 0.00

- Nmembers:

  Number of party members (from PPDB). Standardized; treated as
  time-constant due to data limitations. Range: \[-0.33, 15.02\], mean:
  0.00

## Source

Data compiled from multiple sources:

- **Coalition governments:** Woldendorp, Keman, and Budge (WKB) dataset,
  updated by Seki and Williams (2014)

- **Party ideology:** Comparative Manifesto Project (CMP; Volkens et al.
  2016)

- **Party organization:** Political Party Database (PPDB; Scarrow,
  Poguntke, and Webb 2017)

Missing party-level data imputed using multiple imputation by chained
equations with predictive mean matching.

## Details

This dataset demonstrates multiple-membership multilevel modeling where:

- **Members:** Political parties (identified by `pid`)

- **Groups:** Coalition governments (identified by `gid`)

- **Nesting:** Governments nested within countries (identified by `cid`)

Each coalition government comprises multiple parties, and parties can
participate in multiple governments over time. This creates a
multiple-membership structure where party-level characteristics are
aggregated to the government level using weighting functions specified
in [`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks.

**Sample:** After matching party data across sources and excluding
single-party and caretaker governments, the sample comprises 628
governments formed by 312 unique parties across 29 countries: Australia,
Austria, Belgium, Bulgaria, Croatia, Czech Republic, Denmark, Estonia,
Finland, France, Germany, Greece, Hungary, Ireland, Israel, Italy,
Japan, Latvia, Lithuania, Netherlands, Norway, Poland, Portugal,
Romania, Slovakia, Spain, Sweden, Switzerland, and United Kingdom.

**Measurement notes:**

- Government duration follows the WKB convention: time from investiture
  to termination or new elections

- Early termination events focus on political gridlock (conflict-related
  endings) and exclude terminations within one year of scheduled
  elections

- Party-level variables (`cohesion`, `finance`, `Nmembers`) are
  standardized (mean = 0) for analysis

## References

Seki, K., & Williams, L. K. (2014). Updating the Party Government data
set. *Electoral Studies*, 34, 270-279.

Volkens, A., et al. (2016). The Manifesto Data Collection. Manifesto
Project (MRG/CMP/MARPOR). Version 2016a. Berlin: Wissenschaftszentrum
Berlin fur Sozialforschung.

Scarrow, S. E., Webb, P. D., & Poguntke, T. (Eds.). (2017). *Organizing
Political Parties: Representation, Participation, and Power*. Oxford
University Press.

## See also

[`bml`](https://benrosche.github.io/bml/reference/bml.md) for modeling
examples using this dataset

## Examples

``` r
data(coalgov)

# Explore data structure
str(coalgov)
#> 'data.frame':    1288 obs. of  27 variables:
#>  $ pid        : num  2 6 2 6 2 6 3 5 6 3 ...
#>   ..- attr(*, "label")= chr "Unique party ID"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ pname      : chr  "Social Democratic Labour Party" "Agrarian Party" "Social Democratic Labour Party" "Agrarian Party" ...
#>   ..- attr(*, "label")= chr "Party name"
#>   ..- attr(*, "format.stata")= chr "%-50s"
#>  $ gid        : num  1 1 2 2 3 3 4 4 4 5 ...
#>   ..- attr(*, "label")= chr "Unique government ID"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ gname      : chr  "Erlander III" "Erlander III" "Erlander IV" "Erlander IV" ...
#>   ..- attr(*, "label")= chr "Government name"
#>   ..- attr(*, "format.stata")= chr "%21s"
#>  $ cid        : num  1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "label")= chr "Unique country ID"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ cname      : chr  "Sweden" "Sweden" "Sweden" "Sweden" ...
#>   ..- attr(*, "label")= chr "Country name"
#>   ..- attr(*, "format.stata")= chr "%26s"
#>  $ gstart     : Date, format: "1951-09-30" "1951-09-30" ...
#>  $ gend       : Date, format: "1952-09-21" "1952-09-21" ...
#>  $ n          : num  2 2 2 2 2 2 3 3 3 3 ...
#>   ..- attr(*, "label")= chr "# government parties"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ prime      : num  1 0 1 0 1 0 0 0 1 0 ...
#>   ..- attr(*, "label")= chr "Prime minister party"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ pfam       : hvn_lbll [1:1288] 3, 8, 3, 8, 3, 8, 4, 6, 8, 4, 6, 8, 4, 8, 4, 5, 6, 8...
#>    ..@ label       : chr "Party family"
#>    ..@ format.stata: chr "%13.0g"
#>    ..@ labels      : Named num  1 2 3 4 5 6 7 8 9 10
#>    .. ..- attr(*, "names")= chr [1:10] "ecologist" "communist" "socdem" "liberal" ...
#>  $ rile       : num  -33.4 -4.9 -28.3 1.2 -44.2 1.8 -2.1 2.2 -18.2 -15.2 ...
#>   ..- attr(*, "label")= chr "Right-left position"
#>   ..- attr(*, "format.stata")= chr "%5.2f"
#>  $ ipd        : num  0.25 0 0.25 0 0.25 0 0 0 0 0 ...
#>   ..- attr(*, "label")= chr "Intra-party democracy"
#>   ..- attr(*, "format.stata")= chr "%10.0g"
#>  $ fdep       : num  20.2 1.57 20.2 1.57 20.2 ...
#>   ..- attr(*, "label")= chr "Financial dependency"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ pseatrel   : num  0.577 -0.577 0.618 -0.618 0.696 ...
#>   ..- attr(*, "label")= chr "Party's relative seat share within coalition"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ majority   : num  0 0 0 0 0 0 1 1 1 1 ...
#>   ..- attr(*, "label")= chr "Majority government"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ mwc        : num  1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "label")= chr "Minimal winning coalition"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ hetero     : num  0.49 0.49 0.631 0.631 0.793 ...
#>   ..- attr(*, "label")= chr "SD(rile) of goverment / SD(rile) of parliament"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ investiture: num  1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "label")= chr "Investiture vote"
#>   ..- attr(*, "format.stata")= chr "%10.0g"
#>  $ pmpower    : num  3 3 3 3 3 3 3 3 3 3 ...
#>   ..- attr(*, "label")= chr "Prime ministerial powers"
#>   ..- attr(*, "format.stata")= chr "%10.0g"
#>  $ earlyterm  : num  0 0 0 0 1 1 1 1 1 1 ...
#>   ..- attr(*, "label")= chr "Discretionary early termination"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ govdur     : num  357 357 1466 1466 399 ...
#>   ..- attr(*, "label")= chr "Government duration"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ govmaxdur  : num  357 357 1466 1466 1451 ...
#>   ..- attr(*, "label")= chr "Maximum possible government duration"
#>  $ sim.w      : num  0.326 0.318 0.327 0.319 0.327 ...
#>   ..- attr(*, "label")= chr "Simulated weights"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ sim.y      : num  24.7 24.7 24.4 24.4 23.5 ...
#>   ..- attr(*, "label")= chr "Simulated linear outcome"
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ sim.st     : num  2.42e-10 2.77e-09 5.43e-10 1.92e-10 1.21e-09 ...
#>   ..- attr(*, "label")= chr "Simulated survival time"
#>  $ sim.e      : num  1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "label")= chr "Simulated event status"
table(coalgov$cname)
#> 
#>      Australia        Austria        Belgium Czech Republic        Denmark 
#>             46             45            150             31             76 
#>         France        Germany        Hungary        Ireland         Israel 
#>            182             57             19             34            252 
#>          Italy    Netherlands         Norway         Poland       Portugal 
#>            196             74             40             42             10 
#>          Spain         Sweden United Kingdom 
#>              4             28              2 

# Number of unique units
length(unique(coalgov$gid))   # Governments
#> [1] 402
length(unique(coalgov$pid))   # Parties
#> [1] 194
length(unique(coalgov$cid))   # Countries
#> [1] 18

if (FALSE) { # \dontrun{
# Model: government duration as function of majority status and party characteristics
m1 <- bml(
  Surv(dur_wkb, event_wkb) ~ 1 + majority +
    mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = TRUE) +
    hm(id = id(cid), type = "RE"),
  family = "Weibull",
  data = coalgov
)
summary(m1)
} # }
```
