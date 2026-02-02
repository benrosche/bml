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

A tibble with 2,077 rows and 21 variables. Each row represents a party's
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

- country:

  Three-letter country code (ISO 3166-1 alpha-3)

- pname:

  Full party name

**Government-level variables:**

- pelection:

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

- comp_early:

  Early election indicator for competing risks analysis: 1 = government
  terminated by calling early elections, 0 = otherwise. Sourced from
  WKB. Range: \[0, 1\], mean: 0.04

- comp_replace:

  Nonelectoral replacement indicator for competing risks analysis: 1 =
  government terminated by nonelectoral replacement (cabinet reshuffle
  without elections), 0 = otherwise. Sourced from WKB. Range: \[0, 1\],
  mean: 0.36

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

  Party's proportional seat share within the coalition. Range: \[0.00,
  73.00\], mean: 0.22

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

- `comp_early` and `comp_replace` enable competing risks analysis
  distinguishing between termination mechanisms

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

# Explore structure
head(coalgov)
#> # A tibble: 6 × 21
#>     gid   pid pname   pseat prime   rile cohesion finance Nmembers   cid country
#>   <int> <dbl> <chr>   <dbl> <lgl>  <dbl>    <dbl>   <dbl>    <dbl> <dbl> <chr>  
#> 1     3 11320 Socia… 0.295  TRUE  -1.58     0.670  -0.929   -0.105    11 SWE    
#> 2     3 11810 Agrar… 0.0789 FALSE -0.228   -0.205  -0.983   -0.264    11 SWE    
#> 3     4 11320 Socia… 0.289  TRUE  -1.58     0.670  -0.929   -0.105    11 SWE    
#> 4     4 11810 Agrar… 0.0684 FALSE -0.228   -0.205  -0.983   -0.264    11 SWE    
#> 5     5 11320 Socia… 0.279  TRUE  -1.58     0.670  -0.929   -0.105    11 SWE    
#> 6     5 11810 Agrar… 0.0500 FALSE -0.228   -0.205  -0.983   -0.264    11 SWE    
#> # ℹ 10 more variables: pelection <date>, n <dbl>, dur_wkb <dbl>,
#> #   event_wkb <dbl>, comp_early <dbl>, comp_replace <dbl>, majority <dbl>,
#> #   mwc <dbl>, rile_SD <dbl>, investiture <dbl>
table(coalgov$country)
#> 
#> AUS AUT BEL BGR CHE CZE DEU DNK ESP EST FIN FRA GBR GRC HRV HUN IRL ISR ITA JPN 
#>  46  45 150   7 289  31  57  76   4  27 186 182   2  13  21  19  34 252 196  58 
#> LTU LVA NLD NOR POL PRT ROU SVK SWE 
#>  41  72  74  40  42  10  43  32  28 

# Government statistics
length(unique(coalgov$gid))
#> [1] 628
mean(coalgov$event_wkb, na.rm = TRUE)
#> [1] 0.3957631
summary(coalgov$dur_wkb)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     7.0   236.0   370.0   526.8   724.0  1840.0 

# Party participation patterns
table(table(coalgov$pid))
#> 
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
#> 71 61 38 36 15  8  9  9 10  5  7  4  1  2  1  3  2  1  1  4  1  2  3  1  3  1 
#> 27 29 32 33 35 39 40 41 66 71 72 
#>  1  1  2  1  1  1  1  1  1  1  2 

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
