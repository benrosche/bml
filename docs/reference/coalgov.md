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
#> gropd_df [2,077 × 19] (S3: grouped_df/tbl_df/tbl/data.frame)
#>  $ cid        : num [1:2077] 11 11 11 11 11 11 11 11 11 11 ...
#>   ..- attr(*, "format.stata")= chr "%10.0g"
#>  $ cname      : chr [1:2077] "SWE" "SWE" "SWE" "SWE" ...
#>  $ gid        : int [1:2077] 3 3 4 4 5 5 14 14 14 16 ...
#>  $ pid        : num [1:2077] 11320 11810 11320 11810 11320 ...
#>   ..- attr(*, "format.stata")= chr "%10.0g"
#>  $ pname      : chr [1:2077] "Social Democratic Labour Party" "Agrarian Party" "Social Democratic Labour Party" "Agrarian Party" ...
#>  $ election   : Date[1:2077], format: "1948-09-19" "1948-09-19" ...
#>  $ investiture: num [1:2077] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "format.stata")= chr "%10.0g"
#>  $ dur_wkb    : num [1:2077] 357 357 1466 1466 399 ...
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ event_wkb  : num [1:2077] 0 0 0 0 1 1 1 1 1 1 ...
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ n          : num [1:2077] 2 2 2 2 2 2 3 3 3 3 ...
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ majority   : num [1:2077] 0 0 0 0 0 0 1 1 1 1 ...
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ mwc        : num [1:2077] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "format.stata")= chr "%9.0g"
#>  $ rile_SD    : num [1:2077] 0.0935 0.0935 0.0935 0.0935 0.0935 ...
#>  $ pseat      : num [1:2077] 0.789 0.211 0.809 0.191 0.848 ...
#>  $ prime      : logi [1:2077] TRUE FALSE TRUE FALSE TRUE FALSE ...
#>  $ rile       : num [1:2077] -1.585 -0.228 -1.585 -0.228 -1.585 ...
#>  $ cohesion   : num [1:2077] 0.67 -0.205 0.67 -0.205 0.67 ...
#>  $ finance    : num [1:2077] -0.929 -0.983 -0.929 -0.983 -0.929 ...
#>  $ Nmembers   : num [1:2077] -0.105 -0.264 -0.105 -0.264 -0.105 ...
#>  - attr(*, "groups")= tibble [628 × 2] (S3: tbl_df/tbl/data.frame)
#>   ..$ gid  : int [1:628] 3 4 5 14 16 17 23 28 29 30 ...
#>   ..$ .rows: list<int> [1:628] 
#>   .. ..$ : int [1:2] 1 2
#>   .. ..$ : int [1:2] 3 4
#>   .. ..$ : int [1:2] 5 6
#>   .. ..$ : int [1:3] 7 8 9
#>   .. ..$ : int [1:3] 10 11 12
#>   .. ..$ : int [1:2] 13 14
#>   .. ..$ : int [1:4] 15 16 17 18
#>   .. ..$ : int [1:4] 19 20 21 22
#>   .. ..$ : int [1:4] 23 24 25 26
#>   .. ..$ : int [1:2] 27 28
#>   .. ..$ : int [1:4] 29 30 31 32
#>   .. ..$ : int [1:4] 33 34 35 36
#>   .. ..$ : int [1:4] 37 38 39 40
#>   .. ..$ : int [1:3] 41 42 43
#>   .. ..$ : int [1:3] 44 45 46
#>   .. ..$ : int [1:3] 47 48 49
#>   .. ..$ : int [1:3] 50 51 52
#>   .. ..$ : int [1:3] 53 54 55
#>   .. ..$ : int [1:3] 56 57 58
#>   .. ..$ : int [1:3] 59 60 61
#>   .. ..$ : int [1:3] 62 63 64
#>   .. ..$ : int [1:4] 65 66 67 68
#>   .. ..$ : int [1:2] 69 70
#>   .. ..$ : int [1:2] 71 72
#>   .. ..$ : int [1:2] 73 74
#>   .. ..$ : int [1:2] 75 76
#>   .. ..$ : int [1:2] 77 78
#>   .. ..$ : int [1:3] 79 80 81
#>   .. ..$ : int [1:3] 82 83 84
#>   .. ..$ : int [1:2] 85 86
#>   .. ..$ : int [1:2] 87 88
#>   .. ..$ : int [1:3] 89 90 91
#>   .. ..$ : int [1:2] 92 93
#>   .. ..$ : int [1:3] 94 95 96
#>   .. ..$ : int [1:5] 97 98 99 100 101
#>   .. ..$ : int [1:2] 102 103
#>   .. ..$ : int [1:4] 104 105 106 107
#>   .. ..$ : int [1:4] 108 109 110 111
#>   .. ..$ : int [1:4] 112 113 114 115
#>   .. ..$ : int [1:3] 116 117 118
#>   .. ..$ : int [1:2] 119 120
#>   .. ..$ : int [1:4] 121 122 123 124
#>   .. ..$ : int [1:3] 125 126 127
#>   .. ..$ : int [1:2] 128 129
#>   .. ..$ : int [1:2] 130 131
#>   .. ..$ : int [1:2] 132 133
#>   .. ..$ : int [1:2] 134 135
#>   .. ..$ : int [1:2] 136 137
#>   .. ..$ : int [1:2] 138 139
#>   .. ..$ : int [1:3] 140 141 142
#>   .. ..$ : int [1:2] 143 144
#>   .. ..$ : int [1:5] 145 146 147 148 149
#>   .. ..$ : int [1:4] 150 151 152 153
#>   .. ..$ : int [1:4] 154 155 156 157
#>   .. ..$ : int [1:4] 158 159 160 161
#>   .. ..$ : int [1:4] 162 163 164 165
#>   .. ..$ : int [1:2] 166 167
#>   .. ..$ : int [1:3] 168 169 170
#>   .. ..$ : int [1:2] 171 172
#>   .. ..$ : int [1:4] 173 174 175 176
#>   .. ..$ : int [1:3] 177 178 179
#>   .. ..$ : int [1:2] 180 181
#>   .. ..$ : int [1:3] 182 183 184
#>   .. ..$ : int [1:3] 185 186 187
#>   .. ..$ : int [1:5] 188 189 190 191 192
#>   .. ..$ : int [1:4] 193 194 195 196
#>   .. ..$ : int [1:4] 197 198 199 200
#>   .. ..$ : int [1:4] 201 202 203 204
#>   .. ..$ : int [1:4] 205 206 207 208
#>   .. ..$ : int [1:4] 209 210 211 212
#>   .. ..$ : int [1:5] 213 214 215 216 217
#>   .. ..$ : int [1:5] 218 219 220 221 222
#>   .. ..$ : int [1:5] 223 224 225 226 227
#>   .. ..$ : int [1:4] 228 229 230 231
#>   .. ..$ : int [1:4] 232 233 234 235
#>   .. ..$ : int [1:5] 236 237 238 239 240
#>   .. ..$ : int [1:5] 241 242 243 244 245
#>   .. ..$ : int [1:3] 246 247 248
#>   .. ..$ : int [1:5] 249 250 251 252 253
#>   .. ..$ : int [1:4] 254 255 256 257
#>   .. ..$ : int [1:4] 258 259 260 261
#>   .. ..$ : int [1:4] 262 263 264 265
#>   .. ..$ : int [1:3] 266 267 268
#>   .. ..$ : int [1:4] 269 270 271 272
#>   .. ..$ : int [1:4] 273 274 275 276
#>   .. ..$ : int [1:3] 277 278 279
#>   .. ..$ : int [1:4] 280 281 282 283
#>   .. ..$ : int [1:3] 284 285 286
#>   .. ..$ : int [1:5] 287 288 289 290 291
#>   .. ..$ : int [1:5] 292 293 294 295 296
#>   .. ..$ : int [1:3] 297 298 299
#>   .. ..$ : int [1:3] 300 301 302
#>   .. ..$ : int [1:4] 303 304 305 306
#>   .. ..$ : int [1:4] 307 308 309 310
#>   .. ..$ : int [1:6] 311 312 313 314 315 316
#>   .. ..$ : int [1:5] 317 318 319 320 321
#>   .. ..$ : int [1:5] 322 323 324 325 326
#>   .. ..$ : int [1:4] 327 328 329 330
#>   .. ..$ : int [1:3] 331 332 333
#>   .. .. [list output truncated]
#>   .. ..@ ptype: int(0) 
#>   ..- attr(*, ".drop")= logi TRUE
table(coalgov$cname)
#> 
#> AUS AUT BEL BGR CHE CZE DEU DNK ESP EST FIN FRA GBR GRC HRV HUN IRL ISR ITA JPN 
#>  46  45 150   7 289  31  57  76   4  27 186 182   2  13  21  19  34 252 196  58 
#> LTU LVA NLD NOR POL PRT ROU SVK SWE 
#>  41  72  74  40  42  10  43  32  28 

# Number of unique units
length(unique(coalgov$gid))   # Governments
#> [1] 628
length(unique(coalgov$pid))   # Parties
#> [1] 312
length(unique(coalgov$cid))   # Countries
#> [1] 29

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
