# Building the `examples` vignette

This note explains how the Examples vignette is built, what data it needs, and where that
data has to live. It is written for someone picking the repository up fresh.

The short version: the vignette is **precompiled**. The committed `examples.Rmd` contains no
executable code, so it renders for anyone with no data and no JAGS. You only need the data
and the setup below if you want to change a model and rebuild.

---

## 1. Why the vignette is precompiled

Fitting the fourteen models in this vignette takes about half an hour of JAGS time, and two
of the datasets cannot be redistributed. If the vignette were an ordinary `.Rmd`, every
`R CMD build`, every pkgdown run and every CRAN check would try to refit everything and fail
without the data.

Precompiling solves both problems. The models are fit once, cached, and the results are baked
into a committed artifact.

```
fit_models.R  →  vignettes/.fits/*.rds  →  examples.Rmd.orig  →  precompile.R  →  examples.Rmd
  (JAGS, ~30 min)     (git-ignored)          (the source)        (seconds)      (committed)
```

Three files do the work, all in `vignettes/`:

| File | Role | Needs data? | Needs JAGS? |
|---|---|---|---|
| `fit_models.R` | fits every model once, caches to `.fits/` | yes | yes |
| `examples.Rmd.orig` | the **source** vignette; loads cached fits, never fits anything | no | no |
| `precompile.R` | knits `.orig` into the committed `examples.Rmd` | no | no |

**The committed `examples.Rmd` is a build artifact. Never edit it by hand.** Every code block
in it is inert display text, so an edit there is silently discarded the next time anyone runs
`precompile.R`. Edit `examples.Rmd.orig` instead.

---

## 2. What data the vignette uses

Four datasets, two of which ship and two of which do not.

| Example | Data | Ships? | How it is obtained |
|---|---|---|---|
| 1 | CILS4EU Germany, Wave 1 and 2 | **no** | restricted access, see §3 |
| 2 | Boston 1970 census tracts | yes | `spData` package |
| 3, 7 | `coalgov` | yes | this package |
| 4, 5, 6 | O\*NET occupation-activity structure | **no** | built by a Python script, see §4 |

Examples 2, 3 and 7 therefore rebuild with no setup at all. Examples 1 and 4 through 6 need
the two external datasets.

---

## 3. CILS4EU (Example 1)

Three Feather files, exactly these names:

```
nodedat.feather        Wave 1 student records (smoking, gender)
edgedat.feather        friendship nominations with closeness rank
nodedat-w2.feather     Wave 2 student records (the outcome)
```

Put them in one directory **outside this repository** and point at it with an environment
variable. See §5.

CILS4EU is restricted-access data and must not be committed. `fit_models.R` has no default
path for it, by design. An earlier version fell back to `"data"`, which resolves to the
package's own `data/` directory when run from the repository root, and that directory is
tracked by git and shipped to users. The fallback was removed so a missing environment
variable fails loudly rather than quietly writing restricted data somewhere public.

---

## 4. O\*NET (Examples 4, 5 and 6)

These examples pair O\*NET's occupation-activity structure with activity-level AI exposure
scores and BLS wage and employment series. The dataset is assembled by
`build_onet_ai_dataset.py`.

### 4.1 What the vignette reads

Two Parquet files:

```
long_occupation_member_dwa.parquet     one row per occupation-activity pair
model_occupation_year_dwa.parquet      one row per occupation
```

`fit_models.R` reads both and joins them. The long table carries the membership structure,
the importance weights and each activity's AI exposure. The occupation-level outcomes and
moderators live only in the model table, so both files are required.

**Use the `_dwa` build, not `_task`.** The script writes both. At the task level each task
belongs to exactly one occupation, so there is no multiple membership and the `mm()` term
would add nothing over an ordinary occupation-level regression. At the DWA level the same
work activity recurs across occupations, which is the structure the model needs. The script's
own diagnostic reports this: `multiple_membership_present` is `true` for `_dwa` and `false`
for `_task`, and 97% of activities appear in more than one occupation.

### 4.2 Regenerating the Parquet files

Only needed if the underlying sources change. The script expects these raw files:

```
onet/Task Ratings.txt          task importance ratings, the weights
onet/Tasks to DWAs.txt         links each task to its detailed work activities
onet/Job Zones.txt             occupation preparation band
onet/Work Activities.txt       activity names, used for logging only
OpenAI/full_labelset.tsv       Eloundou et al. task-level AI exposure scores
BLS/oesm22nat/national_M2022_dl.xlsx    base-year wages and employment
BLS/oesm25nat/national_M2025_dl.xlsx    later-year wages and employment
BLS/occupation.xlsx                     employment projections, sheet "Table 1.2"
```

Set the two paths at the top of the script and run it:

```python
"RAW_DIR": "/abs/path/to/onet/raw_data",
"OUT_DIR": "/abs/path/to/onet/processed_data",
```

```bash
python build_onet_ai_dataset.py
```

It writes eight files to `OUT_DIR`, including `data_checks_dwa.json`. That file is worth
opening: it reports the counts the vignette's prose asserts, so the text can be checked
against the build without running any models.

| Vignette says | JSON key |
|---|---|
| 17,537 occupation-activity pairs | `n_member_rows` |
| 894 occupations | `n_occupations` |
| 2,080 activities | `multiple_membership.n_distinct_members` |
| about 20 activities per occupation | `multiple_membership.members_per_occupation_mean` |
| about 8 occupations per activity | `multiple_membership.occupations_per_member_mean` |
| up to 91 | `multiple_membership.occupations_per_member_max` |

`occupations_matched_to_oews` is 854, which is the N reported in Examples 4 and 5. Nothing is
lost beyond the OEWS wage match.

---

## 5. Where to put the data and how to point at it

Keep both datasets outside the repository. Nothing then depends on an ignore rule holding.

```
~/research-data/
├── cils/
│   ├── nodedat.feather
│   ├── edgedat.feather
│   └── nodedat-w2.feather
└── onet/
    ├── raw_data/          # only if regenerating the parquets
    └── processed_data/    # the two _dwa parquets live here
```

Then set two environment variables in `~/.Renviron`, which R reads at startup:

```
BML_CILS_DIR=/abs/path/to/research-data/cils
ONET_DATA_DIR=/abs/path/to/research-data/onet/processed_data
```

No quotes, no spaces around `=`, and absolute paths rather than `~`, which `arrow` does not
expand. Restart R, then confirm:

```r
Sys.getenv(c("BML_CILS_DIR", "ONET_DATA_DIR"))
```

Both are required. `fit_models.R` stops with a named message if either is missing.

Two things about that file are worth knowing in advance, because both fail silently.

`~/.Renviron` means the file in your home directory. If you create it from a shell with the
path quoted, as in `"~/.Renviron"`, the tilde does not expand and you get a literal directory
named `~` inside whatever directory you happened to be in, with the file inside it. R never
looks there. The symptom is simply that both variables come back empty, with nothing to
indicate the file exists in the wrong place. `cat /Users/yourname/.Renviron` confirms it is
where you think it is. The easiest way to avoid the problem is
`Rscript -e 'usethis::edit_r_environ()'`, which opens the right file.

R reads `.Renviron` once, at startup. Editing it in a running session has no effect, so
restart R before checking `Sys.getenv()`.

---

## 6. Requirements

**JAGS** must be installed at system level, plus the `rjags` R package.

For `fit_models.R`: `dplyr`, `spData`, `sf`, `spdep`, `bml`, `arrow`, `loo`.

For `precompile.R`, additionally: `spatialreg`, `ggplot2`, `bayesplot`, `tidyr`, `tibble`,
`posterior`, `knitr`.

`bayesplot` is easy to miss. Without it the posterior predictive figure in Example 1 renders
as an error message rather than a plot, and nothing else fails, so the build looks clean.

---

## 7. Rebuilding

```bash
cd /path/to/bml

# 1. clear stale figures. REQUIRED if any chunk was added or removed.
#    Figures are named by chunk position, so adding a chunk renames every
#    downstream image. Old files keep their old names, the new <img> tags point
#    at new ones, and you get broken images with nothing in the build log.
rm -rf vignettes/examples-figures/

# 2. fit. About 30 minutes. Needs JAGS and both environment variables.
Rscript vignettes/fit_models.R

# 3. precompile. Seconds. No JAGS, no data.
Rscript -e 'source("vignettes/precompile.R")'

# 4. render, to check it
Rscript -e 'rmarkdown::render("vignettes/examples.Rmd")'
```

Step 2 writes 24 files into `vignettes/.fits/`: fourteen fitted models, seven `loo` objects,
one posterior predictive sample, and two raw-data caches. It stops early with a specific
message if a required column is missing from the O\*NET tables, or if the `smax` shape
parameter did not survive caching.

If step 4 fails on pandoc from a plain shell, RStudio bundles a copy:

```bash
export RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools"
```

Knitting in RStudio with Ctrl+Shift+K works too. Note that `.Rmd` files should be knitted
rather than previewed through Quarto, since this is an R Markdown vignette and `quarto
preview` rejects it.

Only steps 1, 3 and 4 are needed for a prose-only change. Step 2 is required only when a
model's specification, data or seed changes.

---

## 8. Checking a rebuild

The prose quotes numbers that the rebuild has to reproduce. If any of these move, something
in the setup differs.

| Check | Expected |
|---|---|
| Ex 1 DIC | 18,700 / 17,577 / 16,964 |
| Ex 2 DIC | 292 / 229 |
| Ex 2 social dissimilarity `b1` | 0.592 [0.403, 0.797] |
| Ex 3 DIC | 3864 / 3868 / 3854 |
| Ex 4 exposure spread `V_ai_exposure` | -0.096 [-0.165, -0.027] |
| Ex 5 `fn[kappa]` | -1.57 [-2.73, -0.40] |
| Ex 6 `Aexp:education` | -0.009 [-0.025, 0.007] |

A few figures in the prose appear in no table, mainly random-effect standard deviations.
`vignettes/check_numbers.R` reads them straight out of the cached fits and prints them, so
they can be verified without refitting anything.

Open the rendered HTML and confirm four figures appear. A missing figure is the failure mode
that step 1 prevents and it is silent in the log.

---

## 9. What is committed and what is not

Committed: `examples.Rmd.orig`, `fit_models.R`, `precompile.R`, `examples.Rmd`, and
`examples-figures/`.

Not committed: `vignettes/.fits/`. Two independent rules keep it out, and both are needed.

- `.gitignore` has `vignettes/.fits/`
- `.Rbuildignore` has `^vignettes/\.fits`

`R CMD build` copies the whole `vignettes/` directory into the source tarball and does not
read `.gitignore`, so without the second rule the cache ships inside the tarball. That matters
because the cache is not only fitted parameters. `cils_raw.rds` and `onet_raw.rds` are full
copies of the source tables, and `ppc_lim.rds` holds the observed outcome vector for 4,002
students.

Worth running once after any change to those rules:

```bash
R CMD build . && tar -tzf bml_*.tar.gz | grep -Ei 'fits|parquet|feather|raw_data'
```

It should return nothing.

---

## 10. A known limitation

`coalgov` ships, so Examples 3 and 7 are fully reproducible by anyone. Examples 1 and 4
through 6 are not, because their data cannot be redistributed. The rendered vignette is
complete for readers, but a user cannot rerun `fit_models.R` themselves.

For the O\*NET examples this is fixable if desired. A trimmed occupation-activity table with
only the columns the models use compresses to well under 1 MB as an `.rda`, which would make
Examples 4 through 6 as reproducible as Example 3. The CILS4EU restriction is not something
the package can work around.
