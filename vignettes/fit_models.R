# =============================================================================
# Fit the models used in examples.Rmd ONCE and cache them as .rds.
#
# Why: the vignette is precompiled (examples.Rmd.orig -> examples.Rmd), but
# re-running precompile.R refits every JAGS model from scratch (~25 min). To make
# cosmetic edits fast, the models are fit here once and saved to vignettes/.fits/;
# examples.Rmd.orig then just readRDS() them, so precompile.R runs in seconds.
#
# Re-run this script only when a model's specification, data, or seed changes.
# Requires JAGS + rjags, plus two datasets that do not ship with the package and
# should live OUTSIDE this repository:
#   * the CILS4EU feathers used by Example 1
#   * the "_dwa" parquets written by build_onet_ai_dataset.py, used by Examples 4-6
# Both paths are required env vars (no defaults, by design -- see the notes at
# each read). Set them once in ~/.Renviron:
#   BML_CILS_DIR=/abs/path/to/cils
#   ONET_DATA_DIR=/abs/path/to/data/processed_data
#
# Run from the repository root:  Rscript vignettes/fit_models.R
# The .fits/ directory is git-ignored; the committed artdefact is examples.Rmd.
# =============================================================================

library(dplyr)
library(spData)
library(sf)
library(spdep)
library(bml)

fits_dir <- "vignettes/.fits"
dir.create(fits_dir, showWarnings = FALSE, recursive = TRUE)

# Shrink a fitted bml object to only what the vignette needs before caching:
# keep the reported parameters' MCMC draws (so monetPlot/mcmcDiag/varDecomp still
# work) and drop the thousands of latent-node draws + unused components. This lets
# every model keep monitor = TRUE while the .rds files stay well under 1 MB.
slim_fit <- function(m) {
  m$pred  <- NULL
  m$re.mm <- NULL
  m$re.hm <- NULL
  if (!is.null(m$jags.out)) {
    m$jags.out$model <- NULL   # external pointer, invalid after reload anyway
    bo   <- m$jags.out$BUGSoutput
    keep <- intersect(c(rownames(m$reg.table), "deviance"), dimnames(bo$sims.array)[[3]])
    base <- unique(sub("\\[.*$", "", keep))
    bo$sims.array <- bo$sims.array[, , keep, drop = FALSE]
    if (!is.null(bo$sims.matrix))
      bo$sims.matrix <- bo$sims.matrix[, intersect(keep, colnames(bo$sims.matrix)), drop = FALSE]
    if (!is.null(bo$sims.list))
      bo$sims.list <- bo$sims.list[names(bo$sims.list) %in% base]
    for (fld in c("mean", "sd", "median"))
      if (!is.null(bo[[fld]])) bo[[fld]] <- bo[[fld]][names(bo[[fld]]) %in% base]
    if (!is.null(bo$summary))
      bo$summary <- bo$summary[rownames(bo$summary) %in% keep, , drop = FALSE]
    bo$last.values <- NULL
    m$jags.out$BUGSoutput <- bo
  }
  m
}

save_fit <- function(obj, name) saveRDS(slim_fit(obj), file.path(fits_dir, paste0(name, ".rds")))

# -----------------------------------------------------------------------------
# Example 1: CILS4EU friendship network (data not shipped)
#
# The feathers are read in an ISOLATED R process and cached as cils_raw.rds:
# the arrow native library can segfault when loaded in the same session as
# sf/GEOS. Delete cils_raw.rds (or the whole .fits/ dir) to force a re-read.
# -----------------------------------------------------------------------------

# No default: the old fallback was "data", which resolves to the package's own
# data/ directory when this is run from the repo root. That is a tracked, shipped
# directory and restricted microdata must never land there, so require the env var.
cils_dir <- Sys.getenv("BML_CILS_DIR", unset = "")
if (!nzchar(cils_dir)) {
  stop("BML_CILS_DIR is not set. Point it at the folder holding nodedat.feather, ",
       "edgedat.feather and nodedat-w2.feather (outside this repository).")
}
cils_rds <- file.path(fits_dir, "cils_raw.rds")

if (!file.exists(cils_rds)) {
  conv <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf("cils_dir <- %s", deparse(cils_dir)),
    "saveRDS(list(",
    '  nodedat    = arrow::read_feather(file.path(cils_dir, "nodedat.feather")),',
    '  edgedat    = arrow::read_feather(file.path(cils_dir, "edgedat.feather")),',
    '  nodedat_w2 = arrow::read_feather(file.path(cils_dir, "nodedat-w2.feather"))',
    sprintf("), %s)", deparse(cils_rds))
  ), conv)
  rscript <- file.path(R.home("bin"),
                       if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (system2(rscript, shQuote(conv)) != 0) {
    stop("Failed to read CILS feathers in an isolated process. Check BML_CILS_DIR.")
  }
}

cils       <- readRDS(cils_rds)
nodedat    <- cils$nodedat
edgedat    <- cils$edgedat
nodedat_w2 <- cils$nodedat_w2

edgedat <- edgedat |>
  mutate(rank = as.integer(rank))

w2_outcome <- nodedat_w2 |> select(youthid, smoking_w2 = student_smoking)
ego_vars   <- nodedat    |> select(youthid, smoking_ego = student_smoking,
                                            gender_ego  = student_gender)
alter_vars <- nodedat    |> select(youthid, smoking_alter = student_smoking)

cils_long <-
  edgedat |>
  inner_join(ego_vars,   by = c("youthid_from" = "youthid")) |>
  inner_join(alter_vars, by = c("youthid_to"   = "youthid")) |>
  inner_join(w2_outcome, by = c("youthid_from" = "youthid")) |>
  group_by(youthid_from) |> mutate(n = n()) |> ungroup() |>
  rename(ego = youthid_from, alter = youthid_to)

mod_lim <-
  bml(
    smoking_w2 ~ smoking_ego + gender_ego +
      mm(
        id   = id(alter, ego),
        vars = bml::vars(smoking_alter),
        w    = w(~ 1/n, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = gaussian(),
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = cils_long
  )
save_fit(mod_lim, "mod_lim")

mod_bf <-
  bml(
    smoking_w2 ~ smoking_ego + gender_ego +
      mm(
        id   = id(alter, ego),
        vars = bml::vars(smoking_alter),
        w    = w(~ rank == min(rank), scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = gaussian(),
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = cils_long
  )
save_fit(mod_bf, "mod_bf")

# Third weighting rule: closeness-graded, keeps every friend (vignette Model 2).
mod_decay <-
  bml(
    smoking_w2 ~ smoking_ego + gender_ego +
      mm(
        id   = id(alter, ego),
        vars = bml::vars(smoking_alter),
        w    = w(~ 1/rank, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = gaussian(),
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = cils_long
  )
save_fit(mod_decay, "mod_decay")

# Cross-validation and posterior-predictive artifacts for the vignette. These need
# the full draws, which slim_fit strips, so compute them here and cache the small
# results: loo objects (~n x 3 pointwise) and a 30-draw yrep matrix.
# (save_fit() does not mutate its argument, so the in-memory fits are still full.)
saveRDS(loo(mod_lim),   file.path(fits_dir, "loo_lim.rds"))
saveRDS(loo(mod_bf),    file.path(fits_dir, "loo_bf.rds"))
saveRDS(loo(mod_decay), file.path(fits_dir, "loo_decay.rds"))
saveRDS(
  list(y = mod_lim$input$y, yrep = posterior_predict(mod_lim, ndraws = 30)),
  file.path(fits_dir, "ppc_lim.rds")
)

# -----------------------------------------------------------------------------
# Example 2: Boston tracts (spData; reconstructed live in the vignette too)
#
# NOTE: examples.Rmd.orig rebuilds this exact frame in a hidden chunk. Any change
# here must be mirrored there or the two drift apart.
# -----------------------------------------------------------------------------

boston <-
  read_sf(system.file("shapes/boston_tracts.gpkg", package = "spData")) |>
  select(CMEDV, NOX, CRIM, RM, DIS, AGE, LSTAT, geom) |>
  st_transform(crs = 5070) |>
  mutate(
    tid     = row_number(),
    lnCMEDV = log(CMEDV),
    across(c(NOX, CRIM, RM, DIS, AGE), ~ as.numeric(scale(.)))
  )

boston_nb <- poly2nb(as_Spatial(boston), row.names = boston$tid)

nb2df <- function(nb) {
  data.frame(
    tid    = rep(seq_along(nb), lengths(nb)),
    tid_nb = unlist(nb)
  )
}

boston_df <-
  nb2df(boston_nb) |>
  group_by(tid) |> mutate(n = n()) |> ungroup() |>
  inner_join(boston, by = "tid") |>
  inner_join(
    as.data.frame(boston) |>
      select(-CMEDV, -lnCMEDV, -geom) |>
      rename_with(~ paste0(.x, "_nb")),
    by = c("tid_nb" = "tid_nb")
  ) |>
  # aggregation variable only: LSTAT enters no formula, just the weight function
  mutate(d_LSTAT = as.numeric(scale(abs(LSTAT - LSTAT_nb))))

mod_bml <-
  bml(
    lnCMEDV ~ NOX + CRIM + RM + DIS + AGE +
      mm(
        id   = id(tid_nb, tid),
        vars = bml::vars(NOX_nb + CRIM_nb + RM_nb + DIS_nb + AGE_nb),
        w    = w(~ 1/n, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = gaussian(),
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = boston_df
  )
save_fit(mod_bml, "mod_bml")

mod_bml_lstat <-
  bml(
    lnCMEDV ~ NOX + CRIM + RM + DIS + AGE +
      mm(
        id   = id(tid_nb, tid),
        vars = bml::vars(NOX_nb + CRIM_nb + RM_nb + DIS_nb + AGE_nb),
        w    = w(~ exp(-b1 * d_LSTAT), scale = TRUE),   # nests 1/n at b1 = 0
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = gaussian(),
    prior    = prior(normal(0, 1), class = "w"),
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = boston_df
  )
save_fit(mod_bml_lstat, "mod_bml_lstat")

# -----------------------------------------------------------------------------
# Example 3: coalition governments (coalgov ships with the package)
# -----------------------------------------------------------------------------

data(coalgov)
coalgov$prime <- as.integer(coalgov$prime)

coalgov <- coalgov |>
  group_by(gid) |>
  mutate(is_smallest = as.integer(pseat == min(pseat) & pid == pid[which.min(pseat)])) |>
  ungroup()

mod_eq <-
  bml(
    Surv(dur_wkb, event_wkb) ~ 1 + majority + mwc +
      mm(
        id   = id(pid, gid),
        vars = vars(finance),
        w    = w(~ 1/n, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = weibull(),
    monitor  = TRUE,
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = coalgov
  )
save_fit(mod_eq, "mod_eq")

mod_pm_n <-
  bml(
    Surv(dur_wkb, event_wkb) ~ 1 + majority + mwc +
      mm(
        id   = id(pid, gid),
        vars = vars(finance),
        w    = w(~ b1 * prime + (1 - b1) * (1/n), scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = weibull(),
    prior    = prior(normal(0, 1), class = "w"),
    monitor  = TRUE,
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = coalgov
  )
save_fit(mod_pm_n, "mod_pm_n")

mod_pm_p <-
  bml(
    Surv(dur_wkb, event_wkb) ~ 1 + majority + mwc +
      mm(
        id   = id(pid, gid),
        vars = vars(finance),
        w    = w(~ b1 * prime + (1 - b1) * pseat, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = weibull(),
    prior    = prior(normal(0, 1), class = "w"),
    monitor  = TRUE,
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = coalgov
  )
save_fit(mod_pm_p, "mod_pm_p")

# -----------------------------------------------------------------------------
# Examples 4-6: O*NET occupation-activity structure (data not shipped)
#
# Written by build_onet_ai_dataset.py into its CONFIG$OUT_DIR. We use the "_dwa"
# build: members are detailed work activities, which recur across occupations, so
# the multiple-membership structure is present. The "_task" build puts each task
# in exactly one occupation -- no multiple membership, and mm() would add nothing
# over an ordinary occupation-level regression.
#
# TWO tables are needed. The long table carries the occupation-activity structure
# (importance weights, member exposure); the occupation-level outcomes and
# moderators live only in the model table. They are joined below.
#
# Both are read in an ISOLATED R process and cached as onet_raw.rds, for the same
# reason as the CILS feathers above: the arrow native library can segfault when
# loaded alongside sf/GEOS, and this script loads both.
# Delete onet_raw.rds (or the whole .fits/ dir) to force a re-read.
# -----------------------------------------------------------------------------

onet_dir <- Sys.getenv("ONET_DATA_DIR", unset = "")
if (!nzchar(onet_dir)) {
  stop("ONET_DATA_DIR is not set. Point it at the folder build_onet_ai_dataset.py ",
       "writes to (its CONFIG$OUT_DIR, e.g. /abs/path/to/data/processed_data).")
}

onet_long_file  <- "long_occupation_member_dwa.parquet"
onet_model_file <- "model_occupation_year_dwa.parquet"
onet_rds        <- file.path(fits_dir, "onet_raw.rds")

if (!file.exists(onet_rds)) {
  conv <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf("onet_dir <- %s", deparse(onet_dir)),
    "saveRDS(list(",
    sprintf("  long  = arrow::read_parquet(file.path(onet_dir, %s)),", deparse(onet_long_file)),
    sprintf("  model = arrow::read_parquet(file.path(onet_dir, %s))",  deparse(onet_model_file)),
    sprintf("), %s)", deparse(onet_rds))
  ), conv)
  rscript <- file.path(R.home("bin"),
                       if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (system2(rscript, shQuote(conv)) != 0) {
    stop("Failed to read O*NET parquets in an isolated process. Check ONET_DATA_DIR ",
         "and that both _dwa parquets are present.")
  }
}

onet_parts <- readRDS(onet_rds)

# The model table is one row per occupation, so this is a clean 1-to-many join
# onto the long rows. Only the columns the vignette models are carried over.
onet_raw <-
  onet_parts$long |>
  inner_join(
    onet_parts$model |>
      select(occupation_id, wage_level, employment_growth,
             education, job_zone, baseline_wage),
    by = "occupation_id"
  )

# bml's id() needs the member and group identifiers. If it rejects the character
# O*NET codes ("15-1252.00", "4.A.2.a.4.i.3"), uncomment these two lines --
# relabelling identifiers does not change the model, only how they are indexed.
# onet_raw$occupation_id <- as.integer(factor(onet_raw$occupation_id))
# onet_raw$member_id     <- as.integer(factor(onet_raw$member_id))

# Fail loudly and specifically if the column names differ from what the vignette
# prose documents, rather than erroring somewhere inside bml().
.onet_req <- c("occupation_id", "member_id", "wage_level", "employment_growth",
               "education", "job_zone", "baseline_wage",
               "ai_exposure", "importance_raw")
.onet_missing <- setdiff(.onet_req, names(onet_raw))
if (length(.onet_missing)) {
  stop("O*NET frame is missing: ", paste(.onet_missing, collapse = ", "),
       ". Check onet_file and the column names against onet_runs.Rmd.")
}

# Two complete-case frames, because Examples 4-5 and Example 6 sit on different
# outcomes.
#
# The vignette's w(~ importance_raw, scale = TRUE) reproduces the pipeline's own
# `w` column exactly: build_onet_ai_dataset.py drops members lacking exposure and
# then renormalises, and both filters below drop whole occupations rather than
# individual members, so the within-occupation normalisation is unchanged.
#
# onet    -- wage_level. ONE frame for both Example 4 and Example 5: Example 5
#            controls on education alone, but it must be estimated on the same rows
#            as Example 4 or the two sections report different N and the tables stop
#            being comparable. Target: N = 854 occupations.
# onet_wg -- employment_growth, for Example 6. Complete-cased on exactly the
#            variables mod_xl_onet uses, so the frame matches the confirmed run.
onet <- onet_raw |>
  filter(!is.na(wage_level), !is.na(education), !is.na(job_zone),
         !is.na(ai_exposure), !is.na(importance_raw))

onet_wg <- onet_raw |>
  filter(!is.na(employment_growth), !is.na(education), !is.na(baseline_wage),
         !is.na(ai_exposure), !is.na(importance_raw))

# Example 4 baseline: weighted mean only, education + job_zone
mod_mean <-
  bml(
    wage_level ~ 1 + education + job_zone +
      mm(
        id   = id(member_id, occupation_id),
        vars = bml::vars(ai_exposure),
        w    = w(~ importance_raw, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = gaussian(),
    iter   = 2000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = onet
  )
save_fit(mod_mean, "mod_mean")

# Example 4: mean + spread (sum + var features, stacked)
mod_spread <-
  bml(
    wage_level ~ 1 + education + job_zone +
      mm(
        id   = id(member_id, occupation_id),
        vars = bml::vars(ai_exposure),
        w    = w(~ importance_raw, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ) +
      mm(
        id   = id(member_id, occupation_id),
        vars = bml::vars(ai_exposure),
        w    = w(~ importance_raw, scale = TRUE),
        fn   = fn("var")
      ),
    family   = gaussian(),
    iter   = 2000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = onet
  )
save_fit(mod_spread, "mod_spread")

# Example 5 baseline: weighted mean only, education alone -- matches the controls
# of mod_smax_onet, so the two columns of the Example 5 table are comparable.
mod_mean_edu <-
  bml(
    wage_level ~ 1 + education +
      mm(
        id   = id(member_id, occupation_id),
        vars = bml::vars(ai_exposure),
        w    = w(~ importance_raw, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = gaussian(),
    iter   = 2000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = onet
  )
save_fit(mod_mean_edu, "mod_mean_edu")

# Example 5: what is the aggregation function? (smax with estimated kappa)
# monitor = TRUE is required: the vignette calls monetPlot() on fn.kappa.1.
mod_smax_onet <-
  bml(
    wage_level ~ 1 + education +
      mm(
        id   = id(member_id, occupation_id),
        vars = bml::vars(ai_exposure),
        w    = w(~ importance_raw, scale = TRUE),
        fn   = fn("smax", kappa = est())
      ) +
      mm(
        id   = id(member_id, occupation_id),
        w    = w(~ importance_raw, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = gaussian(),
    prior    = prior(normal(0, 1), class = "fn"),
    monitor  = TRUE,
    iter   = 2000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = onet
  )
save_fit(mod_smax_onet, "mod_smax_onet")

# Example 6: cross-level interaction via a named block, on employment growth.
# Kept at 5000/500 to match the three-seed confirmation the pasted numbers come
# from. monitor = TRUE because the vignette calls fixef() on this fit.
mod_xl_onet <-
  bml(
    employment_growth ~ 1 + education + baseline_wage + Aexp:education +
      mm(
        name = Aexp,
        id   = id(member_id, occupation_id),
        vars = bml::vars(ai_exposure),
        w    = w(~ importance_raw, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = gaussian(),
    monitor  = TRUE,
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = onet_wg
  )
save_fit(mod_xl_onet, "mod_xl_onet")

# loo objects for the Example 4 and 5 comparisons. As above, these need the full
# draws, so they are computed here from the un-slimmed in-memory fits.
saveRDS(loo(mod_mean),      file.path(fits_dir, "loo_mean.rds"))
saveRDS(loo(mod_spread),    file.path(fits_dir, "loo_spread.rds"))
saveRDS(loo(mod_mean_edu),  file.path(fits_dir, "loo_mean_edu.rds"))
saveRDS(loo(mod_smax_onet), file.path(fits_dir, "loo_smax_onet.rds"))

# slim_fit()'s keep-list is rownames(reg.table) + "deviance". kappa is reported in
# the fn[...] component, which may not be in reg.table, so confirm it survived the
# round-trip: monetPlot() in Example 5 reads it back off the cached object.
.m <- readRDS(file.path(fits_dir, "mod_smax_onet.rds"))
if (!"fn.kappa.1" %in% dimnames(.m$jags.out$BUGSoutput$sims.array)[[3]]) {
  stop("fn.kappa.1 was dropped by slim_fit(); extend its keep-list before precompiling.")
}
rm(.m)

# -----------------------------------------------------------------------------
# Example 7: heterogeneous member effects (explained + residual in one block)
# -----------------------------------------------------------------------------

mod_het <-
  bml(
    Surv(dur_wkb, event_wkb) ~ 1 + majority + mwc +
      mm(
        id   = id(pid, gid),
        vars = vars(finance + finance:prime),
        w    = w(~ 1/n, scale = TRUE),
        fn   = fn("sum"),
        RE   = re(1 + finance)
      ),
    family   = weibull(),
    monitor  = TRUE,
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = coalgov
  )
save_fit(mod_het, "mod_het")

message("All models fit and saved to ", normalizePath(fits_dir))
