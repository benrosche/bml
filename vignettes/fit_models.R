# =============================================================================
# Fit the models used in examples.Rmd ONCE and cache them as .rds.
#
# Why: the vignette is precompiled (examples.Rmd.orig -> examples.Rmd), but
# re-running precompile.R refits every JAGS model from scratch (~25 min). To make
# cosmetic edits fast, the models are fit here once and saved to vignettes/.fits/;
# examples.Rmd.orig then just readRDS() them, so precompile.R runs in seconds.
#
# Re-run this script only when a model's specification, data, or seed changes.
# Requires JAGS + rjags and the CILS4EU feathers used by Example 1 (not shipped):
# point at them via BML_CILS_DIR, e.g.
#   Sys.setenv(BML_CILS_DIR = "C:/path/to/cils/data")
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
cils_dir <- Sys.getenv("BML_CILS_DIR", unset = "data")
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
    warmup = 1000,
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
    warmup = 1000,
    chains = 3,
    seed     = 1,
    data     = cils_long
  )
save_fit(mod_bf, "mod_bf")

# Cross-validation and posterior-predictive artifacts for the vignette. These need
# the full draws, which slim_fit strips, so compute them here and cache the small
# results: loo objects (~n x 3 pointwise) and a 30-draw yrep matrix.
saveRDS(loo(mod_lim), file.path(fits_dir, "loo_lim.rds"))
saveRDS(loo(mod_bf),  file.path(fits_dir, "loo_bf.rds"))
saveRDS(
  list(y = mod_lim$input$y, yrep = posterior_predict(mod_lim, ndraws = 30)),
  file.path(fits_dir, "ppc_lim.rds")
)

# -----------------------------------------------------------------------------
# Example 2: Boston tracts (spData; reconstructed live in the vignette too)
# -----------------------------------------------------------------------------
boston <-
  read_sf(system.file("shapes/boston_tracts.gpkg", package = "spData")) |>
  select(CMEDV, NOX, CRIM, RM, DIS, AGE, geom) |>
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
  )

boston_df2 <- boston_df |>
  mutate(
    d_CRIM = as.numeric(scale(log1p(abs(CRIM - CRIM_nb)))),
    d_AGE  = as.numeric(scale(abs(AGE - AGE_nb)))
  )

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
    seed     = 42,
    data     = boston_df
  )
save_fit(mod_bml, "mod_bml")

mod_bml_w <-
  bml(
    lnCMEDV ~ NOX + CRIM + RM + DIS + AGE +
      mm(
        id   = id(tid_nb, tid),
        vars = bml::vars(NOX_nb + CRIM_nb + RM_nb + DIS_nb + AGE_nb),
        w    = w(~ 1 / (1 + (n - 1) * exp(-(b0 + b1 * d_CRIM + b2 * d_AGE))), scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ),
    family   = gaussian(),
    prior    = prior(normal(0, 1), class = "w"),
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 42,
    data     = boston_df2
  )
save_fit(mod_bml_w, "mod_bml_w")

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
# Example 4: mean vs spread (sum + var features, stacked)
# -----------------------------------------------------------------------------
mod_var <-
  bml(
    Surv(dur_wkb, event_wkb) ~ 1 + majority + mwc +
      mm(
        id   = id(pid, gid),
        vars = vars(finance),
        w    = w(~ 1/n, scale = TRUE),
        fn   = fn("sum"),
        RE   = TRUE
      ) +
      mm(
        id   = id(pid, gid),
        vars = vars(finance),
        w    = w(~ 1/n, scale = TRUE),
        fn   = fn("var")
      ),
    family   = weibull(),
    monitor  = TRUE,
    iter   = 5000,
    warmup = 500,
    chains = 3,
    seed     = 1,
    data     = coalgov
  )
save_fit(mod_var, "mod_var")

# -----------------------------------------------------------------------------
# Example 5: what is the aggregation function? (smax with estimated kappa)
# -----------------------------------------------------------------------------
mod_smax <-
  bml(
    Surv(dur_wkb, event_wkb) ~ 1 + majority + mwc +
      mm(
        id   = id(pid, gid),
        vars = vars(finance),
        w    = w(~ 1/n, scale = TRUE),
        fn   = fn("smax", kappa = est())
      ) +
      mm(
        id   = id(pid, gid),
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
save_fit(mod_smax, "mod_smax")

# -----------------------------------------------------------------------------
# Example 6: cross-level interaction via a named block
# -----------------------------------------------------------------------------
mod_xl <-
  bml(
    Surv(dur_wkb, event_wkb) ~ 1 + majority + mwc + Afin:majority +
      mm(
        name = Afin,
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
save_fit(mod_xl, "mod_xl")

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
