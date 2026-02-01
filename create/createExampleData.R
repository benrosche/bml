# ================================================================================================ #
# Title:   A Multilevel Model for Coalition Governments. Uncovering Party-Level Dependencies Within
#          and Between Governments.
# Author:  Benjamin Rosche (Feb 2026)
# Purpose: This file generates the example data.
# ================================================================================================ #

rm(list = ls())

# ================================================================================================ #
# Load packages and data
# ================================================================================================ #

library(dplyr)
library(tidyr)
library(here)

load(
  "C:/Users/fr2488/OneDrive/GitHub/govsurvival/data/govsurvival.RData"
)

# ================================================================================================ #
# Select relevant variables and save
# ================================================================================================ #

coalgov <-
  dat.l1 |>
  select(
    gid,
    pid,
    pname=partyname,
    pseat,
    prime,
    rile,
    cohesion,
    finance = finance1,
    Nmembers
  ) |>
  left_join(
    dat.l2 |>
      select(
        cid,
        country,
        gid,
        pelection,
        n = gparties,
        dur_wkb,
        event_wkb,
        comp_early,
        comp_replace,
        majority,
        mwc,
        rile_SD = rile.gov_SD,
        investiture
      ),
    by = c("gid")
  )

  save(coalgov, file = "data/coalgov.rda", compress = "xz")
