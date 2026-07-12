# ================================================================================================ #
# Temporal helpers: spell table, decomposition adding-up, turnover
# ================================================================================================ #

toy_panel <- data.frame(
  occ  = rep(c("a", "b"), each = 6),
  task = rep(c("t1", "t2", "t3"), 4),
  yr   = rep(rep(c(2022, 2023), each = 3), 2),
  imp  = c(1, 1, 2, 2, 1, 1, 3, 1, 1, 1, 1, 3),
  xai  = c(0.9, 0.1, 0.5, 0.9, 0.1, 0.5, 0.2, 0.4, 0.9, 0.2, 0.4, 0.9)
)

test_that("bml_delta: components sum to the total change (adding-up identity)", {
  d <- bml_delta(toy_panel, "occ", "task", "yr", "imp", "xai", 2022, 2023)
  expect_equal(d$dA, d$realloc + d$attr_change + d$entry_exit, tolerance = 1e-12)
  # attributes are constant here, so the change is pure reallocation
  expect_equal(d$attr_change, c(0, 0))
  expect_equal(d$entry_exit, c(0, 0))
  expect_equal(d$dA, d$realloc)
})

test_that("bml_realloc: direction and intensity", {
  r <- bml_realloc(toy_panel, "occ", "task", "yr", "imp", "xai", 2022, 2023)
  # occ a: weight moved toward t1 (x = 0.9): positive reallocation
  expect_gt(r$realloc_x[r$occ == "a"], 0)
  expect_true(all(r$realloc_intensity >= 0))
})

test_that("bml_turnover: entry/exit mass and jaccard", {
  # add turnover: in occ "a", t3 replaced by t4 in 2023
  tp <- toy_panel
  tp$task[tp$occ == "a" & tp$yr == 2023 & tp$task == "t3"] <- "t4"
  tv <- bml_turnover(tp, "occ", "task", "yr", "imp", t1 = 2022, t2 = 2023)

  a <- tv[tv$occ == "a", ]
  expect_gt(a$entry_mass, 0)
  expect_gt(a$exit_mass, 0)
  expect_equal(a$jaccard_dist, 1 - 2 / 4)  # 2 stayers of 4 union members
  b <- tv[tv$occ == "b", ]
  expect_equal(b$jaccard_dist, 0)

  # decomposition still adds up under entry/exit
  d <- bml_delta(tp, "occ", "task", "yr", "imp", "xai", 2022, 2023)
  expect_equal(d$dA, d$realloc + d$attr_change + d$entry_exit, tolerance = 1e-12)
})

test_that("weights normalize per snapshot: dw sums to zero", {
  s <- bml:::build_spell_table(toy_panel, "occ", "task", "yr", "imp", "xai", 2022, 2023)
  sums <- tapply(s$dw, s$group, sum)
  expect_true(all(abs(sums) < 1e-12))
})

test_that("equal weights when weight = NULL; missing columns error", {
  d <- bml_delta(toy_panel, "occ", "task", "yr", weight = NULL, x = "xai", 2022, 2023)
  expect_equal(nrow(d), 2)
  expect_error(bml_delta(toy_panel, "occ", "task", "yr", "nope", "xai", 2022, 2023),
               "not found")
  expect_error(bml_delta(toy_panel, "occ", "task", "yr", "imp", "xai", 2022, 1999),
               "No rows")
})
