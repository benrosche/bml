# ================================================================================================ #
# Tests for broom-style methods (tidy, glance) and the functions built on them
# (bmlCompare, coefPlot)
# ================================================================================================ #

# ------------------------------------------------------------------------------------------------ #
# tidy()
# ------------------------------------------------------------------------------------------------ #

test_that("tidy() returns broom-standard columns", {
  m <- create_test_model(monitor = FALSE)

  td <- tidy(m)

  expect_s3_class(td, "tbl_df")
  expect_named(td, c("term", "estimate", "std.error", "conf.low", "conf.high", "component"))
  expect_gt(nrow(td), 0)
})

test_that("tidy() estimates match reg.table", {
  m <- create_test_model(monitor = FALSE)

  td <- tidy(m, component = "all")

  expect_equal(td$term, m$reg.table$Parameter)
  expect_equal(td$estimate, m$reg.table$mean)
  expect_equal(td$conf.low, m$reg.table$lb)
  expect_equal(td$conf.high, m$reg.table$ub)
})

test_that("tidy() component filtering works", {
  m <- create_test_model(monitor = FALSE)

  fixed <- tidy(m, component = "fixed")
  variance <- tidy(m, component = "variance")
  all <- tidy(m, component = "all")

  expect_true(all(fixed$component == "fixed"))
  expect_false(any(grepl("sigma|shape", fixed$term)))
  expect_true(all(variance$component == "variance"))
  expect_gte(nrow(all), nrow(fixed) + nrow(variance))
})

test_that("tidy() with conf.int = FALSE drops interval columns", {
  m <- create_test_model(monitor = FALSE)

  td <- tidy(m, conf.int = FALSE)

  expect_false(any(c("conf.low", "conf.high") %in% names(td)))
})

test_that("tidy() with non-default conf.level requires monitor = TRUE", {
  m <- create_test_model(monitor = FALSE)

  expect_error(
    tidy(m, conf.level = 0.9),
    "monitor = TRUE"
  )
})

test_that("tidy() with non-default conf.level narrows intervals", {
  m <- create_test_model(monitor = TRUE)

  td95 <- tidy(m)
  td80 <- tidy(m, conf.level = 0.8)

  expect_true(all(td80$conf.low >= td95$conf.low))
  expect_true(all(td80$conf.high <= td95$conf.high))
})

# ------------------------------------------------------------------------------------------------ #
# glance()
# ------------------------------------------------------------------------------------------------ #

test_that("glance() returns one-row model summary", {
  m <- create_test_model(monitor = FALSE)

  g <- glance(m)

  expect_s3_class(g, "tbl_df")
  expect_equal(nrow(g), 1)
  expect_true(all(c("nobs", "n.members", "DIC", "family", "n.iter", "n.chains") %in% names(g)))
  expect_true(is.finite(g$DIC))
  expect_gt(g$nobs, 0)
})

# ------------------------------------------------------------------------------------------------ #
# bmlCompare()
# ------------------------------------------------------------------------------------------------ #

test_that("bmlCompare() stacks models with N and DIC", {
  m1 <- create_test_model(monitor = FALSE)
  m2 <- create_test_model(monitor = FALSE)

  cmp <- bmlCompare(m1, m2)

  expect_named(cmp, c("model", "term", "estimate", "conf.low", "conf.high", "N", "DIC"))
  expect_setequal(unique(cmp$model), c("m1", "m2"))
  expect_true(all(is.finite(cmp$DIC)))
})

test_that("bmlCompare() respects user-supplied model names", {
  m1 <- create_test_model(monitor = FALSE)
  m2 <- create_test_model(monitor = FALSE)

  cmp <- bmlCompare(base = m1, extended = m2)

  expect_setequal(unique(cmp$model), c("base", "extended"))
})

test_that("bmlCompare() rejects non-bml objects", {
  expect_error(bmlCompare(lm(speed ~ dist, cars)), "bml objects")
  expect_error(bmlCompare(), "at least one")
})

# ------------------------------------------------------------------------------------------------ #
# coefPlot()
# ------------------------------------------------------------------------------------------------ #

test_that("coefPlot() returns a ggplot for a single model", {
  m <- create_test_model(monitor = FALSE)

  p <- coefPlot(m)

  expect_s3_class(p, "ggplot")
  expect_true("Intercept" %in% levels(p$data$term))
})

test_that("coefPlot() drops intercept and filters parameters", {
  m <- create_test_model(monitor = FALSE)

  p <- coefPlot(m, intercept = FALSE)
  expect_false("Intercept" %in% as.character(p$data$term))

  p2 <- coefPlot(m, parameters = "majority")
  expect_equal(as.character(unique(p2$data$term)), "majority")

  expect_error(coefPlot(m, parameters = "nonexistent"), "not found")
})

test_that("coefPlot() compares multiple models side by side", {
  m1 <- create_test_model(monitor = FALSE)
  m2 <- create_test_model(monitor = FALSE)

  p <- coefPlot(m1, m2)

  expect_s3_class(p, "ggplot")
  expect_setequal(unique(as.character(p$data$model)), c("m1", "m2"))

  # Point ranges are dodged so coefficients sit next to each other
  pr_layer <- p$layers[[which(sapply(p$layers, \(l) inherits(l$geom, "GeomPointrange")))]]
  expect_s3_class(pr_layer$position, "PositionDodge")
})
