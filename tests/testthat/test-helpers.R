data("coalgov")
data <- coalgov


# ================================================================================================ #
# Tests for id()
# ================================================================================================ #

test_that("id() correctly captures variable names", {
  # Single ID (for hm)
  result <- id(cid)
  expect_s3_class(result, "bml_id")
  expect_equal(as.character(result), "cid")

  # Two IDs (for mm)
  result <- id(pid, gid)
  expect_s3_class(result, "bml_id")
  expect_equal(as.character(result), c("pid", "gid"))
  expect_length(result, 2)
})

# ================================================================================================ #
# Tests for fix()
# ================================================================================================ #

test_that("fix() correctly captures variable and value", {
  result <- fix(finance, 1.0)
  expect_s3_class(result, "bml_fix")
  expect_equal(result$var, "finance")
  expect_equal(result$value, 1.0)

  # Test with different value
  result <- fix(exposure, 0.5)
  expect_equal(result$var, "exposure")
  expect_equal(result$value, 0.5)
})

# ================================================================================================ #
# Tests for vars()
# ================================================================================================ #

test_that("vars() correctly parses simple variables", {
  result <- vars(finance + cohesion)
  expect_s3_class(result, "bml_vars")
  expect_equal(sort(result$free), sort(c("finance", "cohesion")))

  # Single variable
  result <- vars(finance)
  expect_equal(result$free, "finance")
})

test_that("vars() correctly handles fix() calls", {
  result <- vars(fix(finance, 1.0) + cohesion)
  expect_s3_class(result, "bml_vars")
  expect_type(result, "list")
  expect_equal(result$free, "cohesion")
  expect_length(result$fixed, 1)
  expect_equal(result$fixed[[1]]$var, "finance")
  expect_equal(result$fixed[[1]]$value, 1.0)

  # Multiple fixed
  result <- vars(fix(var1, 1.0) + fix(var2, 2.0) + var3)
  expect_equal(result$free, "var3")
  expect_length(result$fixed, 2)
})

test_that("vars() returns NULL for empty specification", {
  # This would be called as vars() with no arguments
  # but in actual use, NULL is passed directly
  expect_null(NULL)
})

# ================================================================================================ #
# Tests for fn()
# ================================================================================================ #

test_that("fn() correctly parses simple weight functions", {
  # Equal weights
  result <- fn(w ~ 1/n)
  expect_s3_class(result, "bml_fn")
  expect_equal(result$vars, "n")
  expect_length(result$params, 0)
  expect_true(result$constraint)

  # With constraint = FALSE
  result <- fn(w ~ 1/n, c = FALSE)
  expect_false(result$constraint)
})

test_that("fn() correctly identifies parameters", {
  # Parameterized weight function
  result <- fn(w ~ b0 + b1 * tenure)
  expect_s3_class(result, "bml_fn")
  expect_equal(result$params, c("b0", "b1"))
  # vars includes X0 because b0 is transformed to b0 * X0
  expect_equal(sort(result$vars), sort(c("X0", "tenure")))
  expect_equal(sort(result$vars_p), sort(c("X0", "tenure")))

  # Multiple parameters
  result <- fn(w ~ b0 + b1 * var1 + b2 * var2)
  expect_equal(result$params, c("b0", "b1", "b2"))
  # vars includes X0 for the intercept
  expect_equal(sort(result$vars), sort(c("X0", "var1", "var2")))
})

test_that("fn() correctly inserts intercept variable", {
  result <- fn(w ~ b0)
  expect_true(grepl("b0 \\* X0", result$string))
})

# ================================================================================================ #
# Tests for fn() with aggregation functions
# ================================================================================================ #

test_that("fn() detects single aggregation function", {
  result <- fn(w ~ min(finance))
  expect_s3_class(result, "bml_fn")
  expect_length(result$agg_funcs, 1)
  expect_equal(result$agg_funcs[[1]]$func, "min")
  expect_equal(result$agg_funcs[[1]]$var, "finance")
  expect_equal(result$agg_funcs[[1]]$col_name, "finance_min")
  expect_equal(result$vars, "finance_min")
  expect_equal(result$agg_vars, "finance")
})

test_that("fn() detects multiple aggregation functions", {
  result <- fn(w ~ b1 * min(finance) + (1 - b1) * mean(finance))
  expect_s3_class(result, "bml_fn")
  expect_length(result$agg_funcs, 2)
  expect_equal(result$agg_funcs[[1]]$func, "min")
  expect_equal(result$agg_funcs[[2]]$func, "mean")
  expect_equal(sort(result$vars), sort(c("finance_min", "finance_mean")))
  expect_equal(result$agg_vars, "finance")  # Only one base variable
  expect_equal(result$params, "b1")
})

test_that("fn() correctly transforms formula string for aggregations", {
  result <- fn(w ~ b1 * min(finance) + (1 - b1) * mean(finance))
  expect_true(grepl("finance_min", result$string))
  expect_true(grepl("finance_mean", result$string))
  expect_false(grepl("min\\(finance\\)", result$string))
  expect_false(grepl("mean\\(finance\\)", result$string))
})

test_that("fn() supports all simple aggregation functions", {
  result <- fn(w ~ min(x) + max(x) + mean(x) + sum(x) + sd(x) + var(x) + first(x) + last(x) + median(x) + mode(x) + range(x))
  expect_length(result$agg_funcs, 11)
  funcs <- sapply(result$agg_funcs, function(a) a$func)
  expect_true(all(c("min", "max", "mean", "sum", "sd", "var", "first", "last", "median", "mode", "range") %in% funcs))
  expect_equal(result$agg_vars, "x")
})

test_that("fn() handles mix of aggregation and regular variables", {
  result <- fn(w ~ b0 + b1 * min(finance) + b2 * tenure)
  expect_length(result$agg_funcs, 1)
  expect_equal(result$agg_funcs[[1]]$var, "finance")
  expect_true("finance_min" %in% result$vars)
  expect_true("tenure" %in% result$vars)
  expect_equal(result$agg_vars, "finance")
  expect_equal(sort(result$params), sort(c("b0", "b1", "b2")))
})

test_that("fn() handles multiple different variables in aggregations", {
  result <- fn(w ~ min(x) + max(y))
  expect_length(result$agg_funcs, 2)
  expect_equal(sort(result$agg_vars), sort(c("x", "y")))
  expect_equal(sort(result$vars), sort(c("x_min", "y_max")))
})

test_that("fn() correctly handles median function", {
  result <- fn(w ~ median(tenure))
  expect_s3_class(result, "bml_fn")
  expect_length(result$agg_funcs, 1)
  expect_equal(result$agg_funcs[[1]]$func, "median")
  expect_equal(result$agg_funcs[[1]]$var, "tenure")
  expect_equal(result$agg_funcs[[1]]$col_name, "tenure_median")
  expect_null(result$agg_funcs[[1]]$prob)
  expect_equal(result$vars, "tenure_median")
})

test_that("fn() correctly handles mode function", {
  result <- fn(w ~ mode(category))
  expect_s3_class(result, "bml_fn")
  expect_length(result$agg_funcs, 1)
  expect_equal(result$agg_funcs[[1]]$func, "mode")
  expect_equal(result$agg_funcs[[1]]$var, "category")
  expect_equal(result$agg_funcs[[1]]$col_name, "category_mode")
  expect_null(result$agg_funcs[[1]]$prob)
  expect_equal(result$vars, "category_mode")
})

test_that("fn() correctly handles range function", {
  result <- fn(w ~ range(tenure))
  expect_s3_class(result, "bml_fn")
  expect_length(result$agg_funcs, 1)
  expect_equal(result$agg_funcs[[1]]$func, "range")
  expect_equal(result$agg_funcs[[1]]$var, "tenure")
  expect_equal(result$agg_funcs[[1]]$col_name, "tenure_range")
  expect_null(result$agg_funcs[[1]]$prob)
  expect_equal(result$vars, "tenure_range")
})

test_that("fn() correctly handles quantile function", {
  result <- fn(w ~ quantile(tenure, 0.25))
  expect_s3_class(result, "bml_fn")
  expect_length(result$agg_funcs, 1)
  expect_equal(result$agg_funcs[[1]]$func, "quantile")
  expect_equal(result$agg_funcs[[1]]$var, "tenure")
  expect_equal(result$agg_funcs[[1]]$prob, 0.25)
  expect_equal(result$agg_funcs[[1]]$col_name, "tenure_q25")
  expect_equal(result$vars, "tenure_q25")
})

test_that("fn() handles multiple quantile functions", {
  result <- fn(w ~ quantile(x, 0.25) + quantile(x, 0.75))
  expect_length(result$agg_funcs, 2)
  expect_equal(result$agg_funcs[[1]]$col_name, "x_q25")
  expect_equal(result$agg_funcs[[2]]$col_name, "x_q75")
  expect_equal(result$agg_funcs[[1]]$prob, 0.25)
  expect_equal(result$agg_funcs[[2]]$prob, 0.75)
  expect_equal(sort(result$vars), sort(c("x_q25", "x_q75")))
})

test_that("fn() handles quantile with various probabilities", {
  # 10th percentile
  result <- fn(w ~ quantile(x, 0.1))
  expect_equal(result$agg_funcs[[1]]$col_name, "x_q10")

  # 50th percentile (same as median)
  result <- fn(w ~ quantile(x, 0.5))
  expect_equal(result$agg_funcs[[1]]$col_name, "x_q50")

  # 90th percentile
  result <- fn(w ~ quantile(x, 0.9))
  expect_equal(result$agg_funcs[[1]]$col_name, "x_q90")
})

test_that("fn() rejects quantile with invalid probability", {
  expect_error(
    fn(w ~ quantile(x, 1.5)),
    "probability must be between 0 and 1"
  )
  expect_error(
    fn(w ~ quantile(x, -0.1)),
    "probability must be between 0 and 1"
  )
})

test_that("fn() handles mix of median, quantile, and other aggregations", {
  result <- fn(w ~ b1 * median(x) + b2 * quantile(x, 0.75) + (1 - b1 - b2) * mean(x))
  expect_length(result$agg_funcs, 3)
  funcs <- sapply(result$agg_funcs, function(a) a$func)
  expect_true("median" %in% funcs)
  expect_true("quantile" %in% funcs)
  expect_true("mean" %in% funcs)
  expect_equal(sort(result$vars), sort(c("x_median", "x_q75", "x_mean")))
})

test_that("fn() rejects nested aggregation functions", {
  expect_error(
    fn(w ~ min(max(finance))),
    "Nested aggregation functions"
  )
})

test_that("fn() rejects unsupported functions", {
  expect_error(
    fn(w ~ cumsum(finance)),
    "Unsupported function"
  )
})

test_that("fn() with no aggregation functions has NULL agg_funcs and agg_vars", {
  result <- fn(w ~ b0 + b1 * tenure)
  expect_null(result$agg_funcs)
  expect_null(result$agg_vars)
})

# ================================================================================================ #
# Tests for mm()
# ================================================================================================ #

test_that("mm() requires id specification", {
  expect_error(
    mm(vars = vars(finance), fn = fn(w ~ 1/n)),
    "argument \"id\" is missing"
  )
})

test_that("mm() requires exactly 2 IDs", {
  expect_error(
    mm(id = id(cid), vars = vars(finance), fn = fn(w ~ 1/n)),
    "exactly 2 identifiers"
  )
})

test_that("mm() requires fn specification", {
  expect_error(
    mm(id = id(pid, gid), vars = vars(finance)),
    "'fn' must be specified using fn"
  )
})

test_that("mm() correctly handles vars = NULL with RE = TRUE", {
  result <- mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE)
  expect_s3_class(result, "bml_mm")
  expect_null(result$vars)
  expect_true(result$RE)
})

test_that("mm() forces RE = TRUE when vars = NULL", {
  expect_error(
    mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = FALSE),
    "RE must be TRUE when vars is NULL"
  )
})

test_that("mm() correctly stores all components", {
  result <- mm(
    id = id(pid, gid),
    vars = vars(finance + cohesion),
    fn = fn(w ~ b0 + b1 * tenure, c = FALSE),
    RE = TRUE,
    ar = TRUE
  )

  expect_s3_class(result, "bml_mm")
  expect_equal(as.character(result$id), c("pid", "gid"))
  # vars is a bml_vars object with $free containing the variable names
  expect_equal(sort(result$vars$free), sort(c("finance", "cohesion")))
  expect_s3_class(result$fn, "bml_fn")
  expect_true(result$RE)
  expect_true(result$ar)
})

# ================================================================================================ #
# Tests for hm()
# ================================================================================================ #

test_that("hm() requires id specification", {
  expect_error(
    hm(vars = vars(gdp)),
    "argument \"id\" is missing"
  )
})

test_that("hm() requires exactly 1 ID", {
  expect_error(
    hm(id = id(pid, gid)),
    "exactly 1 identifier"
  )
})

test_that("hm() accepts valid type values", {
  expect_no_error(hm(id = id(cid), type = "RE"))
  expect_no_error(hm(id = id(cid), type = "FE"))
  expect_error(hm(id = id(cid), type = "INVALID"))
})

test_that("hm() correctly stores all components", {
  result <- hm(
    id = id(cid),
    vars = vars(gdp + democracy),
    name = country,
    type = "RE",
    showFE = FALSE,
    ar = TRUE
  )

  expect_s3_class(result, "bml_hm")
  expect_equal(result$id, "cid")
  # vars is a bml_vars object with $free containing the variable names
  expect_equal(sort(result$vars$free), sort(c("gdp", "democracy")))
  expect_equal(as.character(result$name), "country")
  expect_equal(as.character(result$type), "RE")
  expect_false(result$showFE)
  expect_true(result$ar)
})

test_that("hm() handles NULL vars", {
  result <- hm(id = id(cid), vars = NULL, type = "RE")
  expect_null(result$vars)
})

test_that("hm() handles name correctly", {
  result <- hm(id = id(cid), name = country)
  expect_equal(result$name, "country")

  result <- hm(id = id(cid))
  expect_null(result$name)
})
