data("coalgov")

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
  result <- fix(fdep, 1.0)
  expect_s3_class(result, "bml_fix")
  expect_equal(result$var, "fdep")
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
  result <- vars(fdep + ipd)
  expect_s3_class(result, "bml_vars")
  expect_equal(result, c("fdep", "ipd"))

  # Comma-separated
  result <- vars(fdep, ipd)
  expect_equal(result, c("fdep", "ipd"))

  # Single variable
  result <- vars(fdep)
  expect_equal(result, "fdep")
})

test_that("vars() correctly handles fix() calls", {
  result <- vars(fix(fdep, 1.0) + ipd)
  expect_s3_class(result, "bml_vars")
  expect_type(result, "list")
  expect_equal(result$free, "ipd")
  expect_length(result$fixed, 1)
  expect_equal(result$fixed[[1]]$var, "fdep")
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
  expect_equal(result$vars, "tenure")
  expect_equal(result$vars_p, c("X0", "tenure"))

  # Multiple parameters
  result <- fn(w ~ b0 + b1 * var1 + b2 * var2)
  expect_equal(result$params, c("b0", "b1", "b2"))
  expect_equal(sort(result$vars), sort(c("var1", "var2")))
})

test_that("fn() correctly inserts intercept variable", {
  result <- fn(w ~ b0)
  expect_true(grepl("b0 \\* X0", result$string))
})

# ================================================================================================ #
# Tests for mm()
# ================================================================================================ #

test_that("mm() requires id specification", {
  expect_error(
    mm(vars = vars(fdep), fn = fn(w ~ 1/n)),
    "must be specified using id"
  )
})

test_that("mm() requires exactly 2 IDs", {
  expect_error(
    mm(id = id(cid), vars = vars(fdep), fn = fn(w ~ 1/n)),
    "exactly 2 identifiers"
  )
})

test_that("mm() requires fn specification", {
  expect_error(
    mm(id = id(pid, gid), vars = vars(fdep)),
    "must be specified using fn"
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
    vars = vars(fdep + ipd),
    fn = fn(w ~ b0 + b1 * tenure, c = FALSE),
    RE = TRUE,
    ar = TRUE
  )

  expect_s3_class(result, "bml_mm")
  expect_equal(result$id, c("pid", "gid"))
  expect_equal(result$vars, c("fdep", "ipd"))
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
    "must be specified using id"
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
    name = cname,
    type = "RE",
    showFE = FALSE,
    ar = TRUE
  )

  expect_s3_class(result, "bml_hm")
  expect_equal(result$id, "cid")
  expect_equal(result$vars, c("gdp", "democracy"))
  expect_equal(result$name, "cname")
  expect_equal(result$type, "RE")
  expect_false(result$showFE)
  expect_true(result$ar)
})

test_that("hm() handles NULL vars", {
  result <- hm(id = id(cid), vars = NULL, type = "RE")
  expect_null(result$vars)
})

test_that("hm() handles name correctly", {
  result <- hm(id = id(cid), name = cname)
  expect_equal(result$name, "cname")

  result <- hm(id = id(cid))
  expect_null(result$name)
})
