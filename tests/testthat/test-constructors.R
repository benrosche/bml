# ================================================================================================ #
# Constructors: id(), vars(), w(), fn(), est(), re(), fe(), mm(), hm()
# ================================================================================================ #

test_that("id() captures identifiers", {
  x <- id(pid, gid)
  expect_s3_class(x, "bml_id")
  expect_equal(as.character(x), c("pid", "gid"))
})

test_that("vars() supports +, commas, interactions, I(), fix()", {
  v1 <- vars(a + b)
  expect_equal(v1$free, c("a", "b"))

  v2 <- vars(a, b)
  expect_equal(v2$free, c("a", "b"))

  v3 <- vars(a:b)
  expect_equal(all.vars(v3$formula), c("a", "b"))
  expect_true(grepl("a:b", deparse(v3$formula)))

  v4 <- vars(I(a^2))
  expect_true(grepl("I\\(a\\^2\\)", deparse(v4$formula)))

  v5 <- vars(fix(off, 1.0) + a)
  expect_equal(v5$free, "a")
  expect_equal(v5$fixed[[1]]$var, "off")
  expect_equal(v5$fixed[[1]]$value, 1.0)

  expect_null(vars())
})

# ------------------------------------------------------------------------------------------------ #
# w()
# ------------------------------------------------------------------------------------------------ #

test_that("w() takes one-sided formulas and scale", {
  x <- w(~ 1 / n)
  expect_s3_class(x, "bml_w")
  expect_true(x$scale)
  expect_true("n" %in% x$symbols)

  x2 <- w(~ importance, scale = FALSE)
  expect_false(x2$scale)
})

test_that("w() rejects the old two-sided spelling with a migration message", {
  expect_error(w(w ~ 1 / n), "one-sided")
})

test_that("w() detects aggregation shortcuts and rejects nesting/unknown functions", {
  x <- w(~ b1 * min(ten) + (1 - b1) * max(ten))
  aggs <- sapply(x$agg_funcs, function(a) a$func)
  expect_setequal(aggs, c("min", "max"))
  expect_true(all(c("ten_min", "ten_max") %in% x$symbols))

  xq <- w(~ quantile(ten, 0.75))
  expect_equal(xq$agg_funcs[[1]]$col_name, "ten_q75")
  expect_error(w(~ quantile(ten, 1.5)), "between 0 and 1")

  expect_error(w(~ min(max(x))), "Nested")
  expect_error(w(~ foo(x)), "Unsupported")
})

# ------------------------------------------------------------------------------------------------ #
# fn()
# ------------------------------------------------------------------------------------------------ #

test_that("fn() named types and shape parameters", {
  expect_equal(fn("sum")$type, "sum")
  expect_equal(fn("var")$moment, 2)
  expect_equal(fn("var", moment = 3)$moment, 3)
  expect_error(fn("var", moment = 1.5), "integer")

  th <- fn("threshold", c = 0.7)
  expect_equal(th$params$c, 0.7)
  expect_equal(th$params$kappa, 10)
  expect_error(fn("threshold"), "cutpoint")

  expect_s3_class(fn("threshold", c = est())$params$c, "bml_est")
  expect_warning(fn("threshold", c = est(), kappa = est()), "weakly identified")

  expect_error(fn("smax"), "kappa")
  expect_error(fn("smax", kappa = 0), "nonzero")
  expect_error(fn("gmean", p = 0), "nonzero")
  expect_error(fn("gmean"), "p =")
  expect_equal(fn("cov")$type, "cov")

  expect_error(fn("nope"), "type string")
  expect_error(fn("sum", bogus = 1), "Unknown")
})

test_that("fn() rejects the old weights spelling with a migration message", {
  expect_error(fn(w ~ 1 / n), "[Ww]eights moved to w\\(\\)")
})

test_that("fn() DSL: whitelist, E() requirement", {
  f <- fn(~ E((x - E(x))^2))
  expect_equal(f$type, "expr")

  expect_error(fn(~ E(sqrt(x))), "non-whitelisted")
  expect_error(fn(~ x^2), "E\\(")
})

test_that("est() carries bounds", {
  e <- est(lower = 0, upper = 1)
  expect_equal(e$lower, 0)
  expect_equal(e$upper, 1)
})

# ------------------------------------------------------------------------------------------------ #
# re() / fe()
# ------------------------------------------------------------------------------------------------ #

test_that("re() parses the lme4-style LHS grammar; cor defaults to FALSE", {
  r1 <- re(1)
  expect_true(r1$intercept)
  expect_length(r1$slopes, 0)
  expect_false(r1$cor)

  r2 <- re(1 + x)
  expect_true(r2$intercept)
  expect_equal(r2$slopes, "x")
  expect_false(r2$cor)

  r3 <- re(0 + x)
  expect_false(r3$intercept)
  expect_equal(r3$slopes, "x")

  r4 <- re(x - 1)
  expect_false(r4$intercept)

  r5 <- re(1 + x, cor = TRUE)
  expect_true(r5$cor)

  expect_error(re(1, cor = TRUE), "nothing to correlate")
  expect_error(re(1 + x + y, cor = TRUE), "2x2")
})

test_that("fe() mirrors re() with showFE", {
  f1 <- fe(1, showFE = TRUE)
  expect_true(f1$intercept)
  expect_true(f1$showFE)

  f2 <- fe(1 + x)
  expect_equal(f2$slopes, "x")
  expect_false(f2$showFE)
})

# ------------------------------------------------------------------------------------------------ #
# mm()
# ------------------------------------------------------------------------------------------------ #

test_that("mm() defaults: w(~1/n), fn(sum); RE shorthands", {
  m <- mm(id = id(pid, gid), vars = vars(x))
  expect_s3_class(m$w, "bml_w")
  expect_equal(m$fn$type, "sum")
  expect_null(m$RE)

  m2 <- mm(id = id(pid, gid), vars = vars(x), RE = TRUE)
  expect_s3_class(m2$RE, "bml_re")

  m3 <- mm(id = id(pid, gid), vars = vars(x), RE = FALSE)
  expect_null(m3$RE)

  # weights-only sum block defaults to RE (variance decomposition)
  m4 <- mm(id = id(pid, gid))
  expect_s3_class(m4$RE, "bml_re")
})

test_that("mm() composition rules", {
  expect_error(mm(id = id(pid)), "exactly 2")
  expect_error(mm(id = id(pid, gid), vars = vars(x), RE = re(1), FE = fe(1)), "not both")
  expect_error(mm(id = id(pid, gid), vars = vars(x), fn = fn("var"), RE = TRUE),
               "only available with fn\\(\"sum\"\\)")
  expect_error(mm(id = id(pid, gid), vars = vars(x), fn = fn("hhi")), "weights alone")
  expect_error(mm(id = id(pid, gid), fn = fn("var")), "requires member attributes")
  expect_warning(mm(id = id(pid, gid), vars = vars(x), FE = fe(1)), "weakly identified")
  expect_error(mm(id = id(pid, gid), vars = vars(x), ar = TRUE), "moved into")
})

test_that("mm() forces scale = TRUE for emergent types", {
  expect_warning(
    m <- mm(id = id(pid, gid), vars = vars(x), w = w(~ imp, scale = FALSE), fn = fn("var")),
    "normalized"
  )
  expect_true(m$w$scale)
})

test_that("mm() captures names as symbol or string; old c= errors", {
  m1 <- mm(id = id(pid, gid), vars = vars(x), name = Ax)
  expect_equal(m1$name, "Ax")
  m2 <- mm(id = id(pid, gid), vars = vars(x), name = "Ax")
  expect_equal(m2$name, "Ax")
  expect_error(mm(id = id(pid, gid), vars = vars(x), c = TRUE), "scale")
})

# ------------------------------------------------------------------------------------------------ #
# hm()
# ------------------------------------------------------------------------------------------------ #

test_that("hm() defaults to re(1); RE/FE exclusive; removed args error", {
  h <- hm(id = id(cid))
  expect_s3_class(h$RE, "bml_re")

  h2 <- hm(id = id(cid), FE = fe(1, showFE = TRUE))
  expect_null(h2$RE)
  expect_true(h2$FE$showFE)

  expect_error(hm(id = id(cid), RE = re(1), FE = fe(1)), "not both")
  expect_error(hm(id = id(cid), type = "RE"), "redesigned")
  expect_error(hm(id = id(cid), vars = vars(gdp)), "redesigned")
  expect_error(hm(id = id(cid), showFE = TRUE), "redesigned")
  expect_error(hm(id = id(cid, gid)), "exactly 1")
})

test_that("hm() captures labels", {
  h <- hm(id = id(cid), labels = cname)
  expect_equal(h$labels, "cname")
})

# ------------------------------------------------------------------------------------------------ #
# re(ar = ...) and the ar migration
# ------------------------------------------------------------------------------------------------ #

test_that("re(ar = ...) parses TRUE / unquoted column / FALSE", {
  expect_null(re(1)$ar)
  expect_null(re(1, ar = FALSE)$ar)

  r1 <- re(1, ar = TRUE)
  expect_type(r1$ar, "list")
  expect_null(r1$ar$time)

  r2 <- re(1, ar = year)
  expect_equal(r2$ar$time, "year")

  r3 <- re(1, ar = "year")
  expect_equal(r3$ar$time, "year")
})

test_that("re(ar) composition rules", {
  expect_error(re(1 + x, ar = TRUE), "intercept-only")
  expect_error(re(0 + x, ar = TRUE), "intercept-only")
  expect_error(re(1, ar = TRUE, cor = TRUE), "nothing to correlate|cor = TRUE")
})

test_that("mm()/hm() ar argument gives a migration error; ar is carried from re()", {
  expect_error(mm(id = id(pid, gid), vars = vars(x), RE = TRUE, ar = TRUE), "moved into")
  expect_error(hm(id = id(cid), ar = TRUE), "moved into")

  m <- mm(id = id(pid, gid), vars = vars(x), RE = re(1, ar = TRUE))
  expect_type(m$ar, "list")
  m2 <- mm(id = id(pid, gid), vars = vars(x), RE = TRUE)
  expect_null(m2$ar)
  h <- hm(id = id(cid), RE = re(1, ar = yr))
  expect_equal(h$ar$time, "yr")
})
