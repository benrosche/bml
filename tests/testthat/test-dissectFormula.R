# ================================================================================================ #
# dissectFormula: named blocks, interactions, one-parameter rule, validations
# ================================================================================================ #

data("coalgov")

df_gauss <- function(formula) {
  suppressMessages(dissectFormula(formula, "Gaussian", coalgov))
}

test_that("basic dissection: lhs, mainvars, mm, hm", {
  fp <- df_gauss(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE) +
      hm(id = id(cid))
  )
  expect_equal(fp$lhs, "dur_wkb")
  expect_equal(fp$mainvars, c("X0", "majority"))
  expect_length(fp$mm, 1)
  expect_length(fp$hm, 1)
  expect_length(fp$interactions, 0)
})

test_that("auto-names: A_ for sum, V_ for var, H_ for member-paired, C_w for hhi", {
  fp <- df_gauss(
    dur_wkb ~ 1 +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE) +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("var")) +
      mm(id = id(pid, gid), vars = vars(cohesion:pseat), w = w(~ 1/n), RE = FALSE) +
      mm(id = id(pid, gid), w = w(~ pseat), fn = fn("hhi"), RE = FALSE)
  )
  expect_equal(fp$mm[[1]]$feature_labels, "A_cohesion")
  expect_equal(fp$mm[[1]]$name, "A_cohesion")
  expect_equal(fp$mm[[2]]$feature_labels, "V_cohesion")
  expect_equal(fp$mm[[3]]$feature_labels, "H_cohesion_pseat")
  expect_equal(fp$mm[[4]]$feature_labels, "C_w")
})

test_that("one parameter rule: free symbols in w() become parameters, message emitted", {
  expect_message(
    fp <- dissectFormula(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion),
           w = w(~ ilogit(b0 + b1 * pseat)), RE = FALSE),
      "Gaussian", coalgov
    ),
    "Estimating parameters: b0, b1 \\(w\\.1\\)"
  )
  expect_setequal(fp$mm[[1]]$w$params, c("b0", "b1"))
  expect_true("pseat" %in% fp$mm[[1]]$w$vars)
})

test_that("a typo'd column name becomes a parameter and shows in the message", {
  expect_message(
    dissectFormula(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion),
           w = w(~ ilogit(b0 + b1 * pseet)), RE = FALSE),  # pseet: typo
      "Gaussian", coalgov
    ),
    "pseet"
  )
})

test_that("fn shape parameters are resolved with data-driven defaults", {
  fp <- suppressMessages(df_gauss(
    dur_wkb ~ 1 +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
         fn = fn("threshold", c = est(), kappa = 10))
  ))
  ep <- fp$mm[[1]]$fn$est_params
  expect_named(ep, "c")
  expect_match(ep$c$default, "^dunif\\(")
  expect_equal(fp$mm[[1]]$fn$fixed_params$kappa, 10)
})

test_that("DSL: attrs from vars(), free symbols become parameters", {
  # kappa multiplies the attribute inside E(): the confounding warning fires by design
  expect_warning(
    fp <- suppressMessages(df_gauss(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
           fn = fn(~ E(ilogit(kappa * (cohesion - cc)))))
    )),
    "confounded"
  )
  fnobj <- fp$mm[[1]]$fn
  expect_equal(fnobj$attrs, "cohesion")
  expect_setequal(names(fnobj$est_params), c("kappa", "cc"))
})

test_that("DSL: top-level member symbol errors; multiplicative param x attr warns", {
  expect_error(
    df_gauss(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
           fn = fn(~ E(cohesion) + cohesion))
    ),
    "group-level scalar"
  )
  expect_warning(
    suppressMessages(df_gauss(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
           fn = fn(~ E(bb * cohesion + exp(cohesion))))
    )),
    "confounded"
  )
})

test_that("fn arity is validated against the design columns", {
  expect_error(
    df_gauss(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion + pseat), w = w(~ 1/n), fn = fn("var"))
    ),
    "exactly 1"
  )
  expect_error(
    df_gauss(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("cov"))
    ),
    "exactly 2"
  )
})

test_that("gmean requires strictly positive attributes", {
  expect_error(
    df_gauss(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
           fn = fn("gmean", p = 2))
    ),
    "strictly positive"
  )
})

# ------------------------------------------------------------------------------------------------ #
# Named blocks and interactions
# ------------------------------------------------------------------------------------------------ #

test_that("cross-level interaction is routed; macro var must be a main effect", {
  fp <- df_gauss(
    dur_wkb ~ 1 + majority + Ax:majority +
      mm(name = Ax, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE)
  )
  expect_length(fp$interactions, 1)
  expect_equal(fp$interactions[[1]]$type, "macro")
  expect_equal(fp$interactions[[1]]$label, "Ax:majority")
  expect_false("Ax:majority" %in% fp$mainvars)

  expect_error(
    df_gauss(
      dur_wkb ~ 1 + Ax:majority +
        mm(name = Ax, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE)
    ),
    "main effect"
  )
})

test_that("within-group-varying moderator warns (member-paired hint)", {
  expect_warning(
    df_gauss(
      dur_wkb ~ 1 + pseat + Ax:pseat +
        mm(name = Ax, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE)
    ),
    "varies within"
  )
})

test_that("block x block interactions work", {
  fp <- df_gauss(
    dur_wkb ~ 1 + Ax:Vx +
      mm(name = Ax, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE) +
      mm(name = Vx, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("var"))
  )
  expect_length(fp$interactions, 1)
  expect_equal(fp$interactions[[1]]$type, "block")
})

test_that("bare block reference errors; name-column collision errors; duplicates error", {
  expect_error(
    df_gauss(
      dur_wkb ~ 1 + Ax +
        mm(name = Ax, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE)
    ),
    "referenced bare"
  )
  expect_error(
    df_gauss(
      dur_wkb ~ 1 +
        mm(name = majority, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE)
    ),
    "collides"
  )
  expect_error(
    df_gauss(
      dur_wkb ~ 1 +
        mm(name = Ax, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE) +
        mm(name = Ax, id = id(pid, gid), vars = vars(pseat), w = w(~ 1/n), RE = FALSE)
    ),
    "Duplicate"
  )
})

test_that("multi-feature blocks cannot be referenced in interactions", {
  expect_error(
    df_gauss(
      dur_wkb ~ 1 + majority + Axy:majority +
        mm(name = Axy, id = id(pid, gid), vars = vars(cohesion + pseat), w = w(~ 1/n),
           RE = FALSE)
    ),
    "single-feature"
  )
})

# ------------------------------------------------------------------------------------------------ #
# Effects-grammar validations
# ------------------------------------------------------------------------------------------------ #

test_that("mm re() slope without a mean effect in vars() warns", {
  expect_warning(
    df_gauss(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = re(1 + pseat))
    ),
    "no mean effect"
  )
})

test_that("hm slope without a main-formula fixed effect warns", {
  expect_warning(
    df_gauss(
      dur_wkb ~ 1 + hm(id = id(cid), RE = re(1 + majority))
    ),
    "no main-formula fixed effect"
  )
})

test_that("weight/vars overlap warns about the tilted-mean near-equivalence", {
  expect_warning(
    suppressMessages(df_gauss(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion),
           w = w(~ exp(om * cohesion)), RE = FALSE)
    )),
    "smooth-max"
  )
})

test_that("only one RE block per mmid group", {
  expect_error(
    df_gauss(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE) +
        mm(id = id(pid, gid), vars = vars(pseat), w = w(~ 1/n), RE = TRUE)
    ),
    "one mm\\(\\) block per mmid"
  )
})

test_that("missing variables and different mainids are caught", {
  expect_error(
    df_gauss(dur_wkb ~ 1 + mm(id = id(pid, nope), vars = vars(cohesion), w = w(~ 1/n))),
    "not found"
  )
  expect_error(
    df_gauss(dur_wkb ~ 1 + mm(id = id(pid, gid), vars = vars(nope), w = w(~ 1/n))),
    "not found"
  )
  expect_error(
    df_gauss(
      dur_wkb ~ 1 +
        mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE) +
        mm(id = id(pid, cid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE)
    ),
    "same mainid"
  )
})

test_that("family/LHS validation", {
  expect_error(
    dissectFormula(Surv(dur_wkb, event_wkb) ~ 1, "Gaussian", coalgov),
    "one variable"
  )
  expect_error(
    dissectFormula(dur_wkb ~ 1, "Weibull", coalgov),
    "two variables"
  )
})

# ------------------------------------------------------------------------------------------------ #
# Time-indexed AR walk validation
# ------------------------------------------------------------------------------------------------ #

test_that("re(ar = time): time column must exist, be numeric, and be unique within member", {
  expect_error(
    df_gauss(
      dur_wkb ~ 1 + mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
                       RE = re(1, ar = nope))
    ),
    "not found"
  )
  expect_error(
    df_gauss(
      dur_wkb ~ 1 + mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
                       RE = re(1, ar = cname))
    ),
    "numeric"
  )
  # majority repeats within a party across governments -> duplicate times
  expect_error(
    df_gauss(
      dur_wkb ~ 1 + mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
                       RE = re(1, ar = majority))
    ),
    "duplicated time"
  )
  # gid is unique within party: passes validation
  fp <- df_gauss(
    dur_wkb ~ 1 + mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
                     RE = re(1, ar = gid))
  )
  expect_equal(fp$mm[[1]]$RE$ar$time, "gid")
})
