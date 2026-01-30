data("coalgov")

# ================================================================================================ #
# Tests for formatJags() - formats JAGS output for display
# ================================================================================================ #

# Helper function to create mock JAGS output with proper BUGSoutput$summary structure
create_mock_jags_output <- function(params) {
  # Build summary matrix with columns: mean, sd, 2.5%, 25%, 50%, 75%, 97.5%, Rhat, n.eff
  n_params <- length(params)
  summary_mat <- matrix(0, nrow = n_params, ncol = 9)
  rownames(summary_mat) <- names(params)
  colnames(summary_mat) <- c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat", "n.eff")

  for (i in seq_along(params)) {
    val <- params[[i]]
    summary_mat[i, "mean"] <- val
    summary_mat[i, "sd"] <- 0.1
    summary_mat[i, "2.5%"] <- val - 0.2
    summary_mat[i, "25%"] <- val - 0.1
    summary_mat[i, "50%"] <- val
    summary_mat[i, "75%"] <- val + 0.1
    summary_mat[i, "97.5%"] <- val + 0.2
    summary_mat[i, "Rhat"] <- 1.0
    summary_mat[i, "n.eff"] <- 1000
  }

  # Add deviance row
  dev_row <- matrix(c(100, 10, 80, 90, 100, 110, 120, 1.0, 1000), nrow = 1)
  rownames(dev_row) <- "deviance"
  colnames(dev_row) <- colnames(summary_mat)
  summary_mat <- rbind(summary_mat, dev_row)

  list(
    BUGSoutput = list(
      summary = summary_mat,
      DIC = 1000
    ),
    model = "mock model"
  )
}

# Helper to create minimal main structure
create_main <- function(vars = c("X0", "majority"), lhs = "sim.y") {
  list(
    vars = vars,
    vars_fixed = NULL,
    dat = data.frame(mainid = 1:10),
    dat_fixed = NULL,
    fix_values = NULL,
    lhs = lhs,
    formula = ~ 1 + majority
  )
}

# Helper to create minimal Ns structure
create_Ns <- function(n.main = 10, n.umm = 20, mmn = rep(2, 10), n.hm = 5) {
  list(
    n.main = n.main,
    n.umm = n.umm,
    mmn = mmn,
    n.hm = n.hm,
    n.GPN = max(mmn),
    n.HMN = n.hm,
    n.mmblocks = 1
  )
}

test_that("formatJags() creates correct output structure", {
  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "sigma" = 0.8
  ))

  main <- create_main()
  Ns <- create_Ns()

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = NULL,
    main = main,
    hm_blocks = NULL,
    mm = NULL,
    hm = NULL,
    family = "Gaussian",
    cox_intervals = NULL
  )

  expect_type(result, "list")
  expect_true("reg.table" %in% names(result))
  expect_true("w" %in% names(result))
  expect_true("re.mm" %in% names(result))
  expect_true("re.hm" %in% names(result))
  expect_true("pred" %in% names(result))
})

test_that("formatJags() creates regression table with correct columns", {
  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "sigma" = 0.8
  ))

  main <- create_main()
  Ns <- create_Ns()

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = NULL,
    main = main,
    hm_blocks = NULL,
    mm = NULL,
    hm = NULL,
    family = "Gaussian",
    cox_intervals = NULL
  )

  reg_table <- result$reg.table

  expect_true("Parameter" %in% names(reg_table))
  expect_true("mean" %in% names(reg_table))
  expect_true("sd" %in% names(reg_table))
  expect_true("lb" %in% names(reg_table))
  expect_true("ub" %in% names(reg_table))
})

test_that("formatJags() correctly labels parameters", {
  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "sigma" = 0.8
  ))

  main <- create_main()
  Ns <- create_Ns()

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = NULL,
    main = main,
    hm_blocks = NULL,
    mm = NULL,
    hm = NULL,
    family = "Gaussian",
    cox_intervals = NULL
  )

  reg_table <- result$reg.table

  # Should have intercept and covariate labels
  expect_true("Intercept" %in% reg_table$Parameter)
  expect_true("majority" %in% reg_table$Parameter)
})

test_that("formatJags() includes family-specific parameters for Weibull", {
  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "shape" = 1.5
  ))

  main <- create_main(lhs = c("govdur", "earlyterm"))
  Ns <- create_Ns()

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = NULL,
    main = main,
    hm_blocks = NULL,
    mm = NULL,
    hm = NULL,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true("shape" %in% result$reg.table$Parameter)
})

test_that("formatJags() adds metadata attributes", {
  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "sigma" = 0.8
  ))

  main <- create_main()
  Ns <- create_Ns()

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = NULL,
    main = main,
    hm_blocks = NULL,
    mm = NULL,
    hm = NULL,
    family = "Gaussian",
    cox_intervals = NULL
  )

  reg_table <- result$reg.table

  # Check for metadata attributes
  expect_false(is.null(attr(reg_table, "estimate_type")))
  expect_false(is.null(attr(reg_table, "credible_interval")))
  expect_false(is.null(attr(reg_table, "DIC")))
  expect_false(is.null(attr(reg_table, "outcome_family")))
})

test_that("formatJags() correctly reports outcome family for Gaussian", {
  jags_out <- create_mock_jags_output(list("b[1]" = 0.5, "b[2]" = 1.2, "sigma" = 0.8))

  main <- create_main()
  Ns <- create_Ns()

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = NULL,
    main = main,
    hm_blocks = NULL,
    mm = NULL,
    hm = NULL,
    family = "Gaussian",
    cox_intervals = NULL
  )

  expect_true(grepl("Gaussian", attr(result$reg.table, "outcome_family")))
})

test_that("formatJags() reports Cox intervals correctly", {
  jags_out <- create_mock_jags_output(list("b[1]" = 0.5, "b[2]" = 1.2))

  main <- create_main(lhs = c("govdur", "earlyterm"))
  Ns <- create_Ns()

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = NULL,
    main = main,
    hm_blocks = NULL,
    mm = NULL,
    hm = NULL,
    family = "Cox",
    cox_intervals = 10
  )

  # Should mention intervals in outcome description
  expect_true(grepl("10.*interval", attr(result$reg.table, "outcome_family"), ignore.case = TRUE))
})

test_that("formatJags() handles mm() block parameters", {
  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "b.mm.1[1]" = -0.3,
    "sigma.mm" = 0.5,
    "shape" = 1.5
  ))

  main <- create_main(lhs = c("govdur", "earlyterm"))
  Ns <- create_Ns()

  mm_blocks <- list(
    list(
      vars = c("fdep"),
      vars_fixed = NULL,
      dat = data.frame(mmid = 1:20, mainid = rep(1:10, each = 2)),
      dat_fixed = NULL,
      fix_values = NULL,
      wdat = data.frame(mmid = 1:20, mainid = rep(1:10, each = 2), n = 2),
      fn = list(formula = NULL, params = character(0), vars = c("n"), vars_p = character(0)),
      RE = TRUE,
      ar = FALSE
    )
  )
  attr(mm_blocks, "has_RE") <- TRUE
  attr(mm_blocks, "has_vars") <- TRUE

  mm <- list(
    list(
      id = c("pid", "gid"),
      vars = list(formula = ~ 0 + fdep),
      fn = list(formula = NULL, params = character(0), vars = c("n")),
      RE = TRUE,
      ar = FALSE
    )
  )

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = mm_blocks,
    main = main,
    hm_blocks = NULL,
    mm = mm,
    hm = NULL,
    family = "Weibull",
    cox_intervals = NULL
  )

  reg_table <- result$reg.table

  # Should have mm-specific parameters
  expect_true(any(grepl("mm", reg_table$Parameter, ignore.case = TRUE)))
})

test_that("formatJags() handles hm() block parameters", {
  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "sigma.hm.1" = 0.3,
    "shape" = 1.5
  ))

  main <- create_main(lhs = c("govdur", "earlyterm"))
  Ns <- create_Ns()

  hm_blocks <- list(
    list(
      id = "cid",
      vars = NULL,
      vars_fixed = NULL,
      dat = data.frame(hmid = 1:5),
      dat_fixed = NULL,
      fix_values = NULL,
      name = NULL,
      type = "RE",
      showFE = FALSE,
      ar = FALSE
    )
  )

  hm <- list(
    list(
      id = "cid",
      vars = NULL,
      type = "RE",
      ar = FALSE
    )
  )

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = NULL,
    main = main,
    hm_blocks = hm_blocks,
    mm = NULL,
    hm = hm,
    family = "Weibull",
    cox_intervals = NULL
  )

  reg_table <- result$reg.table

  # Should have hm-specific parameters
  expect_true(any(grepl("hm", reg_table$Parameter, ignore.case = TRUE)))
})

test_that("formatJags() includes level specification", {
  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "b.mm.1[1]" = -0.3,
    "sigma.mm" = 0.5,
    "sigma.hm.1" = 0.3,
    "shape" = 1.5
  ))

  main <- create_main(lhs = c("govdur", "earlyterm"))
  Ns <- create_Ns()

  mm_blocks <- list(
    list(
      vars = c("fdep"),
      vars_fixed = NULL,
      dat = data.frame(mmid = 1:20, mainid = rep(1:10, each = 2)),
      dat_fixed = NULL,
      fix_values = NULL,
      wdat = data.frame(mmid = 1:20, mainid = rep(1:10, each = 2), n = 2),
      fn = list(formula = NULL, params = character(0), vars = c("n"), vars_p = character(0)),
      RE = TRUE,
      ar = FALSE
    )
  )
  attr(mm_blocks, "has_RE") <- TRUE
  attr(mm_blocks, "has_vars") <- TRUE

  mm <- list(
    list(
      id = c("pid", "gid"),
      vars = list(formula = ~ 0 + fdep),
      fn = list(formula = NULL, params = character(0), vars = c("n")),
      RE = TRUE,
      ar = FALSE
    )
  )

  hm_blocks <- list(
    list(
      id = "cid",
      vars = NULL,
      vars_fixed = NULL,
      dat = data.frame(hmid = 1:5),
      dat_fixed = NULL,
      fix_values = NULL,
      name = NULL,
      type = "RE",
      showFE = FALSE,
      ar = FALSE
    )
  )

  hm <- list(
    list(
      id = "cid",
      vars = NULL,
      type = "RE",
      ar = FALSE
    )
  )

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = mm_blocks,
    main = main,
    hm_blocks = hm_blocks,
    mm = mm,
    hm = hm,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_false(is.null(attr(result$reg.table, "level_spec")))
  expect_true(grepl("mm", attr(result$reg.table, "level_spec"), ignore.case = TRUE))
  expect_true(grepl("hm", attr(result$reg.table, "level_spec"), ignore.case = TRUE))
})

test_that("formatJags() handles fixed main-level variables", {
  jags_out <- create_mock_jags_output(list(
    "b[1]" = 1.2,
    "sigma" = 0.8
  ))

  main <- list(
    vars = c("majority"),
    vars_fixed = list(list(var = "X0", value = 0)),
    dat = data.frame(mainid = 1:10),
    dat_fixed = data.frame(X0 = rep(1, 10)),
    fix_values = c(0),
    lhs = "sim.y",
    formula = ~ 0 + majority
  )
  Ns <- create_Ns()

  result <- bml:::formatJags(
    jags.out = jags_out,
    monitor = FALSE,
    Ns = Ns,
    mm_blocks = NULL,
    main = main,
    hm_blocks = NULL,
    mm = NULL,
    hm = NULL,
    family = "Gaussian",
    cox_intervals = NULL
  )

  # Should have fixed parameter in output
  expect_true(any(grepl("fixed", result$reg.table$Parameter, ignore.case = TRUE)))
})
