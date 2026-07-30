# ================================================================================================ #
# Compact posterior storage and capability-aware post-estimation
# ================================================================================================ #

test_that("summary storage retains estimates and compact diagnostics only", {
  m <- create_test_model("summary")

  expect_null(m$draws)
  expect_null(m$jags.out)
  expect_identical(m$storage$label, "summary")
  expect_identical(m$storage$ndraws, 0L)
  expect_true(all(c(
    "node", "Parameter", "rhat", "ess_bulk", "ess_tail",
    "mcse_mean", "mcse_sd", "convergence"
  ) %in% names(m$diagnostics)))
  expect_equal(nrow(m$diagnostics), nrow(m$reg.table))

  s <- summary(m)
  expect_true(all(c("Rhat", "ESS_bulk", "ESS_tail") %in% names(s)))
  expect_false(any(c("Rhat", "ESS_bulk", "ESS_tail") %in%
                     names(summary(m, diagnostics = FALSE))))
  expect_output(print(s), "Convergence:")

  td <- tidy(m)
  expect_true(all(c("rhat", "ess_bulk", "ess_tail", "convergence") %in% names(td)))
  gl <- glance(m)
  expect_identical(gl$storage, "summary")
  expect_true(all(c("max_rhat", "min_ess_bulk", "min_ess_tail", "n_flagged") %in%
                    names(gl)))

  expect_error(as_draws(m), "No posterior draws")
  expect_error(varDecomp(m), "monitor = \"parameters\"")
  expect_error(posterior_predict(m), "monitor = \"predictive\"")
})

test_that("parameter storage keeps only reported parameter draws", {
  m <- create_test_model("parameters")
  vars <- posterior::variables(m$draws)

  expect_null(m$jags.out)
  expect_true(all(vars %in% rownames(m$reg.table)))
  expect_false(any(grepl("^(re\\.|w\\.|mu\\[|pred\\[|log_lik\\[)", vars)))
  expect_s3_class(varDecomp(m), "bml_varDecomp")
  expect_s3_class(as_draws_df(m), "draws_df")
  expect_s3_class(mcmcDiag(m, parameters = "b", lag = 10), "bml_mcmcDiag")
  expect_s3_class(suppressWarnings(monetPlot(m, parameter = "b[1]")), "patchwork")
  expect_error(ranef(m), "random_effects")
  expect_error(fitted(m), "fitted")
  expect_error(conditional_effects(m), "parameters.*fitted")
})

test_that("composable capabilities support conditional effects", {
  m <- create_test_model(c("parameters", "fitted"))
  ce <- conditional_effects(m, effects = "majority", resolution = 5)

  expect_named(ce, "majority")
  expect_equal(nrow(ce$majority), 5)
  expect_true(all(c("value", "estimate", "lb", "ub") %in% names(ce$majority)))
})

test_that("predictive storage supports posterior predictive checks", {
  skip_if_not_installed("bayesplot")
  m <- create_test_model("predictive")

  expect_equal(dim(posterior_predict(m, ndraws = 10)),
               c(10, m$input$n.main))
  expect_s3_class(pp_check(m, ndraws = 10), "ggplot")
})

test_that("full storage is canonical, serializable, and capability complete", {
  m <- create_test_model("full")

  expect_s3_class(m$draws, "draws_array")
  expect_null(m$jags.out)
  expect_setequal(m$storage$capabilities, bml:::.bml_monitor_capabilities)
  expect_true(any(grepl("^re\\.mm\\.", posterior::variables(m$draws))))
  expect_true(any(grepl("^mu\\[", posterior::variables(m$draws))))
  expect_true(any(grepl("^pred\\[", posterior::variables(m$draws))))
  expect_true(any(grepl("^log_lik\\[", posterior::variables(m$draws))))

  path <- tempfile(fileext = ".rds")
  suppressWarnings(saveRDS(m, path))
  restored <- readRDS(path)
  expect_equal(as_draws_array(restored), as_draws_array(m))
  expect_s3_class(suppressWarnings(loo(restored)), "loo")
})

test_that("saved diagnostics match posterior calculations", {
  m <- create_test_model("full")
  vars <- rownames(m$reg.table)
  draws <- posterior::subset_draws(m$draws, variable = vars)
  expected <- posterior::summarise_draws(
    draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail,
    mcse_mean = posterior::mcse_mean,
    mcse_sd = posterior::mcse_sd
  )
  observed <- m$diagnostics[match(expected$variable, m$diagnostics$node), ]

  expect_equal(observed$rhat, expected$rhat)
  expect_equal(observed$ess_bulk, expected$ess_bulk)
  expect_equal(observed$ess_tail, expected$ess_tail)
  expect_equal(observed$mcse_mean, expected$mcse_mean)
  expect_equal(observed$mcse_sd, expected$mcse_sd)
})

test_that("storage levels have the expected object-size ordering", {
  summary_fit <- create_test_model("summary")
  parameter_fit <- create_test_model("parameters")
  full_fit <- create_test_model("full")

  expect_lt(as.numeric(object.size(summary_fit)), as.numeric(object.size(parameter_fit)))
  expect_lt(as.numeric(object.size(parameter_fit)), as.numeric(object.size(full_fit)))
})

test_that("raw JAGS output is retained only when explicitly requested", {
  m <- create_test_model("jags")

  expect_null(m$draws)
  expect_false(is.null(m$jags.out))
  expect_identical(m$storage$backend, "R2jags")
  expect_s3_class(as_draws_array(m), "draws_array")
  expect_s3_class(varDecomp(m), "bml_varDecomp")
})

test_that("weibull fits support compact parameter storage", {
  skip_if_no_jags()

  m <- suppressWarnings(suppressMessages(bml(
    Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE),
    data = coalgov,
    family = weibull(),
    monitor = "parameters",
    iter = 500, warmup = 250, chains = 2, seed = 101
  )))

  expect_s3_class(m$draws, "draws_array")
  expect_null(m$jags.out)
  expect_false(any(grepl("^(pred|log_lik)\\[", posterior::variables(m$draws))))
  expect_s3_class(varDecomp(m), "bml_varDecomp")
})
