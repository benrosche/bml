# ================================================================================================ #
# Shared helpers for tests that need a fitted model
# ================================================================================================ #

data("coalgov")

skip_if_no_jags <- function() {
  testthat::skip_on_cran()
  if (!requireNamespace("rjags", quietly = TRUE)) {
    testthat::skip("rjags not available")
  }
}

# Minimal fitted gaussian model with an mm block (cached across tests within a run)
.test_model_cache <- new.env(parent = emptyenv())

create_test_model <- function(monitor = TRUE) {
  skip_if_no_jags()

  key <- paste0("m_", monitor)
  if (!is.null(.test_model_cache[[key]])) return(.test_model_cache[[key]])

  m <- tryCatch({
    suppressWarnings(suppressMessages(bml(
      dur_wkb ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE),
      data = coalgov,
      family = gaussian(),
      iter = 800, warmup = 300, thin = 1, chains = 2,
      seed = 99,
      monitor = monitor
    )))
  }, error = function(e) {
    testthat::skip(paste("Could not fit test model:", e$message))
  })

  .test_model_cache[[key]] <- m
  m
}
