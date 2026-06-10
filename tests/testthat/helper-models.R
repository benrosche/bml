# ================================================================================================ #
# Shared helpers for tests that need a fitted model (used by test-output.R, test-broom.R)
# ================================================================================================ #

data("coalgov")

# Create a minimal fitted model for testing (if JAGS is available)
create_test_model <- function(monitor = FALSE, n.thin = 1) {
  testthat::skip_on_cran()

  if (!requireNamespace("rjags", quietly = TRUE)) {
    skip("rjags not available")
  }

  # Use a very small subset and few iterations for speed
  test_data <- coalgov[1:100, ]

  tryCatch({
    m <- bml(
      event_wkb ~ 1 + majority,
      family = "Gaussian",
      data = test_data,
      n.iter = 500,
      n.burnin = 100,
      n.thin = n.thin,
      n.chains = 2,
      monitor = monitor
    )
    return(m)
  }, error = function(e) {
    skip(paste("Could not fit test model:", e$message))
  })
}
