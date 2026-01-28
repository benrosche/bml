#' Summarize a fitted bml model
#'
#' @description
#' S3 method for summarizing \code{bml} model objects. Returns a formatted table
#' of parameter estimates with posterior means, standard deviations, and credible
#' intervals, along with model information and convergence statistics.
#'
#' @param object A fitted model object of class \code{"bml"} returned by \code{\link{bml}}.
#' @param r Number of decimal places for rounding numeric output. Default: 3.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class \code{"bml_summary"} containing rounded parameter
#'   estimates with the following columns:
#'   \itemize{
#'     \item \code{Parameter}: Labeled parameter names
#'     \item \code{mean}: Posterior mean
#'     \item \code{sd}: Posterior standard deviation
#'     \item \code{lb}: Lower bound of 95\% credible interval
#'     \item \code{ub}: Upper bound of 95\% credible interval
#'   }
#'
#'   The object includes metadata attributes printed above the table:
#'   \itemize{
#'     \item Outcome family and link function
#'     \item Estimate type (posterior mean from MCMC)
#'     \item Credible interval specification (95\% equal-tailed)
#'     \item Level specification (mm and hm block details)
#'     \item DIC (Deviance Information Criterion) for model comparison
#'   }
#'
#' @details
#' The summary method rounds all numeric values for readability while preserving
#' the underlying structure and metadata from the fitted model. All columns remain
#' accessible via standard data frame indexing (e.g., \code{$Parameter}, \code{$mean}).
#'
#' For Cox models with piecewise baseline hazards (when \code{cox_intervals} is
#' specified), the outcome description includes the number of intervals used.
#'
#' @seealso \code{\link{bml}}, \code{\link{monetPlot}}, \code{\link{mcmcDiag}}
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' # Fit model
#' m1 <- bml(
#'   Surv(govdur, earlyterm) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE) +
#'     hm(id = id(cid), type = "RE"),
#'   family = "Weibull",
#'   data = coalgov
#' )
#'
#' # View summary
#' summary(m1)
#'
#' # Summary with more decimal places
#' summary(m1, r = 4)
#'
#' # Access specific columns
#' s <- summary(m1)
#' s$Parameter  # Parameter names
#' s$mean       # Posterior means
#' s$lb         # Lower credible bounds
#' }
#'
#' @exportS3Method summary bml
#' @author Benjamin Rosche <benrosche@@nyu.edu>

summary.bml <- function(object, r=3, ...) {

  rounded_table <- object$reg.table %>% dplyr::mutate(across(where(is.numeric), ~round(.,r)))

  # Preserve metadata
  attr(rounded_table, "estimate_type") <- attr(object$reg.table, "estimate_type")
  attr(rounded_table, "credible_interval") <- attr(object$reg.table, "credible_interval")
  attr(rounded_table, "DIC") <- attr(object$reg.table, "DIC")
  attr(rounded_table, "level_spec") <- attr(object$reg.table, "level_spec")
  attr(rounded_table, "outcome_family") <- attr(object$reg.table, "outcome_family")

  # Add class for custom printing
  class(rounded_table) <- c("bml_summary", class(rounded_table))

  return(rounded_table)

}

#' @exportS3Method print bml_summary
print.bml_summary <- function(x, ...) {
  # Print header with metadata
  estimate_type <- attr(x, "estimate_type")
  ci_info <- attr(x, "credible_interval")
  dic_value <- attr(x, "DIC")
  level_spec <- attr(x, "level_spec")
  outcome_family <- attr(x, "outcome_family")

  if (!is.null(outcome_family)) {
    cat("Outcome family:", outcome_family, "\n")
  }
  if (!is.null(estimate_type)) {
    cat("Estimates:", estimate_type, "\n")
  }
  if (!is.null(ci_info)) {
    cat("Uncertainty:", ci_info, "\n")
  }
  if (!is.null(level_spec)) {
    cat("\nLevel specification:\n")
    cat(level_spec, "\n")
  }
  if (!is.null(outcome_family) || !is.null(estimate_type) || !is.null(ci_info) || !is.null(level_spec)) {
    cat("\n")
  }

  # Remove custom class and print as data frame
  class(x) <- setdiff(class(x), "bml_summary")
  print(x, ...)

  # Print DIC at the bottom
  if (!is.null(dic_value)) {
    cat("\nModel fit: DIC =", dic_value, "\n")
  }
}




