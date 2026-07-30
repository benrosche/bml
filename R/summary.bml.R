#' Summarize a fitted bml model
#'
#' @description
#' S3 method for summarizing \code{bml} model objects. Returns a formatted table
#' of parameter estimates with posterior means, standard deviations, and credible
#' intervals, along with model information and convergence statistics.
#'
#' @param object A fitted model object of class \code{"bml"} returned by \code{\link{bml}}.
#' @param r Number of decimal places for rounding numeric output. Default: 3.
#' @param diagnostics Logical; include per-parameter R-hat and effective sample
#'   sizes in the printed table? Default: \code{TRUE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class \code{"bml_summary"} containing rounded parameter
#'   estimates with the following columns:
#'   \itemize{
#'     \item \code{Parameter}: Labeled parameter names
#'     \item \code{Estimate}: Posterior mean
#'     \item \code{Est.Error}: Posterior standard deviation
#'     \item \code{Q2.5}/\code{Q97.5}: 95\% credible interval
#'     \item \code{Rhat}, \code{ESS_bulk}, and \code{ESS_tail}: convergence
#'       diagnostics when \code{diagnostics = TRUE}
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
#' accessible via standard data frame indexing (e.g., \code{$Parameter},
#' \code{$Estimate}).
#'
#' For Cox models with piecewise baseline hazards (when \code{cox_intervals} is
#' specified), the outcome description includes the number of intervals used.
#'
#' @seealso \code{\link{bml}}, \code{\link{tidy.bml}}, \code{\link{coefPlot}},
#'   \code{\link{monetPlot}}, \code{\link{mcmcDiag}}
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' # Fit model
#' m1 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("sum"), RE = TRUE) +
#'     hm(id = id(cid)),
#'   family = weibull(),
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
#' s$Estimate   # Posterior means
#' s$Q2.5       # Lower credible bounds
#'
#' # Custom posterior summaries (requires monitor = "parameters")
#' draws <- posterior::as_draws_df(m1)
#'
#' # Select specific parameters and compute custom summaries
#' draws |>
#'   dplyr::select(dplyr::starts_with("b[")) |>
#'   tidyr::pivot_longer(everything(), names_to = "param") |>
#'   dplyr::group_by(param) |>
#'   dplyr::summarise(
#'     median = median(value),
#'     mad    = mad(value),
#'     q05    = quantile(value, 0.05),
#'     q95    = quantile(value, 0.95)
#'   )
#' }
#'
#' @exportS3Method summary bml
#' @author Benjamin Rosche <benrosche@@nyu.edu>

summary.bml <- function(object, r = 3, diagnostics = TRUE, ...) {

  diag <- object$diagnostics
  table <- object$reg.table |>
    tibble::rownames_to_column("node")

  if (!is.null(diag) && nrow(diag) > 0L) {
    table <- table |>
      dplyr::left_join(
        diag |>
          dplyr::select(node, rhat, ess_bulk, ess_tail),
        by = "node"
      )
  } else {
    table$rhat <- NA_real_
    table$ess_bulk <- NA_real_
    table$ess_tail <- NA_real_
  }

  rounded_table <- table |>
    dplyr::select(-node) |>
    dplyr::rename(
      Estimate = mean,
      Est.Error = sd,
      Q2.5 = lb,
      Q97.5 = ub,
      Rhat = rhat,
      ESS_bulk = ess_bulk,
      ESS_tail = ess_tail
    )

  if (!isTRUE(diagnostics)) {
    rounded_table <- rounded_table |>
      dplyr::select(-Rhat, -ESS_bulk, -ESS_tail)
  }
  rounded_table <- rounded_table |>
    dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, r)))

  # Preserve metadata
  attr(rounded_table, "estimate_type") <- attr(object$reg.table, "estimate_type")
  attr(rounded_table, "credible_interval") <- attr(object$reg.table, "credible_interval")
  attr(rounded_table, "DIC") <- attr(object$reg.table, "DIC")
  attr(rounded_table, "level_spec") <- attr(object$reg.table, "level_spec")
  attr(rounded_table, "outcome_family") <- attr(object$reg.table, "outcome_family")
  attr(rounded_table, "formula") <- object$input$formula
  attr(rounded_table, "n.main") <- object$input$n.main
  attr(rounded_table, "n.umm") <- object$input$n.umm
  attr(rounded_table, "n.mm") <- object$input$n.mm
  attr(rounded_table, "iter") <- object$input$iter
  attr(rounded_table, "warmup") <- object$input$warmup
  attr(rounded_table, "thin") <- object$input$thin
  attr(rounded_table, "chains") <- object$input$chains
  attr(rounded_table, "seed") <- object$input$seed
  attr(rounded_table, "storage") <- object$storage
  attr(rounded_table, "diagnostic_overview") <- bml_diagnostic_overview(object)

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
  formula <- attr(x, "formula")
  storage <- attr(x, "storage")
  diagnostics <- attr(x, "diagnostic_overview")

  if (!is.null(formula)) {
    cat("Formula:", paste(deparse(formula), collapse = " "), "\n")
  }
  if (!is.null(outcome_family)) {
    cat("Outcome family:", outcome_family, "\n")
  }
  data_parts <- paste(attr(x, "n.main"), "outcome units")
  if (length(attr(x, "n.umm")) == 1L &&
      !is.na(attr(x, "n.umm")) &&
      attr(x, "n.umm") > 0L) {
    data_parts <- c(data_parts, paste(attr(x, "n.umm"), "unique members"))
  }
  if (length(attr(x, "n.mm")) == 1L &&
      !is.na(attr(x, "n.mm")) &&
      attr(x, "n.mm") > 0L) {
    data_parts <- c(data_parts, paste(attr(x, "n.mm"), "membership rows"))
  }
  cat("Data:", paste(data_parts, collapse = "; "), "\n")
  cat(
    "MCMC:", attr(x, "chains"), "chains;",
    attr(x, "iter"), "iterations;",
    attr(x, "warmup"), "warmup; thin =", attr(x, "thin"), "\n"
  )
  cat(
    "Storage:", storage$label, "(", storage$ndraws, "retained draws); seed =",
    attr(x, "seed"), "\n"
  )
  if (!is.null(diagnostics)) {
    cat(
      "Convergence:", diagnostics$convergence,
      "| max R-hat =", format(round(diagnostics$max_rhat, 3), nsmall = 3),
      "| min bulk ESS =", round(diagnostics$min_ess_bulk),
      "| min tail ESS =", round(diagnostics$min_ess_tail), "\n"
    )
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

  class(x) <- setdiff(class(x), "bml_summary")

  # Print by section: coefficients / weights / fn shape params / variance
  comp <- x$component
  if (!is.null(comp)) {
    xp <- x
    xp$component <- NULL

    sections <- list(
      "Coefficients (class \"b\"):" = comp == "fixed",
      "Weight parameters (class \"w\"):" = comp == "weights",
      "Function parameters (class \"fn\"):" = comp == "fn",
      "Variance parameters:" = comp == "random",
      "Other parameters:" = !comp %in% c("fixed", "weights", "fn", "random")
    )

    for (title in names(sections)) {
      rows <- sections[[title]]
      if (any(rows)) {
        cat(title, "\n")
        print(xp[rows, , drop = FALSE], ...)
        cat("\n")
      }
    }
  } else {
    print(x, ...)
  }

  # Print DIC at the bottom
  if (!is.null(dic_value)) {
    cat("Model fit: DIC =", dic_value, "\n")
  }
}
