#' Variance decomposition for fitted bml models
#'
#' @description
#' Computes a posterior variance decomposition and intraclass correlation
#' coefficients (ICCs) from a fitted \code{bml} model. The function automatically
#' discovers all variance components (sigma parameters) in the model, applies
#' weight adjustments for multiple-membership levels, and returns posterior
#' summaries.
#'
#' @param model A fitted model object of class \code{"bml"} returned by \code{\link{bml}}.
#'   Must have been fitted with \code{monitor = TRUE} (the default).
#' @param uncertainty Uncertainty measure to report. One of \code{"sd"} (posterior
#'   standard deviation, the default), \code{"mad"} (median absolute deviation), or
#'   \code{"ci"} (95\% credible interval with lower/upper bounds).
#' @param r Number of decimal places for rounding numeric output. Default: 2.
#'
#' @return A data frame of class \code{"bml_varDecomp"} with one row per variance
#'   component. Always includes \code{Component}, \code{sigma}, and \code{ICC} columns.
#'   Additional columns depend on \code{uncertainty}:
#'   \itemize{
#'     \item \code{"sd"}: \code{sigma_sd} and \code{ICC_sd}
#'     \item \code{"mad"}: \code{sigma_mad} and \code{ICC_mad}
#'     \item \code{"ci"}: \code{sigma_lb}, \code{sigma_ub}, \code{ICC_lb}, \code{ICC_ub}
#'   }
#'
#' @details
#' \strong{Variance decomposition.}
#' The total variance of the outcome is partitioned into additive components,
#' one for each level in the model:
#'
#' \deqn{\mathrm{Var}(y) = \sigma_1^2 + \sigma_2^2 + \ldots}
#'
#' Each component contributes variance \eqn{\sigma^2}, except for multiple-membership (MM) levels,
#' where the effective variance contribution is scaled by the average of the summed squared weights
#' across groups.
#'
#' \deqn{\mathrm{Var}_{\mathrm{mm}} = \sigma_{\mathrm{mm}}^2 \cdot \overline{w^2},
#' \quad \overline{w^2} = \frac{1}{N} \sum_{i=1}^{N} \sum_{k} w_{ik}^2}
#'
#' This weight adjustment accounts for the fact that the member-level variance is
#' distributed across multiple members with potentially unequal influence. With
#' equal weights (\eqn{w_{ik} = 1/n_i}), the effective variance shrinks as group
#' size increases.
#'
#' \strong{Intraclass Correlation Coefficient (ICC).}
#' The ICC for a given level is the proportion of total variance attributable to
#' that level:
#'
#' \deqn{\rho_l = \frac{\sigma_l^2}{\sum_{l'} \sigma_{l'}^2}}
#'
#' Intuitively, the ICC answers: "What fraction of the total variation in the
#' outcome is due to differences between units at this level?" An ICC of 0.30
#' for the country level, for example, means that 30\% of the outcome variation
#' can be attributed to between-country differences.
#'
#' ICCs are computed per posterior draw and then summarized, properly propagating
#' uncertainty from the MCMC samples.
#'
#' \strong{Family-specific handling:}
#' \itemize{
#'   \item \strong{Gaussian / Weibull}: The residual \code{sigma} from the model
#'     is used directly.
#'   \item \strong{Binomial}: There is no residual sigma. The latent logistic
#'     residual variance \eqn{\pi^2/3 \approx 3.29} is used instead.
#'   \item \strong{Cox}: There is no residual variance. ICCs are computed among
#'     the non-residual components only.
#' }
#'
#' @seealso \code{\link{bml}}, \code{\link{summary.bml}}
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' m1 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE) +
#'     hm(id = id(cid), type = "RE"),
#'   family = "Weibull",
#'   data = coalgov
#' )
#'
#' varDecomp(m1)
#' varDecomp(m1, uncertainty = "ci")
#' varDecomp(m1, uncertainty = "mad")
#' }
#'
#' @export
#' @author Benjamin Rosche <benrosche@@nyu.edu>

varDecomp <- function(model, uncertainty = "sd", r = 2) {

  # Validate input
  if (!inherits(model, "bml")) {
    stop("'model' must be a fitted bml object.", call. = FALSE)
  }
  if (is.null(model$jags.out) || is.null(model$jags.out$BUGSoutput$sims.matrix)) {
    stop("Posterior samples not available. Ensure monitor = TRUE when fitting the model.", call. = FALSE)
  }
  uncertainty <- match.arg(uncertainty, c("sd", "mad", "ci"))

  sims <- model$jags.out$BUGSoutput$sims.matrix
  family <- model$input$family
  col_names <- colnames(sims)

  # Discover sigma parameters
  sigma_cols <- col_names[grepl("^sigma(\\.mm\\.\\d+|\\.hm\\.\\d+)?$", col_names)]

  if (length(sigma_cols) == 0 && family != "Binomial" && family != "Cox") {
    stop("No sigma parameters found in posterior samples.", call. = FALSE)
  }

  # Build variance components
  components <- list()

  for (sc in sigma_cols) {
    sigma_samples <- sims[, sc]

    if (grepl("^sigma\\.mm\\.", sc)) {
      # MM-level: needs weight adjustment
      g <- as.integer(sub("^sigma\\.mm\\.", "", sc))

      # Find the mm block with RE = TRUE for this mmid group
      re_block_k <- NULL
      for (k in seq_along(model$input$mm)) {
        block <- model$input$mm[[k]]
        if (block$RE && !is.null(block$mmid_group) && block$mmid_group == g) {
          re_block_k <- k
          break
        }
      }

      if (is.null(re_block_k) || is.null(model$w[[re_block_k]])) {
        warning("Could not find weight matrix for ", sc, ". Skipping.", call. = FALSE)
        next
      }

      # w_bar2: mean over groups of sum-of-squared-weights
      wmat <- model$w[[re_block_k]]
      w_bar2 <- mean(apply(wmat, 1, function(x) sum(x^2, na.rm = TRUE)))

      var_samples <- sigma_samples^2 * w_bar2
      label <- paste0("MM (sigma.mm.", g, ")")

    } else if (grepl("^sigma\\.hm\\.", sc)) {
      # HM-level: plain sigma squared
      k <- as.integer(sub("^sigma\\.hm\\.", "", sc))
      var_samples <- sigma_samples^2
      label <- paste0("HM (sigma.hm.", k, ")")

    } else if (sc == "sigma") {
      # Residual
      var_samples <- sigma_samples^2
      label <- "Residual"
    }

    components[[sc]] <- list(
      label = label,
      sigma_samples = sigma_samples,
      var_samples = var_samples
    )
  }

  # Binomial: add latent-scale residual variance (pi^2/3)
  if (family == "Binomial" && !"sigma" %in% sigma_cols) {
    n_draws <- nrow(sims)
    components[["latent_residual"]] <- list(
      label = "Residual (latent)",
      sigma_samples = rep(pi / sqrt(3), n_draws),
      var_samples = rep(pi^2 / 3, n_draws)
    )
  }

  if (length(components) == 0) {
    stop("No variance components found for decomposition.", call. = FALSE)
  }

  # Compute total variance and ICCs per draw
  var_matrix <- sapply(components, function(comp) comp$var_samples)
  total_var <- rowSums(var_matrix)
  icc_matrix <- var_matrix / total_var

  # Build results with selected uncertainty measure
  results <- data.frame(
    Component = sapply(components, function(c) c$label),
    sigma = sapply(components, function(c) round(mean(c$sigma_samples), r)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (uncertainty == "sd") {
    results$sigma_sd <- sapply(components, function(c) round(stats::sd(c$sigma_samples), r))
    results$ICC <- round(apply(icc_matrix, 2, mean), r)
    results$ICC_sd <- round(apply(icc_matrix, 2, stats::sd), r)
  } else if (uncertainty == "mad") {
    results$sigma_mad <- sapply(components, function(c) round(stats::mad(c$sigma_samples), r))
    results$ICC <- round(apply(icc_matrix, 2, mean), r)
    results$ICC_mad <- round(apply(icc_matrix, 2, stats::mad), r)
  } else if (uncertainty == "ci") {
    results$sigma_lb <- sapply(components, function(c) round(stats::quantile(c$sigma_samples, 0.025), r))
    results$sigma_ub <- sapply(components, function(c) round(stats::quantile(c$sigma_samples, 0.975), r))
    results$ICC <- round(apply(icc_matrix, 2, mean), r)
    results$ICC_lb <- round(apply(icc_matrix, 2, stats::quantile, probs = 0.025), r)
    results$ICC_ub <- round(apply(icc_matrix, 2, stats::quantile, probs = 0.975), r)
  }

  attr(results, "family") <- family
  class(results) <- c("bml_varDecomp", class(results))

  return(results)
}

#' @exportS3Method print bml_varDecomp
print.bml_varDecomp <- function(x, ...) {
  cat("Variance Decomposition\n")
  cat("Family:", attr(x, "family"), "\n\n")

  print_df <- x
  class(print_df) <- "data.frame"
  print(print_df, row.names = FALSE)
}
