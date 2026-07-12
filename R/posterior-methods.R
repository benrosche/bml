# ================================================================================================ #
# Post-estimation methods: the brms / posterior vocabulary
# ================================================================================================ #

# ------------------------------------------------------------------------------------------------ #
# Internal draw access
# ------------------------------------------------------------------------------------------------ #

# The full posterior as [iteration, chain, variable], from R2jags
bml_sims_array <- function(object) {
  if (is.null(object$jags.out)) {
    stop("Posterior draws not stored. Fit the model with monitor = TRUE.", call. = FALSE)
  }
  object$jags.out$BUGSoutput$sims.array
}

# Draws for all indexed columns of one variable (e.g. "pred"), as
# [draws, n_index] with chains stacked
bml_draws_of <- function(object, varname) {
  arr <- bml_sims_array(object)
  vars <- dimnames(arr)[[3]]
  cols <- grep(paste0("^", varname, "(\\[|$)"), vars, value = TRUE)
  if (length(cols) == 0) {
    return(NULL)
  }
  # order by index
  idx <- as.numeric(sub(".*\\[(\\d+)\\]$", "\\1", cols))
  if (!any(is.na(idx))) cols <- cols[order(idx)]
  sub_arr <- arr[, , cols, drop = FALSE]
  matrix(sub_arr, nrow = dim(arr)[1] * dim(arr)[2], ncol = length(cols),
         dimnames = list(NULL, cols))
}

# ------------------------------------------------------------------------------------------------ #
# as_draws (posterior package)
# ------------------------------------------------------------------------------------------------ #

#' @importFrom posterior as_draws
#' @export
posterior::as_draws

#' @importFrom posterior as_draws_df
#' @export
posterior::as_draws_df

#' @importFrom posterior as_draws_matrix
#' @export
posterior::as_draws_matrix

#' @importFrom posterior as_draws_array
#' @export
posterior::as_draws_array

#' Extract posterior draws from a bml model
#'
#' @description
#' Converts the stored MCMC chains into the \pkg{posterior} package's draws
#' formats, unlocking the whole \pkg{posterior}/\pkg{bayesplot} toolchain
#' (\code{summarise_draws()}, \code{rhat()}, \code{ess_bulk()}, ...). Variables
#' keep their internal JAGS node names (\code{b[1]}, \code{b.fn.1[1]},
#' \code{sigma.mm.1}, ...); the mapping to term labels is in
#' \code{model$reg.table} (rownames = node, \code{Parameter} = label).
#'
#' @param x A fitted \code{bml} model (fitted with \code{monitor = TRUE}).
#' @param ... Passed on to the \pkg{posterior} converters.
#'
#' @return A \code{draws_array} (or the format of the called variant).
#' @name as_draws.bml
#' @export
as_draws.bml <- function(x, ...) {
  posterior::as_draws_array(bml_sims_array(x), ...)
}

#' @rdname as_draws.bml
#' @export
as_draws_df.bml <- function(x, ...) {
  posterior::as_draws_df(as_draws.bml(x), ...)
}

#' @rdname as_draws.bml
#' @export
as_draws_matrix.bml <- function(x, ...) {
  posterior::as_draws_matrix(as_draws.bml(x), ...)
}

#' @rdname as_draws.bml
#' @export
as_draws_array.bml <- function(x, ...) {
  as_draws.bml(x, ...)
}

# ------------------------------------------------------------------------------------------------ #
# Coefficient accessors: fixef / ranef / coef
# ------------------------------------------------------------------------------------------------ #

#' @importFrom nlme fixef
#' @export
nlme::fixef

#' @importFrom nlme ranef
#' @export
nlme::ranef

#' Extract fixed coefficients from a bml model
#'
#' @description
#' Returns all class-\code{"b"} coefficients — main-formula terms, block
#' features, and named-feature interactions — labeled by term, one row per
#' coefficient. \code{fn} shape parameters are deliberately excluded (they are
#' not coefficients; see the "Function parameters" section of
#' \code{\link{summary.bml}}).
#'
#' @param object A fitted \code{bml} model.
#' @param ... Unused.
#'
#' @return A matrix with columns \code{Estimate}, \code{Est.Error}, \code{Q2.5},
#'   \code{Q97.5} and term rownames.
#' @export
fixef.bml <- function(object, ...) {
  rt <- object$reg.table
  rows <- rt$component == "fixed"
  out <- as.matrix(rt[rows, c("mean", "sd", "lb", "ub")])
  colnames(out) <- c("Estimate", "Est.Error", "Q2.5", "Q97.5")
  rownames(out) <- rt$Parameter[rows]
  out
}

#' Extract random effects from a bml model
#'
#' @description
#' Returns the posterior means of the unit-level random effects: member effects
#' per \code{mm()} member-id group (intercepts, and slopes when \code{re(1 + x)}
#' was used) and unit effects per \code{hm()} block.
#'
#' @param object A fitted \code{bml} model (fitted with \code{monitor = TRUE}).
#' @param ... Unused.
#'
#' @return A list with elements \code{mm} and \code{hm}.
#' @export
ranef.bml <- function(object, ...) {
  list(mm = object$re.mm, hm = object$re.hm)
}

#' Extract combined coefficients from a bml model
#'
#' @description
#' Following brms semantics loosely: the fixed coefficients plus, per grouping
#' structure, the unit-level effects (fixed part + random deviation).
#'
#' @param object A fitted \code{bml} model.
#' @param ... Unused.
#'
#' @return A list with elements \code{fixef} (matrix) and \code{ranef} (list).
#' @export
coef.bml <- function(object, ...) {
  list(fixef = fixef.bml(object), ranef = ranef.bml(object))
}

# ------------------------------------------------------------------------------------------------ #
# fitted / predict / posterior_predict
# ------------------------------------------------------------------------------------------------ #

#' Posterior expectations of the linear predictor
#'
#' @param object A fitted \code{bml} model (fitted with \code{monitor = TRUE}).
#' @param summary Logical; if \code{TRUE} (default), return a summary matrix
#'   (one row per observation), else the draws matrix (draws x observations).
#' @param ... Unused.
#'
#' @return A matrix of posterior summaries of \code{mu}, or a draws matrix.
#' @export
fitted.bml <- function(object, summary = TRUE, ...) {
  draws <- bml_draws_of(object, "mu")
  if (is.null(draws)) {
    stop("The linear predictor was not monitored. Fit the model with monitor = TRUE.",
         call. = FALSE)
  }
  if (!summary) return(unname(draws))
  out <- cbind(
    Estimate = colMeans(draws),
    Est.Error = apply(draws, 2, stats::sd),
    Q2.5 = apply(draws, 2, stats::quantile, 0.025),
    Q97.5 = apply(draws, 2, stats::quantile, 0.975)
  )
  rownames(out) <- NULL
  out
}

#' Draws from the posterior predictive distribution
#'
#' @description
#' Returns draws of replicated outcomes \code{y_rep} (the monitored \code{pred}
#' node). \code{newdata} is not supported yet: multiple-membership prediction
#' for new units requires the new units' membership map, which is future work.
#'
#' @param object A fitted \code{bml} model (fitted with \code{monitor = TRUE}).
#' @param ndraws Optional number of draws to return (subsampled without
#'   replacement).
#' @param ... Unused.
#'
#' @return A draws x observations matrix.
#' @export
posterior_predict <- function(object, ...) {
  UseMethod("posterior_predict")
}

#' @rdname posterior_predict
#' @export
posterior_predict.bml <- function(object, ndraws = NULL, ...) {
  dots <- list(...)
  if ("newdata" %in% names(dots)) {
    stop("posterior_predict(newdata = ) is not supported yet: multiple-membership ",
         "prediction for new units requires their membership map.", call. = FALSE)
  }
  draws <- bml_draws_of(object, "pred")
  if (is.null(draws)) {
    stop("Posterior predictive draws were not monitored. Fit the model with monitor = TRUE.",
         call. = FALSE)
  }
  draws <- unname(draws)
  if (!is.null(ndraws) && ndraws < nrow(draws)) {
    draws <- draws[sample.int(nrow(draws), ndraws), , drop = FALSE]
  }
  draws
}

#' Predictions from a bml model
#'
#' @description
#' Posterior predictive summaries (one row per observation). For the raw draws
#' use \code{\link{posterior_predict.bml}}; for the linear predictor use
#' \code{\link{fitted.bml}}.
#'
#' @param object A fitted \code{bml} model (fitted with \code{monitor = TRUE}).
#' @param summary Logical; if \code{TRUE} (default), return a summary matrix,
#'   else the draws matrix.
#' @param ... Unused (\code{newdata} is not supported yet).
#'
#' @return A matrix of posterior predictive summaries, or a draws matrix.
#' @export
predict.bml <- function(object, summary = TRUE, ...) {
  draws <- posterior_predict.bml(object, ...)
  if (!summary) return(draws)
  out <- cbind(
    Estimate = colMeans(draws),
    Est.Error = apply(draws, 2, stats::sd),
    Q2.5 = apply(draws, 2, stats::quantile, 0.025),
    Q97.5 = apply(draws, 2, stats::quantile, 0.975)
  )
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------------------------ #
# log_lik / loo / waic
# ------------------------------------------------------------------------------------------------ #

#' Pointwise log-likelihood of a bml model
#'
#' @param object A fitted \code{bml} model. Available for gaussian and
#'   bernoulli families fitted with \code{monitor = TRUE} (the JAGS model
#'   monitors a pointwise \code{log_lik} node).
#' @param ... Unused.
#'
#' @return A draws x observations matrix of pointwise log-likelihood values.
#' @export
log_lik <- function(object, ...) {
  UseMethod("log_lik")
}

#' @rdname log_lik
#' @export
log_lik.bml <- function(object, ...) {
  draws <- bml_draws_of(object, "log_lik")
  if (is.null(draws)) {
    stop("Pointwise log-likelihood not available. It is monitored for gaussian and ",
         "bernoulli families when monitor = TRUE (survival families are not supported yet).",
         call. = FALSE)
  }
  unname(draws)
}

# log_lik as [iteration, chain, observation] for relative_eff
bml_loglik_array <- function(object) {
  arr <- bml_sims_array(object)
  vars <- dimnames(arr)[[3]]
  cols <- grep("^log_lik\\[", vars, value = TRUE)
  if (length(cols) == 0) return(NULL)
  idx <- as.numeric(sub(".*\\[(\\d+)\\]$", "\\1", cols))
  cols <- cols[order(idx)]
  arr[, , cols, drop = FALSE]
}

#' @importFrom loo loo
#' @export
loo::loo

#' @importFrom loo waic
#' @export
loo::waic

#' Efficient approximate leave-one-out cross-validation for bml models
#'
#' @description
#' PSIS-LOO via the \pkg{loo} package, computed from the monitored pointwise
#' log-likelihood. Available for gaussian and bernoulli families fitted with
#' \code{monitor = TRUE}. Compare models with \code{loo::loo_compare()}.
#'
#' @param x A fitted \code{bml} model.
#' @param ... Passed to \code{loo::loo.matrix()}.
#'
#' @return A \code{loo} object.
#' @export
loo.bml <- function(x, ...) {
  ll_arr <- bml_loglik_array(x)
  if (is.null(ll_arr)) {
    stop("Pointwise log-likelihood not available. loo() is supported for gaussian and ",
         "bernoulli families fitted with monitor = TRUE.", call. = FALSE)
  }
  ll_mat <- matrix(ll_arr, nrow = dim(ll_arr)[1] * dim(ll_arr)[2], ncol = dim(ll_arr)[3])
  r_eff <- loo::relative_eff(exp(ll_arr))
  loo::loo(ll_mat, r_eff = r_eff, ...)
}

#' Widely applicable information criterion for bml models
#'
#' @param x A fitted \code{bml} model (see \code{\link{loo.bml}} for availability).
#' @param ... Passed to \code{loo::waic.matrix()}.
#'
#' @return A \code{waic} object.
#' @export
waic.bml <- function(x, ...) {
  ll_arr <- bml_loglik_array(x)
  if (is.null(ll_arr)) {
    stop("Pointwise log-likelihood not available. waic() is supported for gaussian and ",
         "bernoulli families fitted with monitor = TRUE.", call. = FALSE)
  }
  ll_mat <- matrix(ll_arr, nrow = dim(ll_arr)[1] * dim(ll_arr)[2], ncol = dim(ll_arr)[3])
  loo::waic(ll_mat, ...)
}

# ------------------------------------------------------------------------------------------------ #
# pp_check
# ------------------------------------------------------------------------------------------------ #

#' Posterior predictive checks for bml models
#'
#' @description
#' Graphical posterior predictive checks via \pkg{bayesplot} (must be
#' installed). Compares the observed outcome with replicated outcomes from
#' \code{\link{posterior_predict.bml}}.
#'
#' @param object A fitted \code{bml} model (gaussian or bernoulli, fitted with
#'   \code{monitor = TRUE}).
#' @param type A \pkg{bayesplot} PPC type (without the \code{"ppc_"} prefix),
#'   e.g. \code{"dens_overlay"} (default for gaussian), \code{"bars"} (default
#'   for bernoulli), \code{"hist"}, \code{"stat"}.
#' @param ndraws Number of replicated datasets to plot. Default: 30.
#' @param ... Passed to the bayesplot function.
#'
#' @return A ggplot object.
#' @export
pp_check <- function(object, ...) {
  UseMethod("pp_check")
}

#' @rdname pp_check
#' @export
pp_check.bml <- function(object, type = NULL, ndraws = 30, ...) {
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop("pp_check() requires the 'bayesplot' package: install.packages(\"bayesplot\")",
         call. = FALSE)
  }
  y <- object$input$y
  if (is.null(y)) {
    stop("pp_check() is available for gaussian and bernoulli families.", call. = FALSE)
  }
  yrep <- posterior_predict.bml(object, ndraws = ndraws)
  if (is.null(type)) {
    type <- if (object$input$family == "Binomial") "bars" else "dens_overlay"
  }
  ppc_fun <- get(paste0("ppc_", type), envir = asNamespace("bayesplot"))
  ppc_fun(y, yrep, ...)
}

# ------------------------------------------------------------------------------------------------ #
# conditional_effects
# ------------------------------------------------------------------------------------------------ #

#' Conditional effects of main-formula predictors
#'
#' @description
#' Marginal-effect curves for main-formula covariates: the expected outcome on
#' the linear-predictor scale across a grid of the chosen predictor, holding
#' the other main-formula covariates at their means and the block features,
#' random effects, and interactions at their sample-average contribution.
#'
#' @param x A fitted \code{bml} model (fitted with \code{monitor = TRUE}).
#' @param effects Character vector of main-formula term names. Default: all
#'   numeric main-formula terms.
#' @param resolution Grid size. Default: 100.
#' @param ... Unused.
#'
#' @return A named list of data frames (one per effect) with columns
#'   \code{value}, \code{estimate}, \code{lb}, \code{ub}, suitable for plotting.
#' @export
conditional_effects <- function(x, ...) {
  UseMethod("conditional_effects")
}

#' @rdname conditional_effects
#' @export
conditional_effects.bml <- function(x, effects = NULL, resolution = 100, ...) {

  mainvars <- x$input$mainvars
  terms_avail <- setdiff(mainvars, "X0")
  if (length(terms_avail) == 0) {
    stop("No main-formula covariates to compute conditional effects for.", call. = FALSE)
  }
  if (is.null(effects)) effects <- terms_avail
  bad <- setdiff(effects, terms_avail)
  if (length(bad) > 0) {
    stop("Unknown effect(s): ", paste(bad, collapse = ", "),
         ". Available: ", paste(terms_avail, collapse = ", "), call. = FALSE)
  }

  # b draws in mainvars order
  b_draws <- bml_draws_of(x, "b")
  if (is.null(b_draws)) {
    stop("Posterior draws not stored. Fit the model with monitor = TRUE.", call. = FALSE)
  }
  mu_hat <- x$fitted
  if (length(mu_hat) == 0) mu_hat <- NULL

  # Everything in mu outside X.main %*% b (block features, REs, interactions) is
  # absorbed as a constant equal to mean(fitted) - mean(Xbar %*% b); covariates
  # other than the focal one are held at their means.
  out <- list()

  for (eff in effects) {
    p <- length(mainvars)

    if (!is.null(x$input$X.main.ranges) && eff %in% rownames(x$input$X.main.ranges)) {
      rng <- x$input$X.main.ranges[eff, ]
      grid <- seq(rng[1], rng[2], length.out = resolution)
    } else {
      stop("conditional_effects() needs the predictor ranges; refit the model with this ",
           "version of bml (ranges are stored at fit time).", call. = FALSE)
    }

    Xg <- matrix(0, nrow = resolution, ncol = p)
    colnames(Xg) <- mainvars
    if ("X0" %in% mainvars) Xg[, "X0"] <- 1
    for (v in setdiff(terms_avail, eff)) {
      Xg[, v] <- x$input$X.main.means[[v]] %||% 0
    }
    Xg[, eff] <- grid

    mu_draws <- Xg %*% t(b_draws)  # resolution x draws

    # constant offset: average contribution of everything outside X.main %*% b
    offset <- 0
    if (!is.null(mu_hat) && !is.null(x$input$X.main.means)) {
      Xbar <- matrix(unlist(x$input$X.main.means[mainvars]), nrow = 1)
      offset <- mean(mu_hat) - mean(Xbar %*% t(b_draws))
    }

    out[[eff]] <- data.frame(
      value = grid,
      estimate = rowMeans(mu_draws) + offset,
      lb = apply(mu_draws, 1, stats::quantile, 0.025) + offset,
      ub = apply(mu_draws, 1, stats::quantile, 0.975) + offset
    )
  }

  out
}

# ------------------------------------------------------------------------------------------------ #
# make_jagscode / make_jagsdata
# ------------------------------------------------------------------------------------------------ #

#' Generate the JAGS model code / data without fitting
#'
#' @description
#' The bml analogues of \code{brms::make_stancode()} / \code{make_standata()}:
#' build the model exactly as \code{\link{bml}} would, but return the generated
#' JAGS model string (\code{make_jagscode}) or the JAGS data list
#' (\code{make_jagsdata}) instead of fitting.
#'
#' @param formula,data,family,prior As in \code{\link{bml}}.
#' @param monitor As in \code{\link{bml}} (affects which nodes the model defines).
#'
#' @return \code{make_jagscode()}: a character scalar of class
#'   \code{"bml_jagscode"} (prints nicely). \code{make_jagsdata()}: a named list.
#'
#' @examples
#' \dontrun{
#' code <- make_jagscode(
#'   y ~ x + mm(id = id(pid, gid), vars = vars(z), w = w(~ 1/n)),
#'   data = dat, family = gaussian()
#' )
#' code
#' }
#'
#' @export
make_jagscode <- function(formula, data, family = stats::gaussian(), prior = NULL,
                          monitor = TRUE) {
  parts <- suppressMessages(
    bml(formula, data = data, family = family, prior = prior, run = FALSE, monitor = monitor)
  )
  structure(parts$modelstring, class = "bml_jagscode")
}

#' @rdname make_jagscode
#' @export
make_jagsdata <- function(formula, data, family = stats::gaussian(), prior = NULL,
                          monitor = TRUE) {
  parts <- suppressMessages(
    bml(formula, data = data, family = family, prior = prior, run = FALSE, monitor = monitor)
  )
  parts$jags.data
}

#' @export
print.bml_jagscode <- function(x, ...) {
  cat(x, "\n")
  invisible(x)
}
