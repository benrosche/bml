# ================================================================================================ #
# re() / fe(): the unified effects grammar for mm() and hm()
# ================================================================================================ #

# Internal: parse the lme4-style LHS expression (1, 1 + x, 0 + x, x - 1) into
# list(intercept = TRUE/FALSE, slopes = character()).
parse_effects_formula <- function(expr) {

  f <- stats::as.formula(call("~", expr))
  tt <- stats::terms(f)
  intercept <- attr(tt, "intercept") == 1
  slopes <- attr(tt, "term.labels")

  list(intercept = intercept, slopes = slopes)
}

#' Random effects: partial pooling
#'
#' @description
#' Specifies the random-effect structure of an \code{\link{mm}} or
#' \code{\link{hm}} block, mirroring the left-hand side of an lme4/brms
#' \code{( ... | group )} term. The grouping is implied by the block's
#' \code{id()} and therefore not repeated.
#'
#' \preformatted{
#' RE = TRUE            # shorthand for re(1)
#' re(1)                # random intercept
#' re(1 + x)            # intercept + slope on x (INDEPENDENT, the default)
#' re(1 + x, cor = TRUE)# intercept + slope, correlated (lme4/brms `|`) - opt-in
#' re(0 + x)            # slope only (also re(x - 1))
#' re(1, ar = TRUE)     # random walk over participation order
#' re(1, ar = year)     # random walk over calendar time (gap-scaled variance)
#' }
#'
#' The default \code{cor = FALSE} deliberately diverges from brms's correlated
#' default: the correlated draw is a real JAGS cost, and in \code{mm()} blocks
#' the member random effects are only read through the weighted aggregate
#' \code{sum_k w_k u_k}, so the correlation is rarely of interest and often
#' weakly identified. Opt in with \code{cor = TRUE} (supported for a single
#' slope, i.e. the 2x2 case).
#'
#' @param x The effects expression: \code{1}, \code{1 + var}, \code{0 + var},
#'   or \code{var - 1}. For \code{mm()} blocks, slope variables are member-level
#'   columns; for \code{hm()} blocks they are main-level columns varying within
#'   the nesting units.
#' @param cor Logical; if \code{TRUE}, intercept and slope get a joint
#'   bivariate-normal prior (correlated). Default \code{FALSE} (independent).
#' @param ar Autoregressive random walk for the intercept across a unit's
#'   repeated participations/observations. \code{FALSE} (default): independent
#'   effects. \code{TRUE}: random walk over participation order,
#'   \eqn{u_t \sim N(u_{t-1}, \sigma^2)}. An unquoted column name (e.g.
#'   \code{ar = year}): random walk over that (numeric) time variable, with the
#'   step variance scaled by the normalized time gaps,
#'   \eqn{u_t \sim N(u_{t-1}, \sigma^2 \, \Delta_t / \bar\Delta)} — equal
#'   spacing reproduces the \code{ar = TRUE} model, and \eqn{\sigma} keeps its
#'   per-average-step interpretation. Requires an intercept-only \code{re(1)};
#'   time values must be unique within each unit.
#'
#' @return A \code{bml_re} object.
#' @seealso \code{\link{fe}}, \code{\link{mm}}, \code{\link{hm}}
#' @export
re <- function(x = 1, cor = FALSE, ar = FALSE) {
  expr <- substitute(x)
  parsed <- parse_effects_formula(expr)

  # ar: FALSE -> NULL; TRUE -> order walk; unquoted column/string -> time-indexed walk.
  # Work on the unevaluated expression first: a bare column name (ar = year) must not
  # be evaluated as a variable.
  ar_expr <- substitute(ar)
  ar_spec <- if (is.logical(ar_expr)) {
    # literal TRUE/FALSE (including the FALSE default)
    if (isTRUE(ar_expr)) list(time = NULL) else NULL
  } else {
    val <- tryCatch(eval(ar_expr, parent.frame()), error = function(e) e)
    if (isTRUE(val)) {
      list(time = NULL)
    } else if (isFALSE(val) || is.null(val)) {
      NULL
    } else if (is.character(val) && length(val) == 1) {
      list(time = val)
    } else {
      # unresolvable or non-scalar: treat as a literal column name
      list(time = deparse(ar_expr))
    }
  }

  if (!parsed$intercept && length(parsed$slopes) == 0) {
    stop("re() must include an intercept and/or at least one slope.", call. = FALSE)
  }
  if (isTRUE(cor) && length(parsed$slopes) == 0) {
    stop("re(cor = TRUE) requires at least one slope (there is nothing to correlate).",
         call. = FALSE)
  }
  if (isTRUE(cor) && (length(parsed$slopes) > 1 || !parsed$intercept)) {
    stop("re(cor = TRUE) is currently supported for exactly one slope plus the ",
         "intercept (the 2x2 case). Use cor = FALSE for larger structures.", call. = FALSE)
  }
  if (!is.null(ar_spec)) {
    if (length(parsed$slopes) > 0 || !parsed$intercept) {
      stop("re(ar = ...) requires an intercept-only random effect: re(1, ar = ...).",
           call. = FALSE)
    }
    if (isTRUE(cor)) {
      stop("re(ar = ...) cannot be combined with cor = TRUE.", call. = FALSE)
    }
  }

  structure(
    list(intercept = parsed$intercept, slopes = parsed$slopes, cor = isTRUE(cor),
         ar = ar_spec),
    class = "bml_re"
  )
}

#' Fixed effects: no pooling
#'
#' @description
#' The no-pooling peer of \code{\link{re}}: identical formula grammar, but unit
#' effects get flat independent priors instead of a hierarchical prior
#' (Gelman-Hill: "FE = RE with infinite variance"). \code{fe(1)} gives classic
#' unit dummies; \code{fe(1 + x)} additionally gives unit-specific slopes.
#' The first unit is the reference category (its effects are fixed to 0), so
#' \code{fe()} remains identified alongside a global intercept and global
#' slopes.
#'
#' @param x The effects expression: \code{1}, \code{1 + var}, \code{0 + var}.
#' @param showFE Logical; if \code{TRUE}, per-unit estimates are monitored and
#'   reported. Default \code{FALSE}.
#'
#' @return A \code{bml_fe} object.
#' @seealso \code{\link{re}}, \code{\link{hm}}, \code{\link{mm}}
#' @export
fe <- function(x = 1, showFE = FALSE) {
  expr <- substitute(x)
  parsed <- parse_effects_formula(expr)

  if (!parsed$intercept && length(parsed$slopes) == 0) {
    stop("fe() must include an intercept and/or at least one slope.", call. = FALSE)
  }

  structure(
    list(intercept = parsed$intercept, slopes = parsed$slopes, showFE = isTRUE(showFE)),
    class = "bml_fe"
  )
}
