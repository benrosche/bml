# ================================================================================================ #
# Family functions (brms-style) and family validation
# ================================================================================================ #

#' Model families for bml
#'
#' @description
#' Family functions specify the outcome distribution and link function of a
#' \code{\link{bml}} model, following the \pkg{brms}/\code{\link[stats]{glm}}
#' convention of passing a family \emph{function} rather than a string:
#'
#' \itemize{
#'   \item \code{gaussian()}: normal outcome, identity link
#'   \item \code{bernoulli()} (alias \code{binomial()}): binary outcome, logit link
#'   \item \code{weibull()}: Weibull survival model; outcome \code{Surv(time, event)}
#'   \item \code{cox(intervals = NULL)}: Cox proportional hazards model; outcome
#'     \code{Surv(time, event)}. \code{intervals} controls the baseline hazard:
#'     \code{NULL} uses all unique event times (non-parametric); an integer k
#'     uses a piecewise-constant baseline hazard with k intervals (faster).
#' }
#'
#' Family \emph{strings} are accepted as first-class input too (as in brms):
#' \code{family = "gaussian"}, \code{"bernoulli"}, \code{"weibull"}, \code{"cox"}.
#'
#' @param link Link function. Only the canonical links are supported by the JAGS
#'   backend (\code{"identity"} for gaussian, \code{"logit"} for bernoulli,
#'   \code{"log"} for weibull/cox); passing anything else errors.
#' @param intervals For \code{cox()}: \code{NULL} for a non-parametric baseline
#'   hazard (all unique event times), or a positive integer number of intervals
#'   for a piecewise-constant baseline hazard.
#'
#' @return A \code{bml_family} object (a list with elements \code{family} and
#'   \code{link}, plus \code{intervals} for \code{cox()}).
#'
#' @examples
#' \dontrun{
#' bml(y ~ x + mm(...), data = dat, family = gaussian())
#' bml(y01 ~ x + mm(...), data = dat, family = bernoulli())
#' bml(Surv(t, ev) ~ x + mm(...), data = dat, family = weibull())
#' bml(Surv(t, ev) ~ x + mm(...), data = dat, family = cox(intervals = 10))
#' }
#'
#' @name bml-families
NULL

bml_family <- function(family, link, intervals = NULL) {
  structure(
    list(family = family, link = link, intervals = intervals),
    class = "bml_family"
  )
}

#' @rdname bml-families
#' @export
bernoulli <- function(link = "logit") {
  if (!identical(link, "logit")) {
    stop("bernoulli() only supports link = \"logit\" with the JAGS backend.", call. = FALSE)
  }
  bml_family("Binomial", "logit")
}

#' @rdname bml-families
#' @export
weibull <- function(link = "log") {
  if (!identical(link, "log")) {
    stop("weibull() only supports link = \"log\" with the JAGS backend.", call. = FALSE)
  }
  bml_family("Weibull", "log")
}

#' @rdname bml-families
#' @export
cox <- function(intervals = NULL, link = "log") {
  if (!identical(link, "log")) {
    stop("cox() only supports link = \"log\" with the JAGS backend.", call. = FALSE)
  }
  if (!is.null(intervals)) {
    if (!is.numeric(intervals) || length(intervals) != 1 || intervals < 1) {
      stop("cox(intervals = ) must be NULL or a single positive integer.", call. = FALSE)
    }
    intervals <- as.integer(intervals)
  }
  bml_family("Cox", "log", intervals = intervals)
}

# Internal: normalize whatever the user passed as `family` into a bml_family.
# Accepts: bml_family objects, base-R family objects (gaussian(), binomial()),
# and strings ("gaussian", "bernoulli", "binomial", "weibull", "cox";
# capitalized legacy spellings also resolve).
validate_family <- function(family) {

  if (inherits(family, "bml_family")) {
    return(family)
  }

  if (inherits(family, "family")) {
    fam <- tolower(family$family)
    link <- family$link
    if (fam == "gaussian") {
      if (!identical(link, "identity")) {
        stop("gaussian() only supports link = \"identity\" with the JAGS backend.", call. = FALSE)
      }
      return(bml_family("Gaussian", "identity"))
    }
    if (fam == "binomial") {
      if (!identical(link, "logit")) {
        stop("binomial() only supports link = \"logit\" with the JAGS backend.", call. = FALSE)
      }
      return(bml_family("Binomial", "logit"))
    }
    stop("Unsupported family object: '", family$family,
         "'. Supported: gaussian(), bernoulli()/binomial(), weibull(), cox().", call. = FALSE)
  }

  if (is.character(family) && length(family) == 1) {
    fam <- tolower(family)
    out <- switch(
      fam,
      "gaussian"  = bml_family("Gaussian", "identity"),
      "bernoulli" = bml_family("Binomial", "logit"),
      "binomial"  = bml_family("Binomial", "logit"),
      "weibull"   = bml_family("Weibull", "log"),
      "cox"       = bml_family("Cox", "log"),
      NULL
    )
    if (!is.null(out)) return(out)
    stop("Invalid family: '", family, "'. Must be one of: ",
         "gaussian, bernoulli, weibull, cox (function or string).", call. = FALSE)
  }

  stop("'family' must be a family function (e.g. gaussian(), bernoulli(), weibull(), cox()) ",
       "or a family string.", call. = FALSE)
}
