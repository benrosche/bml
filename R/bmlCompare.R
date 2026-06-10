# Internal: validate a list of bml models and name unnamed ones after the
# expressions they were passed as (used by bmlCompare and coefPlot)
nameModels <- function(models, exprs) {
  if (length(models) == 0) {
    stop("Please provide at least one fitted bml model.", call. = FALSE)
  }
  if (!all(sapply(models, inherits, "bml"))) {
    stop("All models must be fitted bml objects.", call. = FALSE)
  }
  nms <- names(models)
  if (is.null(nms)) nms <- rep("", length(models))
  nms[nms == ""] <- exprs[nms == ""]
  names(models) <- nms
  models
}

#' Compare fitted bml models in a tidy table
#'
#' @description
#' Stacks the tidy coefficient tables of one or more fitted \code{\link{bml}}
#' models into a single long table with one row per model and term, including
#' model-level information (N and DIC) for easy comparison.
#'
#' For publication-style side-by-side regression tables, the \pkg{broom}
#' methods \code{\link{tidy.bml}} and \code{\link{glance.bml}} also make
#' \code{modelsummary::modelsummary(list(m1, m2))} work directly on bml models.
#'
#' @param ... Fitted \code{bml} model objects, optionally named
#'   (e.g., \code{bmlCompare(base = m1, weighted = m2)}). Unnamed models are
#'   labeled with the expressions they were passed as (e.g., \code{"m1"}).
#' @param conf.level Width of the equal-tailed credible interval. Default: 0.95.
#'   Other levels require the models to be fitted with \code{monitor = TRUE}
#'   (see \code{\link{tidy.bml}}).
#' @param component Which parameters to include; passed to \code{\link{tidy.bml}}.
#'   Default: \code{"fixed"} (regression coefficients).
#'
#' @return A \code{tibble} with columns \code{model}, \code{term},
#'   \code{estimate}, \code{conf.low}, \code{conf.high}, \code{N} (number of
#'   main-level units), and \code{DIC}.
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' m1 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE),
#'   family = "Weibull",
#'   data = coalgov
#' )
#' m2 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ pseat), RE = TRUE),
#'   family = "Weibull",
#'   data = coalgov
#' )
#'
#' bmlCompare(m1, m2)
#' bmlCompare(equal = m1, seatshare = m2)
#' }
#'
#' @seealso \code{\link{bml}}, \code{\link{tidy.bml}}, \code{\link{glance.bml}},
#'   \code{\link{coefPlot}}
#'
#' @author Benjamin Rosche \email{benrosche@@nyu.edu}
#' @export bmlCompare
bmlCompare <- function(..., conf.level = 0.95, component = "fixed") {

  models <- nameModels(list(...), sapply(substitute(list(...))[-1], deparse))

  purrr::imap_dfr(models, function(m, name) {
    g <- glance(m)
    tidy(m, conf.level = conf.level, component = component) %>%
      dplyr::mutate(model = name, N = g$nobs, DIC = g$DIC) %>%
      dplyr::select(model, term, estimate, conf.low, conf.high, N, DIC)
  })
}
