# ================================================================================================ #
# Broom-style methods: tidy() and glance() for bml objects
# ================================================================================================ #

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' Tidy a fitted bml model
#'
#' @description
#' Turns a fitted \code{\link{bml}} model into a tidy tibble following the
#' \pkg{broom} conventions: one row per parameter with columns \code{term},
#' \code{estimate}, \code{std.error}, \code{conf.low}, and \code{conf.high}.
#' Because the model is Bayesian, \code{estimate} is the posterior mean,
#' \code{std.error} the posterior standard deviation, and
#' \code{conf.low}/\code{conf.high} the bounds of an equal-tailed credible
#' interval.
#'
#' With these methods, ecosystem tools that build on \pkg{broom} work out of
#' the box, e.g. \code{modelsummary::modelsummary(list(m1, m2))} for
#' side-by-side regression tables or \code{dotwhisker::dwplot(m1)} for
#' coefficient plots.
#'
#' @param x A fitted model object of class \code{"bml"} returned by \code{\link{bml}}.
#' @param conf.int Logical. Include credible interval columns
#'   \code{conf.low}/\code{conf.high}? Default: \code{TRUE}.
#' @param conf.level Width of the equal-tailed credible interval. Default: 0.95,
#'   which is read directly from the fitted model. Other levels are computed
#'   from the posterior draws and therefore require the model to be fitted with
#'   \code{monitor = TRUE}.
#' @param component Which parameters to return:
#'   \itemize{
#'     \item \code{"all"} (default): all of the below
#'     \item \code{"fixed"}: regression coefficients (main, mm, and hm level)
#'     \item \code{"random"}: variance/scale parameters (e.g., \code{sigma},
#'       \code{sigma.mm.*}, \code{sigma.hm.*}, Weibull \code{shape}, \code{tau})
#'     \item \code{"weights"}: weight function parameters
#'   }
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{tibble} with columns \code{term}, \code{estimate},
#'   \code{std.error}, \code{conf.low}, \code{conf.high} (if \code{conf.int}),
#'   and \code{component}.
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
#'
#' tidy(m1)                      # regression coefficients
#' tidy(m1, component = "all")   # including weight and variance parameters
#' tidy(m1, conf.level = 0.9)    # 90% credible intervals (requires monitor = TRUE)
#'
#' # Multi-model comparison table (see also bmlCompare()):
#' models <- list(base = m1, weighted = m2)
#' purrr::imap_dfr(models, \(m, name)
#'   tidy(m) |>
#'     dplyr::mutate(model = name, N = glance(m)$nobs, DIC = glance(m)$DIC)
#' )
#' }
#'
#' @seealso \code{\link{bml}}, \code{\link{glance.bml}}, \code{\link{bmlCompare}},
#'   \code{\link{coefPlot}}, \code{\link{summary.bml}}
#'
#' @author Benjamin Rosche \email{benrosche@@nyu.edu}
#' @export
tidy.bml <- function(x, conf.int = TRUE, conf.level = 0.95,
                     component = c("all", "fixed", "random", "weights"), ...) {

  component <- match.arg(component)

  if (is.null(x$reg.table)) {
    stop("No estimates found. Please fit the model with run = TRUE.", call. = FALSE)
  }

  rt <- x$reg.table
  jagsnames <- rownames(rt)

  comp <- dplyr::case_when(
    stringr::str_detect(jagsnames, "^b\\.w\\.") ~ "weights",
    stringr::str_detect(jagsnames, "^(b\\[|b\\.mm\\.|b\\.hm\\.|fix\\.)") ~ "fixed",
    stringr::str_detect(jagsnames, "^(sigma|shape|tau)") ~ "random",
    TRUE ~ "other"
  )

  out <- tibble::tibble(
    term = rt$Parameter,
    estimate = rt$mean,
    std.error = rt$sd,
    conf.low = rt$lb,
    conf.high = rt$ub,
    component = comp
  )

  # Recompute credible intervals from posterior draws for non-default levels
  if (conf.int && !isTRUE(all.equal(conf.level, 0.95))) {
    sims <- x$jags.out$BUGSoutput$sims.matrix
    if (is.null(sims)) {
      stop("Credible levels other than 0.95 require posterior draws. ",
           "Please fit the model with monitor = TRUE.", call. = FALSE)
    }
    alpha <- (1 - conf.level) / 2
    available <- jagsnames %in% colnames(sims)
    qs <- apply(sims[, jagsnames[available], drop = FALSE], 2,
                quantile, probs = c(alpha, 1 - alpha))
    out$conf.low[available] <- qs[1, ]
    out$conf.high[available] <- qs[2, ]
  }

  if (!conf.int) {
    out <- out %>% dplyr::select(-conf.low, -conf.high)
  }

  if (component != "all") {
    keep <- component
    out <- out %>% dplyr::filter(.data$component == keep)
  }

  out
}

#' Glance at a fitted bml model
#'
#' @description
#' Returns a one-row tibble with model-level information following the
#' \pkg{broom} conventions, suitable for model comparison.
#'
#' @param x A fitted model object of class \code{"bml"} returned by \code{\link{bml}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A one-row \code{tibble} with columns:
#'   \itemize{
#'     \item \code{nobs}: Number of main-level (outcome) units
#'     \item \code{n.members}: Number of member-level units across mm blocks
#'     \item \code{DIC}: Deviance Information Criterion
#'     \item \code{family}: Outcome family
#'     \item \code{n.iter}, \code{n.chains}: MCMC settings
#'   }
#'
#' @examples
#' \dontrun{
#' glance(m1)
#' }
#'
#' @seealso \code{\link{bml}}, \code{\link{tidy.bml}}, \code{\link{bmlCompare}}
#'
#' @author Benjamin Rosche \email{benrosche@@nyu.edu}
#' @export
glance.bml <- function(x, ...) {

  if (is.null(x$reg.table)) {
    stop("No estimates found. Please fit the model with run = TRUE.", call. = FALSE)
  }

  tibble::tibble(
    nobs = x$input$n.main,
    n.members = x$input$n.mm,
    DIC = as.numeric(attr(x$reg.table, "DIC")),
    family = x$input$family,
    n.iter = x$input$n.iter,
    n.chains = x$input$n.chains
  )
}
