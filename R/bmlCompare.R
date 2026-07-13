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

#' Compare fitted bml models in a publication-style table
#'
#' @description
#' Builds a side-by-side regression table from one or more fitted
#' \code{\link{bml}} models: one row per term, one column per model, with cells
#' formatted as \code{estimate (conf.low, conf.high)}. Goodness-of-fit rows
#' (e.g. \code{N} and \code{DIC}) are appended at the bottom. The result is a
#' \code{data.frame} of character cells, ready to pass to
#' \code{knitr::kable()}, \code{gt::gt()}, etc.
#'
#' For the underlying numeric values (e.g. to build custom tables or plots) use
#' the \pkg{broom} methods \code{\link{tidy.bml}} and \code{\link{glance.bml}}
#' directly; they also make \code{modelsummary::modelsummary(list(m1, m2))} work
#' on bml models.
#'
#' @param ... Fitted \code{bml} model objects, optionally named
#'   (e.g., \code{bmlCompare(base = m1, weighted = m2)}). The names become the
#'   column headers. Unnamed models are labeled with the expressions they were
#'   passed as (e.g., \code{"m1"}).
#' @param component Which parameters to include; passed to \code{\link{tidy.bml}}.
#'   One of \code{"all"} (default), \code{"fixed"}, \code{"random"}, or
#'   \code{"weights"}.
#' @param terms Optional character vector of term names selecting which rows to
#'   show and in what order (matched against the \code{term} labels from
#'   \code{\link{tidy.bml}}, e.g. \code{"smoking_alter (mm.1)"}). Terms not found
#'   are skipped with a warning. Default \code{NULL}: all terms in
#'   \code{component}, in their natural order.
#' @param labels Optional character vector renaming the rows for display. Must be
#'   the same length as the selected \code{terms} (so supply \code{terms} when
#'   using \code{labels}). Default \code{NULL}: keep the original term labels.
#' @param stats Goodness-of-fit rows to append, in order. Any of \code{"N"}
#'   (main-level units), \code{"n.members"} (member-level units), and
#'   \code{"DIC"}. Default: \code{c("N", "DIC")}. Use \code{character(0)} for none.
#' @param digits Number of decimal places for the estimate and interval bounds.
#'   Default: 3.
#' @param conf.level Width of the equal-tailed credible interval. Default: 0.95.
#'   Other levels require the models to be fitted with \code{monitor = TRUE}
#'   (see \code{\link{tidy.bml}}).
#'
#' @return A \code{data.frame} whose first column \code{Term} holds the term
#'   labels (and the requested \code{stats} labels), followed by one character
#'   column per model containing \code{estimate (conf.low, conf.high)} cells.
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' m1 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("sum"), RE = TRUE),
#'   family = weibull(),
#'   data = coalgov
#' )
#' m2 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ pseat), fn = fn("sum"), RE = TRUE),
#'   family = weibull(),
#'   data = coalgov
#' )
#'
#' bmlCompare("Equal" = m1, "Seat share" = m2)
#'
#' # Select and relabel specific rows
#' bmlCompare(
#'   "Equal" = m1, "Seat share" = m2,
#'   terms  = c("majority", "cohesion (mm.1)"),
#'   labels = c("Majority", "Cohesion")
#' ) |> knitr::kable()
#' }
#'
#' @seealso \code{\link{bml}}, \code{\link{tidy.bml}}, \code{\link{glance.bml}},
#'   \code{\link{coefPlot}}
#'
#' @author Benjamin Rosche \email{benrosche@@nyu.edu}
#' @export bmlCompare
bmlCompare <- function(..., component = "all", terms = NULL, labels = NULL,
                       stats = c("N", "DIC"), digits = 3, conf.level = 0.95) {

  models <- nameModels(list(...), sapply(substitute(list(...))[-1], deparse))
  model_names <- names(models)

  stats <- match.arg(stats, c("N", "n.members", "DIC"), several.ok = TRUE)

  # Long per-model estimates and fit statistics
  long <- purrr::imap_dfr(models, function(m, name) {
    g <- glance(m)
    tidy(m, conf.level = conf.level, component = component) %>%
      dplyr::mutate(model = name, N = g$nobs, n.members = g$n.members, DIC = g$DIC)
  })

  # Term selection and ordering
  all_terms <- unique(long$term)
  if (is.null(terms)) {
    sel_terms <- all_terms
  } else {
    missing <- setdiff(terms, all_terms)
    if (length(missing) > 0) {
      warning("Terms not found and skipped: ", paste(missing, collapse = ", "),
              call. = FALSE)
    }
    sel_terms <- terms[terms %in% all_terms]
  }
  if (length(sel_terms) == 0) {
    stop("No terms to display. Check `terms` and `component`.", call. = FALSE)
  }

  # Display labels for the rows
  disp <- sel_terms
  if (!is.null(labels)) {
    if (length(labels) != length(sel_terms)) {
      stop("`labels` must have the same length as the selected `terms` (",
           length(sel_terms), ").", call. = FALSE)
    }
    disp <- labels
  }

  # Format coefficient cells: "estimate (conf.low, conf.high)"
  fmt <- paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)")
  long$cell <- sprintf(fmt, long$estimate, long$conf.low, long$conf.high)

  out <- data.frame(Term = disp, stringsAsFactors = FALSE, check.names = FALSE)
  for (mn in model_names) {
    out[[mn]] <- vapply(sel_terms, function(tm) {
      v <- long$cell[long$model == mn & long$term == tm]
      if (length(v) == 0) "\u2014" else v[1]   # em dash for terms absent from a model
    }, character(1))
  }

  # Goodness-of-fit rows
  if (length(stats) > 0) {
    stat_rows <- lapply(stats, function(s) {
      vals <- vapply(model_names, function(mn) {
        gv <- long[[s]][long$model == mn][1]
        if (s == "DIC") formatC(gv, format = "f", digits = 0)   # round DIC to integer
        else format(gv, big.mark = "", trim = TRUE)
      }, character(1))
      row <- as.data.frame(c(list(Term = s), as.list(vals)),
                           stringsAsFactors = FALSE, check.names = FALSE)
      names(row) <- c("Term", model_names)
      row
    })
    out <- rbind(out, do.call(rbind, stat_rows))
  }

  rownames(out) <- NULL
  class(out) <- c("bmlCompare", "data.frame")
  out
}

#' @exportS3Method base::print
print.bmlCompare <- function(x, ...) {
  print(`class<-`(x, "data.frame"), row.names = FALSE, ...)
  invisible(x)
}

#' @exportS3Method knitr::knit_print
knit_print.bmlCompare <- function(x, ...) {
  knitr::knit_print(knitr::kable(`class<-`(x, "data.frame")), ...)
}
