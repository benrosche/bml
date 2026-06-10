#' Coefficient plot for fitted bml models
#'
#' @description
#' Creates a coefficient ("dot-and-whisker") plot from one or more fitted
#' \code{\link{bml}} models: posterior means as points and credible intervals
#' as horizontal ranges, with a dashed reference line at zero. When several
#' models are supplied, their coefficients are placed next to each other for
#' each term so the models can be compared directly.
#'
#' @param ... One or more fitted \code{bml} model objects, optionally named
#'   (e.g., \code{coefPlot(base = m1, weighted = m2)}). Names are used in the
#'   legend; unnamed models are labeled with the expressions they were passed
#'   as (e.g., \code{"m1"}).
#' @param parameters Optional character vector of term labels (as shown in
#'   \code{summary()}) to include, in the order they should appear (top to
#'   bottom). Default: all terms of the selected \code{component}.
#' @param intercept Logical. Include the intercept? Default: \code{TRUE}.
#' @param conf.level Width of the equal-tailed credible interval. Default: 0.95.
#'   Other levels require the models to be fitted with \code{monitor = TRUE}
#'   (see \code{\link{tidy.bml}}).
#' @param component Which parameters to plot; passed to \code{\link{tidy.bml}}.
#'   Default: \code{"fixed"} (regression coefficients).
#'
#' @return A \code{ggplot} object that can be further customized with
#'   \pkg{ggplot2} layers or saved with \code{ggplot2::ggsave()}.
#'
#' @details
#' The plot uses a colorblind-friendly palette (Okabe-Ito) to distinguish
#' models and intentionally substantial line widths and point sizes for
#' readability in presentations and publications. All defaults can be
#' overridden by adding \pkg{ggplot2} layers to the returned object, e.g.
#' \code{coefPlot(m1) + ggplot2::scale_colour_brewer(palette = "Dark2")}.
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
#' # Single model
#' coefPlot(m1)
#'
#' # Without intercept
#' coefPlot(m1, intercept = FALSE)
#'
#' # Compare two models side by side
#' coefPlot(equal = m1, seatshare = m2)
#'
#' # Save the plot
#' p <- coefPlot(m1, m2)
#' ggplot2::ggsave("coefplot.pdf", p, width = 6, height = 4)
#' }
#'
#' @seealso \code{\link{bml}}, \code{\link{tidy.bml}}, \code{\link{bmlCompare}},
#'   \code{\link{monetPlot}}
#'
#' @author Benjamin Rosche \email{benrosche@@nyu.edu}
#' @export coefPlot
coefPlot <- function(..., parameters = NULL, intercept = TRUE,
                     conf.level = 0.95, component = "fixed") {

  models <- nameModels(list(...), sapply(substitute(list(...))[-1], deparse))
  n.models <- length(models)

  # Build plot data ------------------------------------------------------------------------------ #

  plot.data <- purrr::imap_dfr(models, function(m, name) {
    tidy(m, conf.level = conf.level, component = component) %>%
      dplyr::mutate(model = name)
  })

  if (!intercept) {
    plot.data <- plot.data %>% dplyr::filter(.data$term != "Intercept")
  }

  if (!is.null(parameters)) {
    missing_terms <- setdiff(parameters, plot.data$term)
    if (length(missing_terms) > 0) {
      stop("Term(s) not found in the model(s): ", paste(missing_terms, collapse = ", "),
           ". Use tidy() to see available terms.", call. = FALSE)
    }
    plot.data <- plot.data %>% dplyr::filter(.data$term %in% parameters)
    term_order <- parameters
  } else {
    term_order <- unique(plot.data$term)
  }

  if (nrow(plot.data) == 0) {
    stop("No terms left to plot.", call. = FALSE)
  }

  # Reverse term order so the first term appears at the top after coord_flip();
  # reverse model levels so the first model appears on top within each term
  plot.data <- plot.data %>%
    dplyr::mutate(
      term = factor(.data$term, levels = rev(term_order)),
      model = factor(.data$model, levels = rev(names(models)))
    )

  # Create plot ----------------------------------------------------------------------------------- #

  # Okabe-Ito colorblind-friendly palette
  palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
               "#E69F00", "#56B4E9", "#F0E442", "#999999")
  palette <- rep_len(palette, n.models)

  position <- if (n.models > 1) ggplot2::position_dodge(width = 0.6) else "identity"

  p <-
    ggplot2::ggplot(plot.data, ggplot2::aes(x = term, y = estimate, colour = model)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      position = position,
      linewidth = 0.9,
      size = 0.7,
      fatten = 2.5
    ) +
    ggplot2::scale_colour_manual(
      values = setNames(palette, names(models)),
      breaks = names(models)  # legend in the order the models were passed
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "",
      y = paste0("Posterior mean (", round(conf.level * 100), "% credible interval)"),
      colour = ""
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.position = if (n.models > 1) "bottom" else "none"
    )

  return(p)
}
