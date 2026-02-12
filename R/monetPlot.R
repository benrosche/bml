#' Visualize posterior distributions with density and trace plots
#'
#' @description
#' Creates a combined diagnostic plot showing both the posterior density and MCMC
#' trace plot for a specified parameter. Helps assess convergence and visualize
#' posterior uncertainty. The plot displays the median and 90\% highest posterior
#' density (HPD) interval.
#'
#' @param bml A fitted model object of class \code{"bml"} returned by \code{\link{bml}}.
#'   Must be fitted with \code{monitor = TRUE} to store MCMC chains.
#'
#' @param parameter Character string specifying the parameter to plot. Must use the
#'   internal parameter name (i.e., row names from \code{bml$reg.table}). Examples:
#'   \code{"b[1]"} (intercept), \code{"b[2]"} (first covariate), \code{"b.mm.1"}
#'   (first mm block coefficient), \code{"sigma.mm"} (mm random effect SD).
#'
#' @param label Optional character string for the parameter label displayed on the
#'   plot. If \code{NULL} (default), uses the internal parameter name.
#'
#' @param r Number of decimal places for displayed quantiles and statistics.
#'   Default: 2.
#'
#' @param yaxis Logical; if \code{TRUE} (default), display axis titles ("Density"
#'   and "Scans"). If \code{FALSE}, omit axis titles for cleaner appearance when
#'   combining multiple plots.
#'
#' @return A \code{ggplot} object (using \code{patchwork}) combining two panels:
#'   \itemize{
#'     \item \strong{Top panel}: Posterior density with shaded 90\% HPD interval.
#'       Solid vertical line at zero, dashed line at posterior median.
#'     \item \strong{Bottom panel}: Trace plot showing MCMC iterations across chains.
#'       Same reference lines as top panel. Helps diagnose convergence and mixing.
#'   }
#'
#' @details
#' \strong{Interpreting the Plot:}
#' \itemize{
#'   \item \strong{Density panel}: Shows the posterior distribution. The dashed line
#'     marks the median (central estimate). Shading indicates the 90\% credible region.
#'   \item \strong{Trace panel}: Shows parameter values across MCMC iterations for
#'     each chain. Good mixing looks like "fuzzy caterpillars" with chains overlapping.
#'     Poor mixing shows trends, stickiness, or separation between chains.
#' }
#'
#' \strong{Convergence Checks:}
#' \itemize{
#'   \item Chains should overlap and explore the same space
#'   \item No sustained trends or drift
#'   \item Rapid mixing (no long autocorrelation)
#' }
#'
#' Use \code{\link{mcmcDiag}} for formal convergence statistics (Gelman-Rubin,
#' Geweke, etc.).
#'
#' @seealso \code{\link{bml}}, \code{\link{mcmcDiag}}, \code{\link{summary.bml}}
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' # Fit model with monitoring enabled
#' m1 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE) +
#'     hm(id = id(cid), type = "RE"),
#'   family = "Weibull",
#'   monitor = TRUE,  # Required for monetPlot
#'   data = coalgov
#' )
#'
#' # Plot intercept
#' monetPlot(m1, parameter = "b[1]", label = "Intercept")
#'
#' # Plot majority coefficient with custom label
#' monetPlot(m1, parameter = "b[2]", label = "Majority Government Effect")
#'
#' # Plot mm coefficient
#' monetPlot(m1, parameter = "b.mm.1", label = "Party Fragmentation")
#'
#' # Plot random effect SD
#' monetPlot(m1, parameter = "sigma.mm.1")
#'
#' # List available parameters
#' rownames(m1$reg.table)
#' }
#'
#' @export monetPlot
#' @author Benjamin Rosche <benrosche@@nyu.edu>

monetPlot <- function(bml, parameter, label=NULL, r=2, yaxis=T) {

  # Checks --------------------------------------------------------------------------------------- #

  if (is.null(bml$jags.out) || is.null(bml$jags.out$BUGSoutput) ||
      is.null(bml$jags.out$BUGSoutput$sims.array)) {
    stop("JAGS output could not be retrieved. Please ensure that monitor = TRUE when fitting the model.", call. = FALSE)
  }
  if(is.null(label)) label = parameter

  # Get mcmclist and posterior stats ------------------------------------------------------------- #

  # Build ggs-compatible data frame manually to avoid ggmcmc::ggs() version issues
  sims <- bml$jags.out$BUGSoutput$sims.array
  param_names <- dimnames(sims)[[3]]
  param_idx <- which(param_names == parameter)
  if (length(param_idx) == 0) {
    stop(paste0("Parameter '", parameter, "' not found. Use rownames(bml$reg.table) to see available parameters."), call. = FALSE)
  }

  n_iter <- dim(sims)[1]
  n_chains <- dim(sims)[2]
  mcmc.ggs <- do.call(rbind, lapply(seq_len(n_chains), function(ch) {
    dplyr::tibble(
      Iteration = seq_len(n_iter),
      Chain     = factor(ch),
      Parameter = factor(label),
      value     = sims[, ch, param_idx]
    )
  }))
  attr(mcmc.ggs, "nChains")     <- n_chains
  attr(mcmc.ggs, "nParameters") <- 1L
  attr(mcmc.ggs, "nIterations") <- n_iter
  attr(mcmc.ggs, "nBurnin")     <- 0L
  attr(mcmc.ggs, "nThin")       <- 1L
  attr(mcmc.ggs, "description") <- "bml"

  p.quantiles <- round(quantile(dplyr::pull(mcmc.ggs, value), c(.05, .5, .95)), r)
  p.mad <- round(mad(dplyr::pull(mcmc.ggs, value)), r)

  # Calculate x-axis limits to trim very low density tails
  # Use 0.1% and 99.9% quantiles with small padding
  xlim_quantiles <- quantile(dplyr::pull(mcmc.ggs, value), c(0.001, 0.999))
  xlim_range <- diff(xlim_quantiles)
  xlim_limits <- c(xlim_quantiles[1] - 0.05 * xlim_range, xlim_quantiles[2] + 0.05 * xlim_range)

  # Create plots --------------------------------------------------------------------------------- #

  if(isTRUE(yaxis)) {
    yaxis <-  "Density"
    xaxis <- paste0("Scans (", n_chains, " chains)")
  } else {
    yaxis <- ""
    xaxis <- ""
  }

  # Density plot
  p1 <-
    ggmcmc::ggs_density(mcmc.ggs, hpd = TRUE) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_vline(xintercept = p.quantiles[2], linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = xlim_limits) +
    ggplot2::labs(x="", y=yaxis, title = paste0("Parameter: ", label)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title.y = ggplot2::element_text(face = "bold", size = 14)
    )

  # Traceplot
  p2 <-
    ggmcmc::ggs_traceplot(mcmc.ggs, original_burnin = FALSE) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(yintercept = p.quantiles[2], linetype = "dashed") +
    ggplot2::coord_flip(ylim = xlim_limits) +
    ggplot2::scale_y_continuous(
      breaks = sort(as.numeric(c(p.quantiles, 0))),
      labels = function(x) {
        labs <- c(
          setNames(
            paste0(as.character(p.quantiles), "\n(", names(p.quantiles), ")"),
            as.character(p.quantiles)
          ),
          "0" = "0"
        )
        labs[as.character(x)]
      }
    ) +
    ggplot2::labs(
      x = xaxis,
      y = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.margin = ggplot2::unit(c(-0.6, 0, 0, 0), "cm"),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(face = "bold", size = 14),
      axis.text.x  = ggplot2::element_text(color = "black", size = 12)
    )

  # Use wrap_plots to avoid S7 operator dispatch issues with /
  return(
    patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 1)) +
      patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")))
  )

}
