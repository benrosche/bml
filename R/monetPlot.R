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
#'   Surv(govdur, earlyterm) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE) +
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
#' monetPlot(m1, parameter = "sigma.mm")
#'
#' # List available parameters
#' rownames(m1$reg.table)
#' }
#'
#' @export monetPlot
#' @author Benjamin Rosche <benrosche@@nyu.edu>

monetPlot <- function(bml, parameter, label=NULL, r=2, yaxis=T) {

  escape_regex <- function(x) gsub("([][{}()+*^$|.?\\\\])", "\\\\\\1", x)
  
  # Checks --------------------------------------------------------------------------------------- #
  
  if(is.null(bml$jags.out)) stop("JAGS output could not be retrieved. Please specify monitor = T when running bml.")
  if(is.null(label)) label = parameter
  
  # Get mcmclist and posterior stats ------------------------------------------------------------- #
  
  mcmc.list <- coda::as.mcmc(bml$jags.out)
  
  mcmc.ggs <- 
    ggmcmc::ggs(
      mcmc.list,
      family = paste0("^", escape_regex(parameter), "$"), 
      par_labels = data.frame(Parameter = parameter, Label = label)
    )
  
  p.quantiles <- round(quantile(dplyr::pull(mcmc.ggs, value), c(.05, .5, .95)), r)
  p.mad <- round(mad(dplyr::pull(mcmc.ggs, value)), r)

  # Create plots --------------------------------------------------------------------------------- #

  if(isTRUE(yaxis)) {
    yaxis <-  "Density"
    xaxis <- paste0("Scans (", max(dplyr::pull(mcmc.ggs, Chain)), " chains)") 
  } else {
    yaxis <- ""
    xaxis <- ""
  }
  
  # Density plot
  p1 <-
    ggmcmc::ggs_density(mcmc.ggs, hpd = TRUE) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_vline(xintercept = p.quantiles[2], linetype = "dashed") +
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
    ggplot2::coord_flip() +
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

  return(
    (p1 / p2) +
      patchwork::plot_layout(heights = c(1, 1)) +
      patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")))
  )
  
}
