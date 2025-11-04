#' @title Create monetPlots for your bml results.
#'
#' @description The function \strong{monetPlot} creates a density plot of the posterior 
#' distribution of your model parameters and the traceplot that led to this density.
#' 
#' @param bml A bml object. bml has to be run with monitor=T
#' @param parameter A string with the parameter name. The internal name has to be used, which are the rownames in the bml reg.table output.
#' @param label String to describe the parameter on the graph's x-axis. Optional. If not specified, the internal parameter name is used.
#' @param r Specify number of decimal places. Default equals 3.
#' @param yaxis Logical. If FALSE, the y-axis title is omitted.
#' @return Returns a plot. The solid vertical is at 0 and the dashed vertical line is the mode of the posterior distributions.
#'
#' @examples data(coalgov)
#' m1 <- bml(Surv(govdur, earlyterm) ~ 1 + majority + mm(id(pid, gid), mmc(fdep), mmw(w ~ 1/n, constraint=T)) + hm(id=cid, name=cname, type=RE, showFE=F),
#'           family="Weibull", monitor=T, data=coalgov)
#' monetPlot(m1, parameter="b.l1")
#'
#' @export monetPlot
#' @author Benjamin Rosche <benrosche@@nyu.edu>

monetPlot <- function(bml, parameter, label=NULL, r=2, yaxis=T) {
  
  library(ggplot2)
  library(ggmcmc)
  library(coda)
  library(patchwork)
  
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
  
  p.quantiles <- mcmc.ggs %>% pull(value) %>% quantile(., c(.05, .5, .95)) %>% round(., r)
  p.mad <- mcmc.ggs %>% pull(value) %>% mad() %>% round(., r)
  
  # Create plots --------------------------------------------------------------------------------- #
  
  if(isTRUE(yaxis)) {
    yaxis <-  "Density"
    xaxis <- paste0("Scans (", mcmc.ggs %>% pull(Chain) %>% max(), " chains)") 
  } else {
    yaxis <- ""
    xaxis <- ""
  }
  
  # Density plot
  p1 <-
    ggs_density(mcmc.ggs, hpd = TRUE) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = p.quantiles[2], linetype = "dashed") +
    labs(x="", y=yaxis, title = paste0("Parameter: ", label)) +
    theme_minimal() +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.text = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title.y = element_text(face = "bold", size = 14)
    )
  
  # Traceplot
  p2 <-
    ggs_traceplot(mcmc.ggs, original_burnin = FALSE) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = p.quantiles[2], linetype = "dashed") +
    coord_flip() +
    scale_y_continuous(
      breaks = c(p.quantiles, 0) %>% as.numeric() %>% sort(),
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
    labs(
      x = xaxis,
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.margin = unit(c(-0.6, 0, 0, 0), "cm"),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text.x  = element_text(color = "black", size = 12)  
    )
  
  return(
    (p1 / p2) +
      plot_layout(heights = c(1, 1)) +  
      plot_annotation(theme = theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))) 
  )
  
}
