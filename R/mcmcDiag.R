#' Summarize MCMC convergence diagnostics
#'
#' Computes common convergence diagnostics for selected parameters from a
#' JAGS/BUGS fit and returns a compact, report-ready table. The diagnostics
#' include Gelman–Rubin \eqn{\hat{R}}, Geweke z-scores, Heidelberger-Welch
#' stationarity p-values, and autocorrelation at lag 50.
#'
#' @param bml.out A model fit object containing JAGS output, typically as returned
#'   by \code{R2jags::jags()}, with component \code{$jags.out$BUGSoutput}.
#' @param parameters Character vector of parameter names (or patterns) to extract.
#'   These may be exact names or patterns (e.g., a prefix like \code{"b"} that
#'   matches \code{"b[1]"}, \code{"b[2]"}, …).
#'
#' @details
#' Internally, the function converts the BUGS/JAGS output to a
#' \code{coda::mcmc.list}, then computes per-chain
#' diagnostics and averages them across chains for each parameter:
#' \itemize{
#'   \item \strong{Gelman–Rubin} (\eqn{\hat{R}}): \code{coda::gelman.diag()}.
#'         Values close to 1 indicate convergence; a common heuristic is
#'         \eqn{\hat{R} \le 1.1}.
#'   \item \strong{Geweke} z-score: \code{coda::geweke.diag()}.
#'         Large absolute values (e.g., \eqn{|z|>2}) suggest lack of convergence.
#'   \item \strong{Heidelberger–Welch} p-value: \code{coda::heidel.diag()} tests
#'         the null of stationarity in the chain segment.
#'   \item \strong{Autocorrelation (lag 50)}: \code{coda::autocorr()} at lag 50,
#'         averaged across chains.
#' }
#' All statistics are rounded to three decimals. The returned table is transposed
#' so that \emph{rows are diagnostics} and \emph{columns are parameters}.
#'
#' @return A \code{data.frame} with one row per diagnostic and one column per
#'   parameter; cell entries are the average diagnostic values across chains.
#'   Row names include: \code{"Gelman/Rubin convergence statistic"},
#'   \code{"Geweke z-score"}, \code{"Heidelberger/Welch p-value"},
#'   \code{"Autocorrelation (lag 50)"}.
#'
#' @seealso \code{\link[coda]{gelman.diag}},
#'   \code{\link[coda]{geweke.diag}}, \code{\link[coda]{heidel.diag}},
#'   \code{\link[coda]{autocorr}}
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' # Fit model
#' m1 <- bml(
#'   Surv(govdur, earlyterm) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE),
#'   family = "Weibull",
#'   monitor = TRUE,
#'   data = coalgov
#' )
#'
#' # Check convergence for main parameters
#' mcmcDiag(m1, parameters = "b")  # All b coefficients
#'
#' # Check specific parameters
#' mcmcDiag(m1, parameters = c("b[1]", "b[2]", "shape"))
#'
#' # Check mm block parameters
#' mcmcDiag(m1, parameters = c("b.mm.1", "sigma.mm"))
#'
#' # Interpreting results:
#' # - Gelman-Rubin < 1.1: Good convergence
#' # - |Geweke z| < 2: No evidence against convergence
#' # - Heidelberger p > 0.05: Chain appears stationary
#' # - Low autocorrelation at lag 50: Good mixing
#' }
#'
#' @author Benjamin Rosche \email{benrosche@@nyu.edu}
#' @export mcmcDiag

mcmcDiag <- function(bml.out, parameters) {

  # Check that JAGS output is available ---------------------------------------------------------- #

  if (is.null(bml.out$jags.out) || is.null(bml.out$jags.out$BUGSoutput) ||
      is.null(bml.out$jags.out$BUGSoutput$sims.array)) {
    stop("JAGS output could not be retrieved. Please ensure that monitor = TRUE when fitting the model.", call. = FALSE)
  }

  # Extract mcmc.list from bml.out --------------------------------------------------------------- #

  crMCMC <- function(bml.out, parameter, regex=T) {
    m <- coda::as.mcmc.list(bml.out$jags.out$BUGSoutput)
    vars <- colnames(m[[1]])
    
    # Escape regex metacharacters
    esc <- function(x) gsub("([][{}()+*^$.|\\?\\\\])", "\\\\\\1", x)
    
    get_matches <- function(p) {
      if (!regex) return(vars[vars %in% p])
      
      # exact match if full name present
      if (p %in% vars) return(p)
      
      # if contains brackets → treat literally
      if (grepl("\\[.*\\]", p)) {
        return(vars[grepl(paste0("^", esc(p), "$"), vars)])
      }
      
      # otherwise treat as prefix (match all indexed variants)
      vars[grepl(paste0("^", esc(p), "\\[[^]]+\\]$"), vars)]
    }
    
    sel <- unique(unlist(lapply(parameter, get_matches)))
    
    if (length(sel) == 0)
      stop("No parameters matched any pattern.", call. = FALSE)
    
    coda::mcmc.list(lapply(m, \(ch) ch[, sel, drop = FALSE]))
  }
  
  
  mcmcl <- crMCMC(bml.out, parameters)
  
  n.chains <- length(mcmcl)
  n.parameters <- mcmcl[[1]] %>% ncol()
  
  # Save different convergence statistics -------------------------------------------------------- #
  
  message("Parameter(s): ", parameters)
  
  gelman_rubin <- 
    coda::gelman.diag(mcmcl, autoburnin = F)$psrf %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("Parameter") %>%
    select(Parameter, "Gelman/Rubin convergence statistic"=2) 
  
  # Gelman/Rubin: "A rule of thumb is that values of 1.1 and less suggests adequate convergence" - Finley (2013): Using JAGS in R with the rjags package
  
  geweke <- 
    coda::geweke.diag(mcmcl) %>% 
    purrr::map(., \(x) x[[1]]) %>%
    as.data.frame(col.names = paste0("Chain ", 1:n.chains)) %>%
    tibble::rownames_to_column("Parameter") %>%
    rowwise() %>%
    summarise(
      Parameter,
      "Geweke z-score" = mean(c_across(starts_with("Chain")), na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Geweke: "The test statistic is a z-score, so |Z| > 2 indicates poor convergence" - Reich (XXXX): Applied Bayesian Analysis
  
  heidel <- 
    coda::heidel.diag(mcmcl) %>% 
    purrr::map(., .f=\(x){ setNames(x[, "pvalue"], rownames(x)) }) %>% 
    as.data.frame(col.names = paste0("Chain ", 1:n.chains)) %>%
    tibble::rownames_to_column("Parameter") %>%
    rowwise() %>%
    summarise(
      Parameter,
      "Heidelberger/Welch p-value" = mean(c_across(starts_with("Chain")), na.rm = TRUE)
    ) %>%
    ungroup()
  
  # "The Heidelberger and Welch diagnostic first tests the null hypothesis that the Markov Chain is 
  # in the stationary distribution and produces p-values for each estimated parameter"
  
  autocorr <- 
    coda::autocorr(mcmcl) %>% 
    purrr::map(., .f=\(x){ if(n.parameters>1) setNames(diag(x[5,,]), colnames(x)) else setNames(x[5,,], colnames(x)) }) %>% 
    as.data.frame(col.names = paste0("Chain ", 1:n.chains)) %>%
    tibble::rownames_to_column("Parameter") %>%
    rowwise() %>%
    summarise(
      Parameter,
      "Autocorrelation (lag 50)" = mean(c_across(starts_with("Chain")), na.rm = TRUE)
    ) %>%
    ungroup()
  
  return(
    gelman_rubin %>%
      left_join(geweke, by = "Parameter") %>%
      left_join(heidel, by = "Parameter") %>%
      left_join(autocorr, by = "Parameter") %>%
      mutate(across(-Parameter, ~ round(., 3))) %>%
      tibble::column_to_rownames("Parameter") %>%
      t() %>%
      as.data.frame() 
  )
  
}
