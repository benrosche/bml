#' @title Bayesian Multiple-Membership Multilevel Models with Parameterizable Weight Functions Using JAGS
#'
#' @description
#' The \strong{bml} package provides a user-friendly interface for fitting Bayesian multiple-membership
#' multilevel models with parameterizable weight functions via JAGS.
#'
#' JAGS must be installed separately: \url{https://sourceforge.net/projects/mcmc-jags/}.
#'
#' @details
#' In addition to hierarchical and cross-classified multilevel models, the \strong{bml} package allows
#' users to fit Bayesian multiple-membership models. Unlike tools such as
#' \code{\link[brms]{brms}} or MLwiN (\url{https://www.bristol.ac.uk/cmm/software/mlwin/}),
#' \strong{bml} lets users specify and estimate models in which membership weights are parameterized
#' through flexible formula syntax. This enables a more nuanced examination of how effects from
#' member-level units aggregate to group level (the micro-macro link).
#'
#' The package automatically generates JAGS code to fit the model and processes the output
#' to facilitate interpretation of model parameters and diagnostics.
#'
#' For accessible introductions to multiple-membership models, see Fielding and Goldstein (2006)
#' and Beretvas (2010). Advanced treatments include Goldstein (2011, Ch. 13),
#' Rasbash and Browne (2001, 2008), Browne et al. (2001), and Leckie (2013).
#' The package and modeling framework are introduced in:
#' Rosche, B. (2025). \emph{A Multilevel Model for Coalition Governments: Uncovering Party-Level
#' Dependencies Within and Between Governments}. \emph{Political Analysis}.
#'
#' @section Formula Components:
#' \itemize{
#'   \item \strong{Outcome (Y):} The dependent variable. For survival models, use \code{Surv(time, event)}.
#'   \item \strong{Intercept (1):} Includes an intercept term; use \code{0} to omit it.
#'   \item \strong{Main-level predictors (X.main):} Variables defined at the main (group) level, separated by \code{+}.
#'   \item \strong{HM-level predictors (X.hm):} Variables defined at the nesting level, separated by \code{+}.
#'   \item \strong{Multiple membership object (\code{mm()}):} Defines how member-level units are
#'         associated with group-level constructs using a user-specified weighting function.
#'         Multiple \code{mm()} objects can be specified with different weight functions.
#'   \item \strong{Hierarchical membership (\code{hm()}):} Specifies nesting of main-level units within
#'         higher-level entities. Cross-classified structures can be modeled by including multiple \code{hm()} objects.
#' }
#'
#' \strong{Important:} The formula parser does not support \code{I()} inside \code{bml()}.
#' Create transformations and interactions in your data object before modeling.
#'
#' @section Multiple Membership Object \code{mm()}:
#' \preformatted{
#' mm(
#'   id   = id(mmid, mainid),
#'   vars = vars(X.mm),
#'   fn   = fn(w ~ 1/n, c = TRUE, ar = FALSE),
#'   RE   = TRUE
#' )
#' }
#'
#' \strong{Components:}
#' \itemize{
#'   \item \code{id(mmid, mainid)}: Specifies identifiers linking each member-level unit (\code{mmid})
#'         to its corresponding group-level entities (\code{mainid}).
#'   \item \code{vars(X.mm)}: Specifies member-level covariates aggregated across memberships.
#'         Use \code{+} to include multiple variables. Set to \code{NULL} for RE-only blocks.
#'   \item \code{fn(w ~ ..., c, ar)}: Defines the weight function (micro-macro link).
#'   \item \code{RE}: Logical; if \code{TRUE}, include random effects for this block.
#'         Automatically \code{TRUE} if \code{vars = NULL}.
#' }
#'
#' \strong{Multiple mm() blocks:}
#' You can specify multiple \code{mm()} blocks with different weight functions:
#' \preformatted{
#' mm(id = id(pid, gid), vars = vars(X1), fn = fn(w ~ 1/n), RE = FALSE) +
#' mm(id = id(pid, gid), vars = vars(X2), fn = fn(w ~ max(n)), RE = FALSE) +
#' mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE)
#' }
#'
#' @section Hierarchical Membership Object \code{hm()}:
#' \preformatted{
#' hm(id = id(hmid), vars = vars(X.hm), name = hmname, type = "RE", showFE = FALSE)
#' }
#'
#' \strong{Components:}
#' \itemize{
#'   \item \code{id = id(hmid)}: Variable identifying nesting-level groups.
#'   \item \code{vars = vars(X.hm)}: Nesting-level variables, or \code{NULL}.
#'   \item \code{name = hmname}: Optional labels for nesting-level units.
#'   \item \code{type}: \code{"RE"} (default) or \code{"FE"}.
#'   \item \code{showFE}: If \code{TRUE} and \code{type = "FE"}, report the fixed effects.
#' }
#'
#' @section Supported Families / Links:
#' \itemize{
#'   \item Gaussian (continuous): \code{family = "Gaussian"}
#'   \item Binomial (logistic): \code{family = "Binomial"}
#'   \item Weibull survival: \code{family = "Weibull"}, outcome: \code{Surv(time, event)}
#'   \item Cox survival: \code{family = "Cox"}, outcome: \code{Surv(time, event)}
#' }
#'
#' @section Priors:
#' Priors can be specified for parameters. With multiple mm() blocks, use indexed names:
#' \preformatted{
#' priors = list(
#'   "b.mm.1 ~ dnorm(0, 0.01)",
#'   "b.w.1 ~ dnorm(0, 0.1)",
#'   "tau.mm ~ dscaled.gamma(25, 1)"
#' )
#' }
#'
#' @param formula A symbolic model formula.
#' @param family A character string specifying the model family.
#' @param priors A named list or character vector mapping parameter groups to JAGS priors.
#' @param inits A list of initial values used for all chains.
#' @param n.iter Total number of MCMC iterations.
#' @param n.burnin Number of burn-in iterations to discard.
#' @param n.thin Thinning rate for MCMC sampling.
#' @param chains Number of MCMC chains to run.
#' @param seed Random seed for reproducibility.
#' @param run Logical; if \code{TRUE} (default), JAGS is executed.
#' @param parallel Logical; if \code{TRUE}, run chains in parallel.
#' @param monitor Logical; if \code{TRUE}, store additional outputs.
#' @param modelfile Character or logical for saving/loading JAGS model.
#' @param data A data frame where each row represents a member-level observation.
#'
#' @return A list of class "bml" containing model outputs.
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' # Single mm() block
#' m1 <- bml(
#'   Surv(govdur, earlyterm) ~ 1 + majority +
#'     mm(
#'       id   = id(pid, gid),
#'       vars = vars(fdep),
#'       fn   = fn(w ~ 1/n, c = TRUE),
#'       RE   = TRUE
#'     ) +
#'     hm(id = id(cid), type = "RE"),
#'   family  = "Weibull",
#'   data    = coalgov
#' )
#'
#' # Multiple mm() blocks with different weight functions
#' m2 <- bml(
#'   Y ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE) +
#'     mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE) +
#'     hm(id = id(cid), type = "RE"),
#'   family  = "Gaussian",
#'   data    = coalgov
#' )
#' }
#'
#' @export
#' @author Benjamin Rosche <benrosche@@nyu.edu>

bml <- function(
  formula,
  family = "Gaussian",
  priors = NULL,
  inits = NULL,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = max(1, floor((n.iter - n.burnin) / 1000)),
  chains = 3,
  seed = NULL,
  run = TRUE,
  parallel = FALSE,
  monitor = TRUE,
  modelfile = FALSE,
  data = NULL
) {

  # Test call:
  # formula <- sim.y ~ 1 + majority + mm( id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1 / n, c = TRUE), RE = TRUE ) + hm(id = id(cid), vars = NULL, type = "RE"); family = "Gaussian"; run=F; data = coalgov

  # ========================================================================================== #
  # 0. Checks
  # ========================================================================================== #

  if (is.null(data)) {
    stop("No data supplied.")
  }

  # ========================================================================================== #
  # 1. Dissect formula
  # ========================================================================================== #

  DIR <- system.file(package = "bml")

  formula_parts <- dissectFormula(formula, family, data)

  mm <- formula_parts$mm
  hm <- formula_parts$hm

  has_mm <- length(mm) > 0
  has_hm <- length(hm) > 0

  # ========================================================================================== #
  # 2. Create data structures
  # ========================================================================================== #

  data_parts <- createData(data, formula_parts)

  data      <- data_parts$data
  mm_blocks <- data_parts$mm_blocks
  main      <- data_parts$main
  hm_blocks <- data_parts$hm_blocks

  # ========================================================================================== #
  # 3. Create/edit JAGS modelstring
  # ========================================================================================== #

  modelstring <- editModelstring(
    family,
    priors,
    mm_blocks,
    main,
    hm_blocks,
    mm,
    hm,
    DIR,
    monitor,
    modelfile
  )

  # Save or read modelstring
  if (isTRUE(modelfile)) {
    modelfile_path <- file.path(getwd(), "modelstring.txt")
    readr::write_file(modelstring, modelfile_path)
    message("JAGS model saved to: ", modelfile_path)
  } else if (
    !isFALSE(modelfile) && length(modelfile) > 0 && is.character(modelfile)
  ) {
    tryCatch(
      {
        modelstring <- readr::read_file(modelfile)
      },
      error = function(e) {
        stop("Could not find/read model file in ", modelfile)
      }
    )
  }

  # ========================================================================================== #
  # 4. Transform data into JAGS format
  # ========================================================================================== #

  jags_vars <- createJagsVars(
    data,
    family,
    mm_blocks,
    main,
    hm_blocks,
    mm,
    hm,
    monitor,
    modelfile,
    chains,
    inits
  )

  ids <- jags_vars$ids
  Ns <- jags_vars$Ns
  Xs <- jags_vars$Xs
  Ys <- jags_vars$Ys
  jags.params <- jags_vars$jags.params
  jags.inits <- jags_vars$jags.inits
  jags.data <- jags_vars$jags.data

  # ========================================================================================== #
  # 5. Run JAGS
  # ========================================================================================== #

  if (run) {
    # Get seed
    if (is.null(seed)) {
      seed <- round(runif(1, 0, 1000))
    }

    if (parallel) {
      # Run parallel -------------------------------------------------------------------- #

      parallelfile <- tempfile(fileext = ".jags")
      on.exit(unlink(parallelfile), add = TRUE)
      writeLines(modelstring, parallelfile)
      jags.out <- do.call(
        R2jags::jags.parallel,
        list(
          data = jags.data,
          inits = jags.inits[1],
          n.chains = chains,
          parameters.to.save = jags.params,
          n.iter = n.iter,
          n.burnin = n.burnin,
          n.thin = n.thin,
          jags.seed = seed,
          model.file = parallelfile
        )
      )
    } else {
      # Run sequentially ---------------------------------------------------------------- #

      set.seed(seed)
      jags.out <- R2jags::jags(
        data = jags.data,
        inits = jags.inits,
        n.chains = chains,
        parameters.to.save = jags.params,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin,
        model.file = textConnection(modelstring)
      )
    }

    # Format JAGS output ---------------------------------------------------------------- #

    formatted <- formatJags(
      jags.out,
      monitor,
      Ns,
      mm_blocks,
      main,
      hm_blocks,
      mm,
      hm
    )

    reg.table <- formatted$reg.table
    w <- formatted$w
    re.mm <- formatted$re.mm
    re.hm <- formatted$re.hm
    pred <- formatted$pred

    # Prepare additional information ---------------------------------------------------- #

    # Collect mm block info
    mm_info <- lapply(seq_along(mm), function(k) {
      list(
        vars = mm_blocks[[k]]$vars,
        fn = mm_blocks[[k]]$fn,
        RE = mm_blocks[[k]]$RE
      )
    })

    # Collect hm block info
    hm_info <- if (has_hm) {
      lapply(seq_along(hm), function(k) {
        list(
          id = hm_blocks[[k]]$id,
          vars = hm_blocks[[k]]$vars,
          type = hm_blocks[[k]]$type,
          showFE = hm_blocks[[k]]$showFE
        )
      })
    } else {
      list()
    }

    # Save info on input
    input <- list(
      family = family,
      priors = priors,
      inits = inits,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      chains = chains,
      parallel = parallel,
      seed = seed,
      monitor = monitor,
      modelfile = modelfile,
      run = run,
      lhs = main$lhs,
      mainvars = main$vars,
      mm = mm_info,
      hm = hm_info,
      n.umm = Ns$n.umm,
      n.mm = Ns$n.mm,
      n.main = Ns$n.main,
      n.hm = Ns$n.hm,
      n.mmblocks = Ns$n.mmblocks
    )

    # Create return --------------------------------------------------------------------- #

    out <- list(
      reg.table = reg.table,
      w = w,
      re.mm = re.mm,
      re.hm = re.hm,
      pred = pred,
      input = input,
      jags.out = if (isTRUE(monitor)) jags.out else NULL
    )

    class(out) <- "bml"

    return(out)
  } else {
    message("Data and model have been created without any errors.")
    invisible(list(
      modelstring = modelstring,
      jags.data = jags.data,
      jags.params = jags.params
    ))
  }
}
