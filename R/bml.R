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
#' lower-level units aggregate to higher levels (the micro–macro link).
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
#'   \item \strong{Level-2 predictors (X.L2):} Variables defined at level 2, separated by \code{+}.
#'   \item \strong{Level-3 predictors (X.L3):} Variables defined at level 3, separated by \code{+}.
#'   \item \strong{Multiple membership object (\code{mm()}):} Defines how lower-level (level-1) units are
#'         associated with higher-level (level-2) constructs using a user-specified weighting function
#'         (see details below).
#'   \item \strong{Hierarchical membership (\code{hm()}):} Specifies nesting of level-2 units within
#'         level-3 entities. Cross-classified structures can be modeled by including multiple \code{hm()} objects.
#' }
#'
#' \strong{Important:} The formula parser does not support \code{I()} inside \code{bml()}.
#' Create transformations and interactions in your data object before modeling.
#'
#' @section Multiple Membership Object \code{mm()}:
#' \preformatted{
#' mm(
#'   id  = id(l1id, l2id),
#'   mmc = mmc(X.L1),     # level-1 covariates to aggregate
#'   mmw = mmw(w ~ 1 / N, constraint = 1, ar = FALSE)
#' )
#' }
#'
#' \strong{Components:}
#' \itemize{
#'   \item \code{id(l1id, l2id)}: Specifies identifiers linking each level-1 unit (\code{l1id})
#'         to its corresponding level-2 entities (\code{l2id}).
#'   \item \code{mmc(X.L1)}: Specifies level-1 covariates aggregated across memberships.
#'         No intercept is allowed and interactions must be precomputed.
#'         If omitted, only the random-effects structure is modeled.
#'   \item \code{mmw(w ~ ..., constraint, ar)}: Defines the weight function (micro–macro link).
#' }
#'
#' \strong{Weight function details:}
#' \itemize{
#'   \item The right-hand side specifies how weights (\code{w}) are estimated. The function may be nonlinear
#'         and can include variables and parameters, provided it is identifiable and bounded.
#'   \item Must follow standard formula syntax.
#'   \item Parameters should be named \code{b1}, \code{b2}, \code{b3}, etc.
#'   \item Variables without explicit parameters are interpreted as \code{1 * X}, similar to an offset term.
#'   \item Each parameter must have a corresponding prior, and priors must be compatible with the specified
#'         weight function (see Tips below).
#'   \item \strong{Default:} If \code{mmw()} is omitted, equal weights (\code{w ~ 1 / N}) are assumed.
#'         Note that \code{N} is created internally.
#'   \item \strong{Generalized logistic example (Rosche, 2025):}
#'         \code{w ~ 1 / (1 + (N - 1) * exp(-(b1 * X1 + b2 * X2 + ...)))}.
#' }
#'
#' \strong{Constraint option:}
#' \itemize{
#'   \item \code{constraint = TRUE} (default): Weights for each level-2 entity sum to 1.
#'   \item \code{constraint = FALSE}: Weights are not constrained to sum to 1. This may lead to
#'         different scaling of aggregated covariates relative to their level-1 counterparts.
#' }
#'
#' \strong{Autoregressive random effects option:}
#' \itemize{
#'   \item \code{ar = TRUE}: Allows level-1 random effects to vary across level-2 entities as a random walk,
#'         with each effect centered on the corresponding entity’s effect from the previous instance.
#'   \item \code{ar = FALSE} (default): Assumes a single random effect per level-1 unit.
#' }
#'
#' @section Hierarchical Membership Object \code{hm()}:
#' \preformatted{
#' hm(id = l3id, name = l3name, type = RE, showFE = FALSE)
#' }
#'
#' \strong{Components:}
#' \itemize{
#'   \item \code{id = l3id}: Variable identifying level-3 groups.
#'   \item \code{name = l3name}: Optional labels for level-3 units.
#'   \item \code{type}: \code{"RE"} (default) or \code{"FE"}.
#'         \itemize{
#'           \item If \code{"FE"} is selected, each level-3 unit has its own intercept, and any level-3 predictors
#'                 are omitted. The first \code{l3id} serves as the reference category.
#'           \item If \code{"RE"} is selected, a normally distributed random-effects term is estimated.
#'         }
#'   \item \code{showFE}: If \code{TRUE} and \code{type = "FE"}, report the fixed effects; otherwise omit (default).
#' }
#'
#' @section Supported Families / Links:
#' \itemize{
#'   \item Gaussian (continuous): \code{family = "Gaussian"}
#'   \item Binomial (logistic): \code{family = "Binomial"}
#'   \item Weibull survival: \code{family = "Weibull"}, outcome: \code{Surv(time, event)}
#'   \item Cox survival: \code{family = "Cox"}, outcome: \code{Surv(time, event)}
#'   \item Not yet implemented: \code{family = "CondLogit"}
#' }
#'
#' @section Priors:
#' Priors can be specified for the following parameters: \code{b.l1}, \code{b.l2}, \code{b.l3}, \code{b.w},
#' \code{tau.l1}, \code{tau.l2}, and \code{tau.l3}.
#'
#' Supply a list of character strings, for example:
#' \preformatted{
#' priors = list("b.l1~dnorm(0,0.01)", "tau.l1~dscaled.gamma(25,1)")
#' }
#' For more details on priors in JAGS, see the
#' \href{https://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/jags_user_manual.pdf}{JAGS User Manual}.
#'
#' @section Tips:
#' \itemize{
#'   \item \strong{JAGS error on \code{w[i]} nodes:} The weight function must produce values consistent with the priors of all parameters.
#'         Negative or unbounded weights can push parameters outside their prior support, causing errors.
#'         Ensure weights sum to 1 and remain properly bounded. Standardizing variables (where appropriate) can help.
#'   \item \strong{Weight regressors:} Including regressors in the weight function increases model complexity and data demands.
#'         Start with moderately informative priors, for example: \code{priors = list("b.w" = "dnorm(0,0.1"))},
#'         and relax precision gradually as the model stabilizes.
#'   \item \strong{Input data:} The input data frame must represent level-1 units (one row per level-1 observation).
#'         All transformed and interaction variables should be precomputed before fitting the model.
#' }
#'
#' @param formula A symbolic model formula. The left-hand side specifies the outcome
#'   (e.g., \code{Y} or \code{survival::Surv(time, event)}), and the right-hand side includes fixed
#'   effects and membership structures, for example:
#'   \code{1 + X.L2 + X.L3 + mm(id(l1id, l2id), mmc(X.L1), mmw(w ~ 1/N)) + hm(id = l3id, type = RE)}.
#'   Use standard formula syntax; \code{I()} is not supported inside \code{bml()}—create transformed
#'   and interaction variables in the data before modeling.
#'
#' @param family A character string specifying the model family. Supported options are
#'   \code{"Gaussian"}, \code{"Binomial"}, \code{"Weibull"}, and \code{"Cox"}.
#'   (Not yet implemented: \code{"CondLogit"}.)
#'
#' @param priors A named list or character vector mapping parameter groups to JAGS priors.
#'   For example: \code{priors = list("b.l1" = "dnorm(0,0.01"))}.
#'   See “More details on priors” in the package documentation.
#'
#' @param inits A list of initial values used for all chains. If \code{NULL}, suitable
#'   initial values are chosen automatically by JAGS.
#'
#' @param n.iter Total number of MCMC iterations.
#' @param n.burnin Number of burn-in iterations to discard.
#' @param n.thin Thinning rate for MCMC sampling.
#' @param chains Number of MCMC chains to run.
#' @param seed Random seed for reproducibility.
#' @param run Logical; if \code{TRUE} (default), JAGS is executed to estimate the model.
#'
#' @param monitor Logical; if \code{TRUE}, additional components (weights, random effects,
#'   predictions, and raw JAGS output) are stored in the result object.
#'
#' @param modelfile Character or logical. If \code{TRUE}, the generated JAGS model is saved
#'   to \code{bml/temp/modelstring.txt}. If a file path is provided, \code{bml()} will
#'   use that file instead of generating new JAGS code and only build the data structure.
#'   Use \code{.libPaths()} to locate package storage directories.
#'
#' @param data A data frame where each row represents a level-1 observation (the unit of analysis).
#'
#' @return
#' A list containing the following elements:
#' \itemize{
#'   \item \code{reg.table} — Regression results summary.
#'   \item \code{w} — Estimated weights.
#'   \item \code{re.l1} — Level-1 random effects (if specified).
#'   \item \code{re.l3} — Level-3 random effects (if specified).
#'   \item \code{pred} — Posterior predictions (linear predictor for Gaussian models;
#'         survival time for Weibull models).
#'   \item \code{input} — Internally created model variables.
#'   \item \code{jags.out} — Raw, unformatted JAGS output.
#' }
#'
#' If \code{monitor = FALSE}, only \code{reg.table} is returned.
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' # Fit a Weibull survival model with multiple membership
#' m1 <- bml(
#'   Surv(govdur, earlyterm) ~
#'     1 +
#'     mm(id(pid, gid),
#'        mmc(fdep),
#'        mmw(w ~ 1/n)) +
#'     majority +
#'     hm(id = cid, type = RE),
#'   family  = "Weibull",
#'   monitor = TRUE,
#'   data    = coalgov
#' )
#'
#' # Inspect model outputs
#' m1$reg.table  # Regression summary
#' m1$w          # Estimated weights
#' m1$re.l1      # Level-1 random effects
#' m1$re.l3      # Level-3 random effects
#' m1$pred       # Posterior predictions
#' m1$input      # Internally generated variables
#' jags.out <- m1$jags.out  # Raw JAGS output
#'
#' # Summaries and visualization
#' summary(m1)
#' monetPlot(m1, "b.l1")  # Inspect posterior distribution of parameters
#' }
#'
#' @seealso
#' \code{\link[brms]{brm}}, \href{http://www.bristol.ac.uk/cmm/software/mlwin/}{MLwiN}
#'
#' @keywords Bayesian multilevel multiple-membership JAGS
#'
#' @export
#' @author Benjamin Rosche <benrosche@@nyu.edu>
#'
#' @references
#' Rosche, B. (2025). \emph{A Multilevel Model for Coalition Governments: Uncovering Dependencies
#' Within and Across Governments Due to Parties}.
#' \url{https://doi.org/10.31235/osf.io/4bafr}


bml <- function(formula, family="Gaussian", priors=NULL, inits=NULL, n.iter = 1000, n.burnin = 500, n.thin = max(1, floor((n.iter - n.burnin) / 1000)), chains=3, seed=NULL, run=T, parallel=F, monitor=T, modelfile=F, data=NULL) {

  # formula = sim.y ~ 1 + majority + mm(id(pid, gid), mmc(ipd), mmw(w ~ 1/n^exp(-(b0 + b1*rile.gov_SD)), c=T)); family = "Gaussian";  priors=c("b.w~dunif(0,1)", "b.l1~dnorm(0,1)", "tau.l2~dscaled.gamma(50,2)"); inits=NULL; n.iter=100; n.burnin=10; n.thin = max(1, floor((n.iter - n.burnin) / 1000)); chains = 3; seed = 123; run = T; parallel = F; monitor = T; modelfile = F; data = coalgov %>% rename(rile.gov_SD=hetero)
  # source("./R/dissectFormula.R"); source("./R/createData.R"); source("./R/editModelstring.R"); source("./R/createJagsVars.R"); source("./R/formatJags.R");

  # ---------------------------------------------------------------------------------------------- #
  # 0. Checks
  # ---------------------------------------------------------------------------------------------- #

  if(is.null(data)) stop("No data supplied.")

  # ---------------------------------------------------------------------------------------------- #
  # 1. Dissect formula
  # ---------------------------------------------------------------------------------------------- #

  DIR <- system.file(package = "bml")

  c(ids, vars, l1, l3) %<-% dissectFormula(formula, family, data) # updated (Jan 2025)

  # ---------------------------------------------------------------------------------------------- #
  # 2. Disentangle vars and data into l1 to l3
  # ---------------------------------------------------------------------------------------------- #

  c(data, level1, level2, level3, weight) %<-% createData(data, ids, vars, l1, l3) # updated (Feb 2025)

  # Remove varlist
  rm(vars)

  # ---------------------------------------------------------------------------------------------- #
  # 3. Create/edit jags modelstring
  # ---------------------------------------------------------------------------------------------- #

  modelstring <- editModelstring(family, priors, l1, l3, level1, level2, level3, weight, DIR, monitor, modelfile) # updated (Feb 2025)

  # Save or read modelstring
  if(isTRUE(modelfile)) {
    modelfile_path <- file.path(getwd(), "modelstring.txt")
    readr::write_file(modelstring, modelfile_path) # save model to file
    message("JAGS model saved to: ", modelfile_path)
  } else if(!isFALSE(modelfile) & length(modelfile)>0) {
    tryCatch({
      modelstring <- readr::read_file(modelfile) # read model from file
    }, error = function(e) {
      stop("Could not find/read model file in ", modelfile)
    })
  }

  # ---------------------------------------------------------------------------------------------- #
  # 4. Transform data into JAGS format
  # ---------------------------------------------------------------------------------------------- #

  c(ids, Ns, Xs, Ys, jags.params, jags.inits, jags.data) %<-% createJagsVars(data, family, level1, level2, level3, weight, ids, l1, l3, monitor, modelfile, chains, inits) # updated (Feb 2025)

  list2env(c(ids, Ns, Xs, Ys), envir=environment())

  # ---------------------------------------------------------------------------------------------- #
  # 5. Run JAGS
  # ---------------------------------------------------------------------------------------------- #

  if(run==T) {

    # Get seed
    if(is.null(seed)) seed <- round(runif(1, 0, 1000))

    if(parallel) {

      # Run parallel ----------------------------------------------------------------------------- #

      parallelfile <- tempfile(fileext = ".jags") # crate temp modelfile for parallel execution
      on.exit(unlink(parallelfile), add = TRUE) # ensure it will be deleted again after function call
      writeLines(modelstring, parallelfile) # save modelstring to parallelfile
      jags.out <- do.call(jags.parallel, list(data=jags.data, inits = jags.inits[1], n.chains = chains, parameters.to.save = jags.params, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, jags.seed = seed, model.file = parallelfile))

      # Three peculiarities about jags.parallel:
      # - It cannot read the model from textConnection(modelstring)
      # - It cannot read variables from the global environment - do.call needs to be used
      # - There seems to be a bug in that it wants just one list element of inits instead of n.chains number of list elements

    } else {

      # Run sequentially ------------------------------------------------------------------------- #

      set.seed(seed)
      jags.out <- jags(data=jags.data, inits = jags.inits, n.chains = chains, parameters.to.save = jags.params, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = textConnection(modelstring))

    }

    # Format JAGS output ------------------------------------------------------------------------- #

    c(reg.table, w, re.l1, re.l3, pred) %<-% formatJags(jags.out, monitor, Ns, l1, l3, level1, level2, level3, weight)

    # Prepare additional information ------------------------------------------------------------- #

    # Save info on input
    input <-
      c(
        list(
          "family"=family, "priors"=priors, "inits"=inits,
          "n.iter"=n.iter, "n.burnin"=n.burnin, "n.thin"=n.thin, "chains"=chains, "parallel"=parallel, "seed"=seed,
          "monitor"=monitor, "modelfile"=modelfile, "run"=run,
          "lhs" = level2$lhs, "l1vars"=level1$vars, "l2vars"=level2$vars, "l3vars"=level3$vars,
          "n.ul1"=Ns$n.ul1, "n.l1"=Ns$n.l1, "n.l2"=Ns$n.l2, "n.l3"=Ns$n.l3
        ),
        c(l1, l3)
      )

    # Create return ------------------------------------------------------------------------------ #

    out <- list("reg.table"=reg.table, "w"=w, "re.l1"=re.l1, "re.l3"=re.l3, "pred"=pred, "input"=input, "jags.out"=if(isTRUE(monitor)) jags.out else c())

    class(out) <- "bml"

    return(out)

  } else {

    message("Data and model have been created without any errors.")

  }

}
