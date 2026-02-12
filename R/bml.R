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
#' \code{\link[brms]{brms}} or \href{https://www.bristol.ac.uk/cmm/software/mlwin/}{MLwiN},
#' \strong{bml} lets users specify and estimate models in which membership weights are parameterized
#' through flexible formula syntax. This enables a more nuanced examination of how effects from
#' member-level units aggregate to group level (the micro-macro link).
#'
#' The package automatically generates JAGS code to fit the model and processes the output
#' to facilitate interpretation of model parameters and diagnostics.
#'
#' The package and modeling framework are introduced in:
#' Rosche, B. (2026). \emph{A Multilevel Model for Coalition Governments: Uncovering Party-Level
#' Dependencies Within and Between Governments}. \emph{Political Analysis}.
#'
#' For accessible introductions to multiple-membership models, see Fielding and Goldstein (2006)
#' and Beretvas (2010). Advanced treatments include Goldstein (2011, Ch. 13),
#' Rasbash and Browne (2001, 2008), Browne et al. (2001), and Leckie (2013).
#'
#' @section Formula Components:
#' \itemize{
#'   \item \strong{Outcome (Y):} The dependent variable. For survival models, use \code{Surv(time, event)}.
#'   \item \strong{Intercept:} Follows standard R formula conventions (like \code{lm()}):
#'     \itemize{
#'       \item \code{y ~ x}: Includes intercept by default
#'       \item \code{y ~ 1 + x}: Explicitly includes intercept (same as default)
#'       \item \code{y ~ 0 + x} or \code{y ~ -1 + x}: Excludes intercept
#'     }
#'   \item \strong{Main-level predictors (X.main):} Variables defined at the main (group) level, separated by \code{+}.
#'   \item \strong{HM-level predictors (X.hm):} Variables defined at the nesting level, separated by \code{+}.
#'   \item \strong{Multiple membership object (\code{mm()}):} Defines how member-level units are
#'         associated with group-level constructs using a user-specified weighting function.
#'         Multiple \code{mm()} objects can be specified with different weight functions.
#'   \item \strong{Hierarchical membership (\code{hm()}):} Specifies nesting of main-level units within
#'         higher-level entities. Cross-classified structures can be modeled by including multiple \code{hm()} objects.
#' }
#'
#' \strong{Formula Features:} The main formula and \code{vars()} specifications support standard R formula syntax:
#' \itemize{
#'   \item \strong{Interactions:} Use \code{*} for main effects plus interaction, or \code{:} for interaction only.
#'         Example: \code{y ~ a * b} expands to \code{y ~ a + b + a:b}.
#'   \item \strong{Transformations:} Use \code{I()} for arithmetic operations.
#'         Example: \code{y ~ I(x^2)} or \code{y ~ I(a + b)}.
#' }
#'
#' These features work in:
#' \itemize{
#'   \item \strong{Main formula:} \code{y ~ 1 + a * b + I(x^2)}
#'   \item \strong{mm() vars:} \code{vars(a * b)} or \code{vars(I(x^2))}
#'   \item \strong{hm() vars:} \code{vars(a:b)} or \code{vars(I(log(x)))}
#' }
#'
#' \strong{Note on weight functions:} The \code{fn()} weight function in \code{mm()} does NOT support
#' interactions or \code{I()} transformations. Users must pre-create any needed transformed variables
#' in their data before using them in weight functions. For example, instead of \code{fn(w ~ b1 * x^2)},
#' first create \code{data$x_sq <- data$x^2} and use \code{fn(w ~ b1 * x_sq)}.
#'
#' \strong{Note on intercepts:} Intercept syntax (\code{1}, \code{0}, \code{-1}) only applies to the main formula.
#' Numeric literals in \code{vars()} are ignored (e.g., \code{vars(1 + x)} is equivalent to \code{vars(x)}).
#'
#' @section Multiple Membership Object \code{mm()}:
#' \preformatted{
#' mm(
#'   id   = id(mmid, mainid),
#'   vars = vars(X.mm),
#'   fn   = fn(w ~ 1/n, c = TRUE),
#'   RE   = TRUE,
#'   ar   = FALSE
#' )
#' }
#'
#' \strong{Components:}
#' \itemize{
#'   \item \code{id(mmid, mainid)}: Specifies identifiers linking each member-level unit (\code{mmid})
#'         to its corresponding group-level entities (\code{mainid}).
#'   \item \code{vars(X.mm)}: Specifies member-level covariates aggregated across memberships.
#'         Use \code{+} to include multiple variables. Supports interactions (\code{*}, \code{:})
#'         and transformations (\code{I()}). Set to \code{NULL} for RE-only blocks.
#'   \item \code{fn(w ~ ..., c)}: Defines the weight function (micro-macro link).
#'         The \code{c} parameter controls weight normalization: when \code{c = TRUE} (default),
#'         weights are normalized to sum to 1 within each group
#'         (\eqn{\tilde{w}_{ik} = w_{ik} / \sum_{k} w_{ik}}).
#'         Set \code{c = FALSE} for unnormalized weights (e.g., when aggregating sums).
#'         Note: Does not support interactions or \code{I()} - pre-create transformed variables.
#'   \item \code{RE}: Logical; if \code{TRUE}, include random effects for this block.
#'         Automatically \code{TRUE} if \code{vars = NULL}.
#'   \item \code{ar}: Logical; if \code{TRUE}, member-level random effects evolve as a
#'         random walk across repeated participations in groups. This captures dynamics
#'         where a member's unobserved heterogeneity changes over time. Default: \code{FALSE}.
#' }
#'
#' \strong{Multiple mm() blocks:}
#' You can specify multiple \code{mm()} blocks with different weight functions.
#' However, \code{RE = TRUE} can only be specified for one \code{mm()} block.
#' \preformatted{
#' mm(id = id(pid, gid), vars = vars(X.mm.1), fn = fn(w ~ 1/n), RE = FALSE) +
#' mm(id = id(pid, gid), vars = vars(X.mm.2), fn = fn(w ~ pseat == max(pseat)), RE = FALSE) +
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
#'         Supports interactions (\code{*}, \code{:}) and transformations (\code{I()}).
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
#'   "tau.mm.1 ~ dscaled.gamma(25, 1)"
#' )
#' }
#'
#' @param formula A symbolic model formula. See 'Formula Components' section for details.
#'   The general structure is: \code{outcome ~ 1 + predictors + mm(...) + hm(...)}.
#'   For survival models, use \code{Surv(time, event)} on the left-hand side.
#'
#' @param family Character string specifying the outcome distribution and link function.
#'   Options:
#'   \itemize{
#'     \item \code{"Gaussian"}: Normal distribution with identity link (continuous outcomes)
#'     \item \code{"Binomial"}: Binomial distribution with logit link (binary outcomes)
#'     \item \code{"Weibull"}: Weibull survival model (requires \code{Surv(time, event)} outcome)
#'     \item \code{"Cox"}: Cox proportional hazards model (requires \code{Surv(time, event)} outcome)
#'   }
#'
#' @param priors Named list or character vector of JAGS prior specifications.
#'   Parameter names follow JAGS naming conventions:
#'   \itemize{
#'     \item \strong{Main level:} \code{b[x]} for main-equation coefficients
#'       (e.g., \code{"b[1] ~ dnorm(0, 0.01)"} for the intercept)
#'     \item \strong{HM level:} \code{b.hm.k[x]} for hm block \code{k} coefficients,
#'       \code{tau.hm.k} for hm block \code{k} random effect precision
#'     \item \strong{MM level:} \code{b.mm.k[x]} for mm block \code{k} coefficients,
#'       \code{b.w.k[x]} for weight function parameters in mm block \code{k},
#'       \code{tau.mm.g} for mm random effect precision (indexed by member-group
#'       ID group \code{g})
#'     \item \strong{Other:} \code{shape} (Weibull shape parameter),
#'       \code{lambda0[k]} (Cox baseline hazard intervals)
#'   }
#'   \strong{Note:} Priors on variance components must be specified on the
#'   \emph{precision} scale (\code{tau = 1/sigma^2}), not the standard deviation,
#'   since JAGS parameterizes normal distributions using precision.
#'   Example: \code{list("b.mm.1 ~ dnorm(0, 0.01)", "tau.mm.1 ~ dgamma(2, 0.1)")}.
#'   Default priors are weakly informative.
#'
#' @param inits List of initial values for MCMC chains. Applied to all chains. If
#'   \code{NULL}, JAGS generates initial values automatically. Weight function
#'   parameters (\code{b.w.k}) are always initialized at 0 by default to prevent
#'   numerical instability (e.g., \code{ilogit} with extreme inputs). User-supplied
#'   inits override these defaults.
#'
#' @param n.iter Total number of MCMC iterations per chain. Default: 1000.
#'   Increase for better convergence (e.g., 10000-50000 for production models).
#'
#' @param n.burnin Number of burn-in iterations to discard at the start of each
#'   chain. Default: 500. Should be sufficient for chains to reach stationarity.
#'
#' @param n.thin Thinning rate: save every k-th iteration to reduce autocorrelation.
#'   Default: \code{max(1, floor((n.iter - n.burnin) / 1000))} (targets ~1000 samples).
#'   Increase if posterior samples show high autocorrelation.
#'
#' @param n.chains Number of MCMC chains. Default: 3. Use 3-4 chains to assess
#'   convergence via Gelman-Rubin diagnostics.
#'
#' @param seed Integer random seed for reproducibility. If \code{NULL}, results
#'   will vary across runs.
#'
#' @param run Logical; if \code{TRUE} (default), JAGS is executed and the model
#'   is fitted. If \code{FALSE}, returns the model specification without fitting
#'   (useful for inspecting generated JAGS code or data structures).
#'
#' @param parallel Logical; if \code{TRUE}, run MCMC chains in parallel using
#'   multiple cores. Requires parallel backend setup. Default: \code{FALSE}.
#'
#' @param monitor Logical; if \code{TRUE}, store full MCMC chains and additional
#'   outputs for diagnostic plots. Required for \code{\link{monetPlot}} and
#'   \code{\link{mcmcDiag}}. Default: \code{TRUE}.
#'
#' @param modelfile Logical or character path:
#'   \itemize{
#'     \item \code{FALSE} (default): JAGS code generated internally
#'     \item \code{TRUE}: Save generated JAGS code to \code{modelstring.txt} in working directory
#'     \item Character path: Read JAGS code from specified file instead of generating
#'   }
#'
#' @param cox_intervals For Cox models only: controls baseline hazard flexibility
#'   and computational efficiency.
#'   \itemize{
#'     \item \code{NULL} (default): Non-parametric baseline hazard using all unique
#'       event times (maximum flexibility, slower for large datasets)
#'     \item Integer k: Piecewise constant baseline hazard with k intervals
#'       (faster, suitable for datasets with many unique event times). Recommended:
#'       k = 10-20 for most applications.
#'   }
#'
#' @param data Data frame in member-level (long) format where each row represents
#'   a member-level observation. Must contain all variables referenced in the
#'   formula, including identifiers specified in \code{id()}.
#'
#' @return A list of class \code{"bml"} containing:
#'   \itemize{
#'     \item \code{reg.table}: Data frame of posterior summaries with columns
#'       \code{Parameter}, \code{mean}, \code{sd}, \code{lb}, \code{ub}
#'       (95\% credible interval bounds).
#'     \item \code{w}: List of weight matrices (one per \code{mm()} block). Each
#'       matrix has rows = groups and columns = members within each group.
#'     \item \code{re.mm}: List of member-level random effects (one per mmid group
#'       with \code{RE = TRUE}). Vector for standard RE, matrix for autoregressive RE.
#'     \item \code{re.hm}: List of nesting-level random effects (one per \code{hm()}
#'       block with \code{type = "RE"}).
#'     \item \code{pred}: Vector of predicted values (posterior means) for each group.
#'     \item \code{input}: List of model specifications including \code{family},
#'       \code{mm} and \code{hm} block info, sample sizes, and MCMC settings.
#'     \item \code{jags.out}: Full R2jags output object (if \code{monitor = TRUE};
#'       \code{NULL} otherwise). Contains posterior samples, MCMC chains, and
#'       convergence diagnostics.
#'   }
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' # Basic multiple-membership model
#' # Parties (pid) within governments (gid), nested in countries (cid)
#' m1 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(
#'       id   = id(pid, gid),
#'       vars = vars(cohesion),
#'       fn   = fn(w ~ 1/n, c = TRUE),
#'       RE   = TRUE
#'     ) +
#'     hm(id = id(cid), type = "RE"),
#'   family = "Weibull",
#'   data   = coalgov
#' )
#'
#' # View results
#' summary(m1)
#' monetPlot(m1, "b[2]")  # Plot for majority coefficient
#'
#' # Multiple mm() blocks with different weight functions
#' m2 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion),
#'        fn = fn(w ~ cohesion == max(cohesion)), RE = FALSE) +
#'     mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE),
#'   family = "Weibull",
#'   data   = coalgov
#' )
#'
#' # Cox model with piecewise baseline hazard (faster for large datasets)
#' m3 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE),
#'   family = "Cox",
#'   cox_intervals = 10,  # Use 10 intervals instead of all unique times
#'   data   = coalgov
#' )
#'
#' # Parameterized weight function
#' # ilogit() bounds raw weights between 0 and 1; c = TRUE normalizes to sum to 1
#' m4 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(
#'       id   = id(pid, gid),
#'       vars = vars(cohesion),
#'       fn   = fn(w ~ ilogit(b0 + b1 * pseat), c = TRUE),
#'       RE   = FALSE
#'     ),
#'   family = "Weibull",
#'   data   = coalgov
#' )
#'
#' # Fixed coefficients (offsets)
#' m5 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + fix(majority, 1) + # Fix majority coefficient to 1.0
#'     mm(
#'       id   = id(pid, gid),
#'       vars = vars(rile),
#'       fn   = fn(w ~ 1/n, c = TRUE),
#'       RE   = FALSE
#'     ),
#'   family = "Weibull",
#'   data   = coalgov
#' )
#'
#' # Custom priors
#' m6 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE),
#'   family = "Weibull",
#'   priors = list(
#'     "b[1] ~ dnorm(0, 0.01)",       # Intercept prior
#'     "b.mm.1 ~ dnorm(0, 0.1)",      # MM coefficient prior
#'     "tau.mm.1 ~ dgamma(2, 0.5)"    # MM precision prior
#'   ),
#'   data   = coalgov
#' )
#'
#' # Cross-classified model (multiple hm blocks)
#' # Governments are cross-classified by country and election year
#' m7 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE) +
#'     hm(id = id(cid), type = "RE") +
#'     hm(id = id(year), type = "RE"),
#'   family = "Weibull",
#'   data   = coalgov |> mutate(year = format(election, "%Y") |> as.integer())
#' )
#' }
#'
#' @seealso
#' \code{\link{summary.bml}} for model summaries,
#' \code{\link{monetPlot}} for posterior visualization,
#' \code{\link{mcmcDiag}} for convergence diagnostics,
#' \code{\link{mm}}, \code{\link{hm}} for model specification helpers
#'
#' @references
#' Rosche, B. (2026). A Multilevel Model for Coalition Governments: Uncovering
#' Party-Level Dependencies Within and Between Governments. \emph{Political Analysis}.
#'
#' Browne, W. J., Goldstein, H., & Rasbash, J. (2001). Multiple membership
#' multiple classification (MMMC) models. \emph{Statistical Modelling}, 1(2), 103-124.
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
  n.chains = 3,
  seed = NULL,
  run = TRUE,
  parallel = FALSE,
  monitor = TRUE,
  modelfile = FALSE,
  cox_intervals = NULL,
  data = NULL
) {

  # Test call:
  # formula <- sim.y ~ 1 + majority + mm( id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1 / n, c = TRUE), RE = TRUE ) + hm(id = id(cid), vars = NULL, type = "RE"); family = "Gaussian"; run=F; data = coalgov

  # ========================================================================================== #
  # 0. Checks
  # ========================================================================================== #

  if (is.null(data)) {
    stop("No data supplied.")
  }
  data <- dplyr::ungroup(data)

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
  # 1b. Validate data
  # ========================================================================================== #

  # Collect all RHS variables that need to be checked for missing values
  vars_to_check <- c()

  # Main-level variables - use all.vars() to get base variable names
  # This handles interactions (a*b, a:b) and I() transformations
  if (!is.null(formula_parts$main_formula)) {
    main_base_vars <- all.vars(formula_parts$main_formula)
    vars_to_check <- c(vars_to_check, main_base_vars)
  }

  # Fixed main-level variables
  if (!is.null(formula_parts$mainvars_fixed)) {
    fixed_vars <- sapply(formula_parts$mainvars_fixed, function(x) x$var)
    fixed_vars <- fixed_vars[fixed_vars != "X0"]
    vars_to_check <- c(vars_to_check, fixed_vars)
  }

  # MM-level variables and weight function variables
  weight_vars <- c()
  if (has_mm) {
    for (m in mm) {
      # mm vars (free)
      if (!is.null(m$vars)) {
        if (is.list(m$vars) && !is.null(m$vars$free)) {
          vars_to_check <- c(vars_to_check, m$vars$free)
          # Fixed vars
          if (!is.null(m$vars$fixed)) {
            vars_to_check <- c(vars_to_check, sapply(m$vars$fixed, function(x) x$var))
          }
        } else {
          vars_to_check <- c(vars_to_check, as.character(m$vars))
        }
      }
      # Weight function variables
      if (!is.null(m$fn$vars) && length(m$fn$vars) > 0) {
        vars_to_check <- c(vars_to_check, m$fn$vars)
        weight_vars <- c(weight_vars, m$fn$vars)
      }
    }
  }

  # HM-level variables
  if (has_hm) {
    for (h in hm) {
      if (!is.null(h$vars)) {
        if (is.list(h$vars) && !is.null(h$vars$free)) {
          vars_to_check <- c(vars_to_check, h$vars$free)
          if (!is.null(h$vars$fixed)) {
            vars_to_check <- c(vars_to_check, sapply(h$vars$fixed, function(x) x$var))
          }
        } else {
          vars_to_check <- c(vars_to_check, as.character(h$vars))
        }
      }
    }
  }

  # Check for missing values in RHS variables
  vars_to_check <- unique(vars_to_check)
  for (v in vars_to_check) {
    if (v %in% names(data) && any(is.na(data[[v]]))) {
      n_missing <- sum(is.na(data[[v]]))
      stop("Missing values detected in variable '", v, "' (", n_missing, " observations).\n",
           "Please remove or impute missing values before fitting the model.")
    }
  }

  # Check if weight function variables are constant within mainid groups
  # Only check if parameters are being estimated (not for simple aggregation like w ~ 1/n)
  if (has_mm && length(weight_vars) > 0) {
    # Check if any mm block has parameters to estimate
    has_params <- any(sapply(mm, function(m) length(m$fn$params) > 0))

    if (has_params) {
      mainid_var <- mm[[1]]$id[2]  # mainid is the second element of id()
      weight_vars <- unique(weight_vars)

      for (wv in weight_vars) {
        if (wv %in% names(data)) {
          # Check variance within each mainid group
          var_by_group <- data %>%
            dplyr::group_by(.data[[mainid_var]]) %>%
            dplyr::summarise(var = stats::var(.data[[wv]], na.rm = TRUE), .groups = "drop")

          # If all variances are 0 (or NA for single-member groups), the variable is constant
          non_na_vars <- var_by_group$var[!is.na(var_by_group$var)]
          if (length(non_na_vars) > 0 && all(non_na_vars == 0)) {
            warning("Weight function variable '", wv, "' is constant across members within groups.")
          }
        }
      }
    }
  }

  # Check for duplicate member-group combinations
  if (has_mm) {
    mmid_var <- mm[[1]]$id[1]    # mmid is the first element of id()
    mainid_var <- mm[[1]]$id[2]  # mainid is the second element of id()

    duplicates <- data %>%
      dplyr::group_by(.data[[mmid_var]], .data[[mainid_var]]) %>%
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::ungroup()

    if (nrow(duplicates) > 0) {
      n_dups <- nrow(duplicates)
      stop("Duplicate member-group combinations detected (", n_dups, " rows). ",
           "Each member (", mmid_var, ") should appear only once per group (", mainid_var, ").")
    }
  }

  # ========================================================================================== #
  # 2. Create data structures
  # ========================================================================================== #

  data_parts <- createData(data, formula_parts)

  data      <- data_parts$data
  mm_blocks <- data_parts$mm_blocks
  main      <- data_parts$main
  hm_blocks <- data_parts$hm_blocks

  # ========================================================================================== #
  # 3. Create JAGS modelstring
  # ========================================================================================== #

  modelstring <- createModelstring(
    family,
    priors,
    mm_blocks,
    main,
    hm_blocks,
    mm,
    hm,
    DIR,
    monitor,
    modelfile,
    cox_intervals
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
    n.chains,
    inits,
    cox_intervals
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
          n.chains = n.chains,
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
        n.chains = n.chains,
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
      hm,
      family,
      cox_intervals
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
        RE = mm_blocks[[k]]$RE,
        mmid_group = mm_blocks[[k]]$mmid_group
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
      n.chains = n.chains,
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
