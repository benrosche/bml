#' @title Bayesian Multilevel Models for Micro-Macro Analysis (Multiple Membership) Using JAGS
#'
#' @description
#' The \strong{bml} package estimates micro-to-macro regressions: Bayesian
#' multiple-membership multilevel models in which the aggregation from member-level
#' records to group-level outcomes is an explicit, estimable object — weights
#' (\code{\link{w}}), aggregation functions (\code{\link{fn}}: mean, variance,
#' concentration, thresholds, smooth max/min, CES means, covariance, or a
#' user-written reduction), member/unit random and fixed effects
#' (\code{\link{re}}, \code{\link{fe}}), and cross-level or feature-by-feature
#' interactions via named blocks.
#'
#' JAGS must be installed separately: \url{https://sourceforge.net/projects/mcmc-jags/}.
#'
#' @details
#' The general formula structure is
#' \preformatted{
#' outcome ~ 1 + predictors + mm(...) + hm(...)
#' }
#' where \code{\link{mm}} defines a multiple-membership block (the micro-macro
#' link) and \code{\link{hm}} a hierarchical nesting level. Multiple \code{mm()}
#' blocks stack features (e.g. mean + variance + concentration); multiple
#' \code{hm()} blocks give cross-classification. For survival families use
#' \code{Surv(time, event)} on the left-hand side.
#'
#' Named \code{mm()} blocks can be referenced in main-formula interactions:
#' \preformatted{
#' Y ~ education + Ax:education +
#'   mm(name = Ax, id = id(task, occ), vars = vars(x), w = w(~ importance), fn = fn("sum"))
#' }
#' Both cross-level (feature x macro variable) and block x block (feature x
#' feature) interactions are supported; the macro variable must also appear as
#' a main effect, and a bare (non-interacted) block reference is an error.
#'
#' The package is introduced in Rosche (2026), \emph{Political Analysis}.
#'
#' @param formula A symbolic model formula; see Details.
#' @param data Data frame in member-level (long) format: one row per membership.
#'   Must contain all variables referenced in the formula, including the
#'   identifiers in \code{id()}.
#' @param family Model family: a family function (\code{gaussian()},
#'   \code{\link{bernoulli}()}, \code{\link{weibull}()}, \code{\link{cox}()}) or
#'   a string (\code{"gaussian"}, \code{"bernoulli"}, \code{"weibull"},
#'   \code{"cox"}). Cox baseline-hazard intervals travel with the family:
#'   \code{cox(intervals = 10)}.
#' @param prior Prior specifications built with \code{\link{prior}} and combined
#'   with \code{c()}; see \code{\link{get_prior}} for the settable parameters.
#'   Raw JAGS strings (e.g. \code{"b.fn.1[1] ~ dnorm(0, 0.1)"}) pass through
#'   untranslated as an escape hatch.
#' @param inits List of initial values for MCMC chains (applied to all chains).
#'   Weight and \code{fn} shape parameters are started from dispersed values
#'   drawn separately for each chain, as R-hat assumes; user values override and
#'   are applied as given.
#' @param iter Total number of MCMC iterations per chain. Default: 1000.
#' @param warmup Number of warmup (burn-in) iterations discarded from the start
#'   of each chain. Default: 500.
#' @param thin Thinning rate. Default targets ~1000 retained draws.
#' @param chains Number of MCMC chains. Default: 3.
#' @param cores Number of cores; \code{cores > 1} runs chains in parallel.
#'   Default: 1.
#' @param n.adapt Number of JAGS adaptation iterations used to tune the samplers
#'   before warmup. Tuning is frozen once adaptation ends, so models with
#'   non-conjugate nodes (estimated \code{w()} or \code{fn()} parameters) benefit
#'   from a longer adaptation phase. Default: 1000.
#' @param seed Integer random seed for reproducibility.
#' @param run Logical; if \code{FALSE}, returns the generated JAGS model string,
#'   data, and monitors without fitting (see also \code{\link{make_jagscode}}).
#' @param monitor Character vector selecting retained posterior capabilities.
#'   The default, \code{"summary"}, retains posterior summaries and convergence
#'   diagnostics but no draws. Add \code{"parameters"},
#'   \code{"random_effects"}, \code{"weights"}, \code{"fitted"},
#'   \code{"predictive"}, or \code{"log_lik"} as needed. \code{"full"} stores
#'   all family-supported draws once in a compact
#'   \code{posterior::draws_array}; \code{"jags"} explicitly retains the raw
#'   R2jags object.
#' @param modelfile Logical or character path: \code{TRUE} saves the generated
#'   JAGS code to \code{modelstring.txt}; a path reads JAGS code from that file
#'   instead of generating it.
#' @param ... Not used; catches removed legacy arguments (\code{n.iter},
#'   \code{n.burnin}, \code{n.thin}, \code{n.chains}, \code{parallel},
#'   \code{priors}, \code{cox_intervals}) with a migration message.
#'
#' @return A list of class \code{"bml"} with elements \code{reg.table}
#'   (posterior summaries; labeled terms), \code{w} (weight matrices per block),
#'   \code{re.mm}/\code{re.hm} (random effects), \code{pred} (posterior
#'   predictive means), \code{fitted} (posterior means of the linear predictor),
#'   \code{input} (model metadata), \code{diagnostics}, \code{storage},
#'   \code{draws} (requested compact posterior draws), and \code{jags.out}
#'   (non-\code{NULL} only when \code{monitor = "jags"}).
#'
#' @examples
#' \dontrun{
#' data(coalgov)
#'
#' # Additive aggregation with member random effects, country nesting
#' m1 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("sum"), RE = TRUE) +
#'     hm(id = id(cid)),
#'   data = coalgov,
#'   family = weibull()
#' )
#' summary(m1)
#'
#' # Estimated weights (b0, b1 are free parameters by the one-parameter rule)
#' m2 <- bml(
#'   Surv(dur_wkb, event_wkb) ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion),
#'        w = w(~ ilogit(b0 + b1 * pseat), scale = TRUE),
#'        fn = fn("sum")),
#'   data = coalgov,
#'   family = weibull()
#' )
#'
#' # Emergent features: variance + concentration, stacked
#' m3 <- bml(
#'   sim.y ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("var")) +
#'     mm(id = id(pid, gid), w = w(~ pseat, scale = TRUE), fn = fn("hhi")),
#'   data = coalgov,
#'   family = gaussian()
#' )
#'
#' # Structured priors
#' m4 <- bml(
#'   sim.y ~ 1 + majority +
#'     mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("sum"), RE = TRUE),
#'   data = coalgov,
#'   family = gaussian(),
#'   prior = c(
#'     prior(normal(0, 5), class = "b"),
#'     prior(exponential(1), class = "sd")
#'   )
#' )
#' }
#'
#' @seealso
#' \code{\link{mm}}, \code{\link{hm}}, \code{\link{w}}, \code{\link{fn}},
#' \code{\link{re}}, \code{\link{fe}}, \code{\link{prior}}, \code{\link{get_prior}},
#' \code{\link{summary.bml}}, \code{\link{fixef.bml}}, \code{\link{ranef.bml}},
#' \code{\link{loo.bml}}, \code{\link{varDecomp}}, \code{\link{bmlCompare}},
#' \code{\link{make_jagscode}}
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
  data,
  family = stats::gaussian(),
  prior = NULL,
  inits = NULL,
  iter = 1000,
  warmup = 500,
  thin = max(1, floor((iter - warmup) / 1000)),
  chains = 3,
  cores = 1,
  n.adapt = 1000,
  seed = NULL,
  run = TRUE,
  monitor = "summary",
  modelfile = FALSE,
  ...
) {

  # ========================================================================================== #
  # 0. Checks
  # ========================================================================================== #

  dots <- list(...)
  legacy <- intersect(names(dots),
                      c("n.iter", "n.burnin", "n.thin", "n.chains", "parallel",
                        "priors", "cox_intervals"))
  if (length(legacy) > 0) {
    stop("Removed argument(s): ", paste(legacy, collapse = ", "), ". The bml API follows ",
         "brms now:\n",
         "  n.iter -> iter, n.burnin -> warmup, n.thin -> thin, n.chains -> chains,\n",
         "  parallel -> cores (integer), priors -> prior (see ?prior),\n",
         "  cox_intervals -> family = cox(intervals = ...)", call. = FALSE)
  }
  if (length(dots) > 0) {
    stop("Unknown argument(s): ", paste(names(dots), collapse = ", "), call. = FALSE)
  }

  if (missing(data) || is.null(data)) {
    stop("No data supplied.", call. = FALSE)
  }
  data <- dplyr::ungroup(data)

  fam <- validate_family(family)
  family_str <- fam$family
  cox_intervals <- fam$intervals
  monitor_spec <- normalize_monitor(monitor, family_str)

  # ========================================================================================== #
  # 1. Dissect formula
  # ========================================================================================== #

  formula_parts <- dissectFormula(formula, family_str, data)

  mm_spec <- formula_parts$mm
  hm_spec <- formula_parts$hm

  has_mm <- length(mm_spec) > 0
  has_hm <- length(hm_spec) > 0

  # ========================================================================================== #
  # 1b. Validate data
  # ========================================================================================== #

  vars_to_check <- c()

  if (!is.null(formula_parts$main_formula)) {
    vars_to_check <- c(vars_to_check, all.vars(formula_parts$main_formula))
  }

  if (!is.null(formula_parts$mainvars_fixed)) {
    fixed_vars <- sapply(formula_parts$mainvars_fixed, function(x) x$var)
    vars_to_check <- c(vars_to_check, fixed_vars[fixed_vars != "X0"])
  }

  weight_vars <- c()
  if (has_mm) {
    for (m in mm_spec) {
      if (!is.null(m$vars)) {
        vars_to_check <- c(vars_to_check, m$vars$free)
        if (!is.null(m$vars$fixed)) {
          vars_to_check <- c(vars_to_check, sapply(m$vars$fixed, function(x) x$var))
        }
      }
      wv <- intersect(m$w$vars, names(data))
      vars_to_check <- c(vars_to_check, wv, m$w$agg_vars)
      weight_vars <- c(weight_vars, wv, m$w$agg_vars)
      if (!is.null(m$RE)) vars_to_check <- c(vars_to_check, m$RE$slopes)
      if (!is.null(m$FE)) vars_to_check <- c(vars_to_check, m$FE$slopes)
    }
  }

  if (has_hm) {
    for (h in hm_spec) {
      if (!is.null(h$RE)) vars_to_check <- c(vars_to_check, h$RE$slopes)
      if (!is.null(h$FE)) vars_to_check <- c(vars_to_check, h$FE$slopes)
    }
  }

  vars_to_check <- unique(vars_to_check)
  for (v in vars_to_check) {
    if (v %in% names(data) && any(is.na(data[[v]]))) {
      n_missing <- sum(is.na(data[[v]]))
      stop("Missing values detected in variable '", v, "' (", n_missing, " observations).\n",
           "Please remove or impute missing values before fitting the model.", call. = FALSE)
    }
  }

  # Weight covariates that are constant within groups cannot inform estimated weights
  if (has_mm && length(weight_vars) > 0) {
    has_params <- any(sapply(mm_spec, function(m) length(m$w$params) > 0))
    if (has_params) {
      mainid_var <- mm_spec[[1]]$id[2]
      for (wv in unique(weight_vars)) {
        if (wv %in% names(data)) {
          var_by_group <- data %>%
            dplyr::group_by(.data[[mainid_var]]) %>%
            dplyr::summarise(var = stats::var(.data[[wv]], na.rm = TRUE), .groups = "drop")
          non_na_vars <- var_by_group$var[!is.na(var_by_group$var)]
          if (length(non_na_vars) > 0 && all(non_na_vars == 0)) {
            warning("Weight variable '", wv, "' is constant across members within groups.",
                    call. = FALSE)
          }
        }
      }
    }
  }

  # Duplicate member-group combinations
  if (has_mm) {
    for (m in mm_spec) {
      mmid_var <- m$id[1]
      mainid_var <- m$id[2]
      dat_m <- data[!is.na(data[[mmid_var]]), , drop = FALSE]
      dup_n <- dat_m %>%
        dplyr::count(.data[[mmid_var]], .data[[mainid_var]]) %>%
        dplyr::filter(n > 1) %>%
        nrow()
      if (dup_n > 0) {
        stop("Duplicate member-group combinations detected. Each member (", mmid_var,
             ") should appear only once per group (", mainid_var, ").", call. = FALSE)
      }
    }
  }

  # ========================================================================================== #
  # 2. Resolve priors
  # ========================================================================================== #

  prior_table <- build_prior_table(formula_parts, family_str)
  prior_spec <- resolve_priors(prior, prior_table)

  # ========================================================================================== #
  # 3. Create data structures
  # ========================================================================================== #

  data_parts <- createData(data, formula_parts)

  data      <- data_parts$data
  mm_blocks <- data_parts$mm_blocks
  main      <- data_parts$main
  hm_blocks <- data_parts$hm_blocks
  interactions <- data_parts$interactions

  # ========================================================================================== #
  # 4. Create JAGS modelstring
  # ========================================================================================== #

  modelstring <- createModelstring(
    family_str,
    prior_spec,
    mm_blocks,
    main,
    hm_blocks,
    interactions,
    monitor_spec,
    cox_intervals
  )

  if (isTRUE(modelfile)) {
    modelfile_path <- file.path(getwd(), "modelstring.txt")
    readr::write_file(modelstring, modelfile_path)
    message("JAGS model saved to: ", modelfile_path)
  } else if (!isFALSE(modelfile) && length(modelfile) > 0 && is.character(modelfile)) {
    tryCatch(
      {
        modelstring <- readr::read_file(modelfile)
      },
      error = function(e) {
        stop("Could not find/read model file in ", modelfile, call. = FALSE)
      }
    )
  }

  # ========================================================================================== #
  # 5. Transform data into JAGS format
  # ========================================================================================== #

  jags_vars <- createJagsVars(
    data,
    family_str,
    mm_blocks,
    main,
    hm_blocks,
    interactions,
    monitor_spec,
    chains,
    inits,
    cox_intervals
  )

  Ns <- jags_vars$Ns
  Xs <- jags_vars$Xs
  Ys <- jags_vars$Ys
  jags.params <- jags_vars$jags.params
  jags.inits <- jags_vars$jags.inits
  jags.data <- jags_vars$jags.data

  # ========================================================================================== #
  # 6. Run JAGS
  # ========================================================================================== #

  if (run) {
    if (is.null(seed)) {
      seed <- round(runif(1, 0, 1000))
    }

    if (cores > 1) {
      # R2jags::jags.parallel is deliberately not used: it forwards no n.adapt, so adaptation
      # would silently stay at R2jags's default of 100 iterations however n.adapt is set here.
      # Each worker runs one chain and calls the inits function itself, so every chain also
      # gets its own dispersed starting values. The chain outputs are recombined exactly the
      # way jags.parallel does it.
      parallelfile <- tempfile(fileext = ".jags")
      on.exit(unlink(parallelfile), add = TRUE)
      writeLines(modelstring, parallelfile)

      cluster <- parallel::makeCluster(min(cores, chains), methods = FALSE)
      on.exit(parallel::stopCluster(cluster), add = TRUE)

      chain_fits <- parallel::parLapply(
        cl = cluster,
        X = seed + seq_len(chains),
        fun = function(chain_seed, jags.data, jags.inits, jags.params, model.file,
                       iter, warmup, thin, n.adapt) {
          set.seed(chain_seed)
          R2jags::jags(
            data = jags.data,
            inits = jags.inits,
            parameters.to.save = jags.params,
            model.file = model.file,
            n.chains = 1,
            n.iter = iter,
            n.burnin = warmup,
            n.thin = thin,
            n.adapt = n.adapt,
            jags.seed = chain_seed,
            progress.bar = "none"
          )
        },
        jags.data = jags.data,
        jags.inits = jags.inits,
        jags.params = jags.params,
        model.file = parallelfile,
        iter = iter,
        warmup = warmup,
        thin = thin,
        n.adapt = n.adapt
      )

      template <- chain_fits[[1]]$BUGSoutput$sims.array
      sims.array <- array(
        NA_real_,
        dim = c(dim(template)[1], chains, dim(template)[3]),
        dimnames = list(NULL, NULL, dimnames(template)[[3]])
      )
      for (ch in seq_len(chains)) {
        sims.array[, ch, ] <- chain_fits[[ch]]$BUGSoutput$sims.array[, 1, ]
      }

      as.bugs.array2 <- utils::getFromNamespace("as.bugs.array2", "R2jags")
      jags.out <- list(
        model = lapply(chain_fits, function(chain_fit) chain_fit$model),
        BUGSoutput = as.bugs.array2(
          sims.array,
          model.file = parallelfile,
          program = "jags",
          DIC = TRUE,
          n.iter = iter,
          n.burnin = warmup,
          n.thin = thin
        ),
        parameters.to.save = jags.params,
        model.file = parallelfile,
        n.iter = iter,
        DIC = TRUE
      )
      class(jags.out) <- c("rjags.parallel", "rjags")
    } else {
      set.seed(seed)
      jags.out <- R2jags::jags(
        data = jags.data,
        inits = jags.inits,
        n.chains = chains,
        parameters.to.save = jags.params,
        n.iter = iter,
        n.burnin = warmup,
        n.thin = thin,
        n.adapt = n.adapt,
        model.file = textConnection(modelstring)
      )
    }

    # Format JAGS output ---------------------------------------------------------------- #

    formatted <- formatJags(
      jags.out,
      monitor_spec,
      Ns,
      Xs,
      mm_blocks,
      main,
      hm_blocks,
      interactions,
      family_str,
      cox_intervals
    )

    # Prepare metadata ------------------------------------------------------------------- #

    mm_info <- lapply(seq_along(mm_blocks %||% list()), function(k) {
      list(
        name = mm_blocks[[k]]$name,
        feature_labels = mm_blocks[[k]]$feature_labels,
        vars = mm_blocks[[k]]$vars,
        fn = list(type = mm_blocks[[k]]$fn$type,
                  est_params = names(mm_blocks[[k]]$fn$est_params)),
        w = list(params = mm_blocks[[k]]$w$params, scale = mm_blocks[[k]]$w$scale),
        RE = mm_blocks[[k]]$RE,
        FE = mm_blocks[[k]]$FE,
        mmid_group = mm_blocks[[k]]$mmid_group
      )
    })

    hm_info <- lapply(seq_along(hm_blocks %||% list()), function(k) {
      list(
        id = hm_blocks[[k]]$id,
        RE = hm_blocks[[k]]$RE,
        FE = hm_blocks[[k]]$FE,
        labels = hm_blocks[[k]]$labels
      )
    })

    input <- list(
      formula = formula,
      family = family_str,
      family_obj = fam,
      prior = prior,
      prior_table = prior_table,
      inits = inits,
      iter = iter,
      warmup = warmup,
      thin = thin,
      chains = chains,
      cores = cores,
      seed = seed,
      monitor = monitor_spec$requested,
      modelfile = modelfile,
      run = run,
      lhs = main$lhs,
      mainvars = main$vars,
      X.main.ranges = if (!is.null(Xs$X.main)) {
        r <- t(apply(Xs$X.main, 2, range))
        colnames(r) <- c("min", "max")
        r
      } else NULL,
      X.main.means = if (!is.null(Xs$X.main)) as.list(colMeans(Xs$X.main)) else NULL,
      mm = mm_info,
      hm = hm_info,
      interactions = interactions,
      y = if (family_str %in% c("Gaussian", "Binomial")) Ys$Y else NULL,
      n.umm = Ns$n.umm,
      n.mm = Ns$n.mm,
      n.main = Ns$n.main,
      n.hm = Ns$n.hm,
      n.mmblocks = Ns$n.mmblocks
    )

    diagnostics <- bml_parameter_diagnostics(
      jags.out = jags.out,
      reg.table = formatted$reg.table,
      chains = chains
    )

    retained_variables <- bml_retained_draw_variables(
      jags.out = jags.out,
      reg.table = formatted$reg.table,
      monitor_spec = monitor_spec
    )
    draws <- if (isTRUE(monitor_spec$raw_jags)) {
      NULL
    } else {
      bml_select_draws(jags.out, retained_variables)
    }
    raw_variables <- dimnames(bml_raw_sims_array(jags.out))[[3]]
    stored_variables <- if (isTRUE(monitor_spec$raw_jags)) {
      raw_variables
    } else if (!is.null(draws)) {
      posterior::variables(draws)
    } else {
      character()
    }
    n_draws <- if (isTRUE(monitor_spec$raw_jags)) {
      prod(dim(bml_raw_sims_array(jags.out))[1:2])
    } else if (!is.null(draws)) {
      posterior::ndraws(draws)
    } else {
      0L
    }

    storage <- list(
      requested = monitor_spec$requested,
      capabilities = monitor_spec$capabilities,
      label = bml_storage_label(monitor_spec),
      backend = if (isTRUE(monitor_spec$raw_jags)) "R2jags" else
        if (is.null(draws)) "summary" else "posterior::draws_array",
      variables = stored_variables,
      ndraws = as.integer(n_draws)
    )

    out <- list(
      reg.table = formatted$reg.table,
      w = formatted$w,
      re.mm = formatted$re.mm,
      re.hm = formatted$re.hm,
      pred = formatted$pred,
      fitted = formatted$fitted,
      input = input,
      diagnostics = diagnostics,
      storage = storage,
      draws = draws,
      jags.out = if (isTRUE(monitor_spec$raw_jags)) jags.out else NULL
    )

    class(out) <- "bml"

    return(out)
  } else {
    message("Data and model have been created without any errors.")
    invisible(list(
      modelstring = modelstring,
      jags.data = jags.data,
      jags.params = jags.params,
      jags.inits = jags.inits,
      storage = monitor_spec
    ))
  }
}
