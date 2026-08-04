# ================================================================================================ #
# fn(): the aggregation function f in theta^{micro,f}(M_it), plus est() and the expression DSL
# ================================================================================================ #

#' Mark a shape parameter as estimated
#'
#' @description
#' Used inside \code{\link{fn}}'s named-shortcut arguments to say "estimate this
#' parameter" rather than fixing it to a number:
#' \code{fn("threshold", c = est(), kappa = 10)}. In the formula DSL a free
#' symbol alone declares a parameter, so \code{est()} is only needed where a
#' value slot must distinguish estimate-vs-fix.
#'
#' @param lower,upper Optional bounds; if both given, the default prior is
#'   \code{dunif(lower, upper)}. Otherwise the parameter gets a
#'   parameter-specific default (e.g. the observed range of \code{x} for a
#'   threshold cutpoint). Override via \code{prior(..., class = "fn")}.
#'
#' @return A \code{bml_est} marker object.
#' @seealso \code{\link{fn}}, \code{\link{prior}}
#' @export
est <- function(lower = NULL, upper = NULL) {
  structure(list(lower = lower, upper = upper), class = "bml_est")
}

is_est <- function(x) inherits(x, "bml_est")

# Named fn types and their attribute arity
FN_TYPES <- c("sum", "var", "hhi", "effn", "entropy", "threshold", "smax", "gmean", "cov")
FN_ARITY <- c(sum = NA, var = 1, hhi = 0, effn = 0, entropy = 0,
              threshold = 1, smax = 1, gmean = 1, cov = 2)

# Whitelisted operations for the expression DSL (transpilable to JAGS)
DSL_OPS <- c("+", "-", "*", "/", "^", "(", "E", "exp", "log", "ilogit", "pow")

#' The aggregation function of an mm() block
#'
#' @description
#' Selects the function \code{f} that reduces the weighted member records to a
#' group-level feature — the \code{f} in the framework's
#' \eqn{\theta^{micro,f}(M_{it})}. \code{fn("sum")} (the default) is the
#' additive weighted mean; the other named types are "emergent" features of the
#' whole member set. \code{fn()} also accepts a one-sided formula written in a
#' restricted expression DSL, so users can define their own reductions.
#'
#' \strong{Named types:}
#' \preformatted{
#' fn("sum")                          # A_x  = sum_k w_k x_k          (default)
#' fn("var")                          # V_x  = sum_k w_k (x_k - A_x)^2
#' fn("var", moment = 3)              # higher central moments
#' fn("hhi")                          # C    = sum_k w_k^2  (weights only)
#' fn("effn")                         # 1 / C
#' fn("entropy")                      # -sum_k w_k log w_k
#' fn("threshold", c = 0.7, kappa = 10)  # T(c) = sum_k w_k ilogit(kappa (x_k - c))
#' fn("threshold", c = est())         # estimate the cutpoint
#' fn("smax", kappa = 5)              # (1/kappa) log sum_k w_k exp(kappa x_k)
#' fn("smax", kappa = est())          # estimate the aggregation function itself:
#'                                    #   kappa<0 -> min, kappa->0 -> mean, kappa>0 -> max
#' fn("gmean", p = est())             # power/CES mean (sum_k w_k x_k^p)^(1/p); x > 0
#' fn("cov")                          # C_xz = sum_k w_k (x_k - A_x)(z_k - A_z); two attributes
#' }
#'
#' \strong{Expression DSL:} \code{w()} normalizes weights to a probability
#' measure over the group's members; \code{E(e)} is the expectation
#' \eqn{\sum_k w_k e_k} under that measure. Member-level quantities are the
#' \code{vars()} attributes and the reserved symbol \code{w}; anything wrapped
#' in \code{E()} is a group scalar. Whitelisted operations: \code{+ - * / ^},
#' \code{exp}, \code{log}, \code{ilogit}, \code{pow}. Any symbol that is not a
#' data column or reserved word is a \emph{free parameter} with a default prior
#' (one rule, shared with \code{\link{w}}); the build messages the detected
#' parameters.
#' \preformatted{
#' fn(~ E((x - E(x))^2))                       # variance, written out
#' fn(~ E(ilogit(kappa * (x - c))))            # threshold with free c, kappa
#' fn(~ E((x - E(x)) * (z - E(z))))            # covariance
#' }
#'
#' \strong{Identification:} internal parameters reach the outcome only through
#' \code{b * feature}; if the feature's coefficient is ~0 they are unidentified.
#' Dispersion functions need real within-group spread; threshold/tail functions
#' need mass of \code{x} in the region they read. A free parameter that
#' multiplies an attribute inside a nonlinear \code{E()} is confounded with the
#' feature's coefficient — treat inner parameters as shape parameters
#' (cutpoints, sharpness), not slopes.
#'
#' @param type A type string (one of \code{"sum"}, \code{"var"}, \code{"hhi"},
#'   \code{"effn"}, \code{"entropy"}, \code{"threshold"}, \code{"smax"},
#'   \code{"gmean"}, \code{"cov"}) or a one-sided formula (\code{~ E(...)}).
#' @param ... Shape parameters for the named types: \code{moment} (var),
#'   \code{c} and \code{kappa} (threshold), \code{kappa} (smax), \code{p}
#'   (gmean). Each is a number (fixed) or \code{\link{est}()} (estimated).
#'
#' @return A \code{bml_fn} object.
#' @seealso \code{\link{mm}}, \code{\link{w}}, \code{\link{est}}, \code{\link{prior}}
#' @export
fn <- function(type = "sum", ...) {

  dots <- list(...)

  # ------------------------------------------------------------------ #
  # Formula: the expression DSL
  # ------------------------------------------------------------------ #
  if (inherits(type, "formula")) {
    if (length(type) != 2) {
      if (identical(deparse(type[[2]]), "w")) {
        stop("Weights moved to w(). The old fn(w ~ ...) spelling was removed;\n",
             "  write w = w(~ ...) for weights and fn = fn(\"...\") for the aggregation ",
             "function.", call. = FALSE)
      }
      stop("fn() formulas must be one-sided: fn(~ E(...)).", call. = FALSE)
    }
    expr <- type[[2]]
    validate_dsl_ops(expr)
    if (!dsl_contains_E(expr)) {
      stop("fn() expressions must contain at least one E(...) reduction ",
           "(the top level must be a group-level scalar).", call. = FALSE)
    }
    if (length(dots) > 0 && is.null(names(dots))) {
      stop("Additional fn() arguments must be named.", call. = FALSE)
    }
    return(structure(
      list(type = "expr", expr = expr, moment = NULL, params = dots),
      class = "bml_fn"
    ))
  }

  # ------------------------------------------------------------------ #
  # Old weights spelling: removed with a pointed error
  # ------------------------------------------------------------------ #
  if (is.character(type) && length(type) == 1 && type %in% FN_TYPES) {
    fn_type <- type
  } else {
    stop("fn() takes a type string (", paste(FN_TYPES, collapse = ", "),
         ") or a one-sided formula (~ E(...)).\n",
         "  Note: weights moved to w(). The old fn(w ~ ...) spelling was removed;\n",
         "  write w = w(~ ...) for weights and fn = fn(\"...\") for the aggregation function.",
         call. = FALSE)
  }

  # ------------------------------------------------------------------ #
  # Named types: collect shape parameters
  # ------------------------------------------------------------------ #
  params <- list()
  moment <- NULL

  if (fn_type == "var") {
    moment <- dots$moment %||% 2
    if (is_est(moment)) stop("fn(\"var\"): 'moment' must be a fixed integer >= 2.", call. = FALSE)
    if (!is.numeric(moment) || moment < 2 || moment != round(moment)) {
      stop("fn(\"var\"): 'moment' must be an integer >= 2.", call. = FALSE)
    }
    dots$moment <- NULL
  } else if (fn_type == "threshold") {
    if (is.null(dots$c)) {
      stop("fn(\"threshold\") requires a cutpoint: c = <number> or c = est().", call. = FALSE)
    }
    params$c <- dots$c
    params$kappa <- dots$kappa %||% 10
    dots$c <- NULL; dots$kappa <- NULL
    if (is_est(params$c) && is_est(params$kappa)) {
      warning("Estimating both c and kappa in fn(\"threshold\") is only weakly identified ",
              "(a shallower slope plus a shifted cutpoint mimics a sharper one). ",
              "Consider fixing kappa.", call. = FALSE)
    }
  } else if (fn_type == "smax") {
    if (is.null(dots$kappa)) {
      stop("fn(\"smax\") requires kappa = <number> or kappa = est().", call. = FALSE)
    }
    params$kappa <- dots$kappa
    dots$kappa <- NULL
    if (is.numeric(params$kappa) && params$kappa == 0) {
      stop("fn(\"smax\"): kappa must be nonzero (kappa -> 0 is the mean; use fn(\"sum\")).",
           call. = FALSE)
    }
  } else if (fn_type == "gmean") {
    if (is.null(dots$p)) {
      stop("fn(\"gmean\") requires p = <number> or p = est().", call. = FALSE)
    }
    params$p <- dots$p
    dots$p <- NULL
    if (is.numeric(params$p) && params$p == 0) {
      stop("fn(\"gmean\"): p must be nonzero (p -> 0 is the geometric mean; ",
           "use fn(~ E(log(x))) for the log-mean).", call. = FALSE)
    }
  }

  if (length(dots) > 0) {
    stop("Unknown fn(\"", fn_type, "\") argument(s): ",
         paste(names(dots), collapse = ", "), call. = FALSE)
  }

  structure(
    list(type = fn_type, expr = NULL, moment = moment, params = params),
    class = "bml_fn"
  )
}

# ------------------------------------------------------------------------------------------------ #
# DSL utilities
# ------------------------------------------------------------------------------------------------ #

validate_dsl_ops <- function(expr) {
  walk <- function(e) {
    if (is.call(e)) {
      op <- as.character(e[[1]])
      if (!op %in% DSL_OPS) {
        stop("fn() expression uses non-whitelisted operation '", op, "'.\n",
             "  Allowed: ", paste(setdiff(DSL_OPS, "("), collapse = ", "), call. = FALSE)
      }
      for (i in seq_along(e)[-1]) walk(e[[i]])
    }
  }
  walk(expr)
  invisible(TRUE)
}

dsl_contains_E <- function(expr) {
  found <- FALSE
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(as.character(e[[1]]), "E")) found <<- TRUE
      for (i in seq_along(e)[-1]) walk(e[[i]])
    }
  }
  walk(expr)
  found
}

# All symbols appearing in a DSL expression
dsl_symbols <- function(expr) {
  syms <- character(0)
  walk <- function(e) {
    if (is.name(e)) {
      syms <<- c(syms, as.character(e))
    } else if (is.call(e)) {
      for (i in seq_along(e)[-1]) walk(e[[i]])
    }
  }
  walk(expr)
  unique(syms)
}

# Decompose a DSL expression into a dependency-ordered list of E-nodes plus a
# top-level group-scalar expression. E(...) calls are replaced by placeholder
# symbols .E1, .E2, ... (inner nodes first, so bodies may reference earlier
# placeholders).
fn_expr_graph <- function(expr) {

  enodes <- list()

  transform <- function(e) {
    if (is.call(e)) {
      if (identical(as.character(e[[1]]), "E")) {
        if (length(e) != 2) {
          stop("E() takes exactly one argument.", call. = FALSE)
        }
        body <- transform(e[[2]])
        id <- length(enodes) + 1
        enodes[[id]] <<- list(id = id, body = body)
        return(as.name(paste0(".E", id)))
      }
      for (i in seq_along(e)[-1]) e[[i]] <- transform(e[[i]])
      return(e)
    }
    e
  }

  top <- transform(expr)
  list(enodes = enodes, top = top)
}

# Validate the resolved DSL: symbols classified into attributes (member-level
# data columns), the reserved weight symbol w, and free parameters.
# - top level must be a group scalar (no member-level symbols outside E)
# - warn when a free parameter enters multiplicatively with an attribute
validate_dsl_resolved <- function(graph, attrs, params) {

  member_syms <- c(attrs, "w")

  # Top level must be group scalar
  top_syms <- dsl_symbols(graph$top)
  bad <- intersect(top_syms, member_syms)
  if (length(bad) > 0) {
    stop("fn() expression: member-level quantities (", paste(bad, collapse = ", "),
         ") appear outside E(). The top level must be a group-level scalar - wrap ",
         "member-level terms in E(...).", call. = FALSE)
  }

  # Multiplicative param x attribute warning (confounded with the feature coefficient)
  warn_mult <- FALSE
  walk <- function(e) {
    if (is.call(e)) {
      op <- as.character(e[[1]])
      if (op %in% c("*", "/") && length(e) == 3) {
        l <- dsl_symbols(e[[2]]); r <- dsl_symbols(e[[3]])
        if ((any(l %in% params) && any(r %in% attrs)) ||
            (any(r %in% params) && any(l %in% attrs))) {
          warn_mult <<- TRUE
        }
      }
      for (i in seq_along(e)[-1]) walk(e[[i]])
    }
  }
  for (node in graph$enodes) walk(node$body)

  if (warn_mult) {
    warning("fn() expression: a free parameter multiplies an attribute inside E(). ",
            "Such parameters are confounded with the feature's coefficient (joint scaling ",
            "non-identification). Treat inner parameters as shape parameters (cutpoints, ",
            "sharpness), not slopes; consider fixing the parameter.", call. = FALSE)
  }

  invisible(TRUE)
}

# Deparse a DSL (sub)expression into a JAGS-compatible string with symbol
# substitutions applied. `subs` is a named character vector symbol -> JAGS code.
dsl_deparse_jags <- function(e, subs) {
  if (is.name(e)) {
    s <- as.character(e)
    if (s %in% names(subs)) return(subs[[s]])
    return(s)
  }
  if (is.numeric(e) || is.logical(e)) return(format(e))
  if (is.call(e)) {
    op <- as.character(e[[1]])
    args <- lapply(as.list(e)[-1], dsl_deparse_jags, subs = subs)
    if (op == "(") return(paste0("(", args[[1]], ")"))
    if (op == "-" && length(args) == 1) return(paste0("(-", args[[1]], ")"))
    if (op %in% c("+", "-", "*", "/")) {
      return(paste0("(", args[[1]], " ", op, " ", args[[2]], ")"))
    }
    if (op == "^") return(paste0("pow(", args[[1]], ", ", args[[2]], ")"))
    if (op %in% c("exp", "log", "ilogit")) return(paste0(op, "(", args[[1]], ")"))
    if (op == "pow") return(paste0("pow(", args[[1]], ", ", args[[2]], ")"))
    stop("Internal error: unexpected DSL operation '", op, "'.", call. = FALSE)
  }
  stop("Internal error: unexpected DSL node.", call. = FALSE)
}

# Evaluate a DSL graph in R for one group (used for the precompute fast path).
# `member_env`: list with attribute vectors and w (normalized weights).
# `param_values`: named list of fixed parameter values.
dsl_eval_group <- function(graph, member_env, param_values) {

  base_funs <- list(
    ilogit = stats::plogis,
    pow = function(a, b) a^b
  )

  evals <- list()
  for (node in graph$enodes) {
    env_list <- c(member_env, param_values, evals, base_funs)
    e_i <- eval(node$body, envir = env_list, enclos = baseenv())
    evals[[paste0(".E", node$id)]] <- sum(member_env$w * e_i)
  }
  env_list <- c(param_values, evals, base_funs)
  eval(graph$top, envir = env_list, enclos = baseenv())
}

# ------------------------------------------------------------------------------------------------ #
# Shape-parameter helpers (shared by codegen and prior table)
# ------------------------------------------------------------------------------------------------ #

# Extract the estimated shape parameters of a resolved fn:
# returns named list param -> list(default = <JAGS dist>, est = bml_est or NULL)
fn_est_params <- function(fnobj, dsl_params = character(0), x_range = NULL) {

  out <- list()

  default_for <- function(name, marker) {
    if (is_est(marker) && !is.null(marker$lower) && !is.null(marker$upper)) {
      return(sprintf("dunif(%s, %s)", format(marker$lower), format(marker$upper)))
    }
    switch(
      name,
      "c" = if (!is.null(x_range)) {
        sprintf("dunif(%s, %s)", format(x_range[1]), format(x_range[2]))
      } else "dnorm(0, 0.01)",
      # threshold sharpness is positive; smax kappa is real-valued (kappa < 0 = min,
      # kappa > 0 = max) so the adjudication needs support on both sides of 0
      "kappa" = if (identical(fnobj$type, "smax")) "dnorm(0, 0.04)" else "dgamma(2, 0.2)",
      "p" = "dnorm(1, 0.1)",        # centered on the arithmetic mean
      "dnorm(0, 0.01)"              # DSL free symbols
    )
  }

  init_for <- function(name, marker) {
    if (is_est(marker) && !is.null(marker$lower) && !is.null(marker$upper)) {
      return(mean(c(marker$lower, marker$upper)))
    }
    switch(
      name,
      "c" = if (!is.null(x_range)) mean(x_range) else 0,
      # smax kappa must not initialize at 0 (the 1/kappa limit); 1 is a mild max-lean
      "kappa" = if (identical(fnobj$type, "smax")) 1 else 5,
      "p" = 1,
      0
    )
  }

  # Parameters that sit at a singularity when they hit zero, so dispersed starting values
  # must be kept away from it: smax is (1/kappa) log(...) and gmean is (...)^(1/p).
  excludes_zero <- function(name) {
    (identical(fnobj$type, "smax") && identical(name, "kappa")) ||
      (identical(fnobj$type, "gmean") && identical(name, "p"))
  }

  # named-type est() markers
  for (pn in names(fnobj$params %||% list())) {
    if (is_est(fnobj$params[[pn]])) {
      out[[pn]] <- list(default = default_for(pn, fnobj$params[[pn]]),
                        init = init_for(pn, fnobj$params[[pn]]),
                        exclude_zero = excludes_zero(pn),
                        est = fnobj$params[[pn]])
    }
  }

  # DSL free symbols are parameters by the one-parameter rule
  for (pn in dsl_params) {
    out[[pn]] <- list(default = default_for(pn, NULL), init = init_for(pn, NULL),
                      exclude_zero = FALSE, est = NULL)
  }

  out
}

# Auto-generated feature label(s) for a block
fn_feature_labels <- function(fnobj, attr_cols) {
  clean <- function(x) gsub(":", "_", x)
  switch(
    fnobj$type,
    "sum" = ifelse(grepl(":", attr_cols),
                   paste0("H_", clean(attr_cols)),
                   paste0("A_", attr_cols)),
    "var" = if ((fnobj$moment %||% 2) == 2) paste0("V_", attr_cols[1])
            else paste0("M", fnobj$moment, "_", attr_cols[1]),
    "hhi" = "C_w",
    "effn" = "effN",
    "entropy" = "H_w",
    "threshold" = paste0("T_", attr_cols[1]),
    "smax" = paste0("smax_", attr_cols[1]),
    "gmean" = paste0("M_", attr_cols[1]),
    "cov" = paste0("C_", attr_cols[1], "_", attr_cols[2]),
    "expr" = "f",
    "F"
  )
}
