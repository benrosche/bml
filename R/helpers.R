# ================================================================================================ #
# Helper functions for bml formula syntax
# ================================================================================================ #

#' Specify identifier variables for multiple-membership and hierarchical structures
#'
#' @description
#' Helper function used within \code{\link{mm}} and \code{\link{hm}} to specify
#' the identifier variables that define memberships and nesting structures. In
#' multiple-membership models, \code{id()} links member-level units (e.g., task
#' IDs) to group-level units (e.g., occupation IDs) in long format: one row per
#' membership, any number of members per group. In hierarchical models,
#' \code{id()} specifies the nesting-level identifier (e.g., country ID).
#'
#' @param ... Unquoted variable names from your data:
#'   \itemize{
#'     \item For \code{mm()}: Two identifiers \code{id(mmid, mainid)} where
#'       \code{mmid} identifies member-level units and \code{mainid} identifies
#'       group-level units
#'     \item For \code{hm()}: One identifier \code{id(hmid)} where \code{hmid}
#'       identifies nesting-level units
#'   }
#'
#' @return A \code{bml_id} object containing the variable names as character strings.
#'
#' @seealso \code{\link{mm}}, \code{\link{hm}}, \code{\link{bml}}
#'
#' @examples
#' \dontrun{
#' # Multiple-membership: parties (pid) within governments (gid)
#' id(pid, gid)
#'
#' # Hierarchical: governments within countries
#' id(cid)
#' }
#'
#' @export
id <- function(...) {
  ids <- as.character(match.call(expand.dots = FALSE)$...)
  structure(ids, class = "bml_id")
}

#' Fix a coefficient to a known value
#'
#' @description
#' Specify a covariate whose coefficient should be held constant at a fixed value
#' rather than estimated from the data. This is useful for offset variables or
#' when you want to impose theoretical constraints. Fixed coefficients are
#' handled efficiently by pre-computing their contribution in R before passing
#' data to JAGS.
#'
#' @param var Unquoted variable name from your data
#' @param value Numeric value for the coefficient (e.g., \code{1.0} for a standard offset)
#'
#' @return A \code{bml_fix} object that can be used within \code{\link{vars}}.
#'
#' @seealso \code{\link{vars}}, \code{\link{mm}}
#'
#' @examples
#' \dontrun{
#' # Fix a coefficient to 1.0 (standard offset)
#' fix(exposure, 1.0)
#'
#' # Use within vars() for multiple-membership models
#' vars(fix(population, 0.5) + income + education)
#' }
#'
#' @export
fix <- function(var, value) {
  var_name <- deparse(substitute(var))
  structure(list(var = var_name, value = as.numeric(value)), class = "bml_fix")
}

#' Specify member-level attributes for mm() blocks
#'
#' @description
#' Helper function used within \code{\link{mm}} to specify the member-level
#' attributes (the effect covariates \eqn{x_{kt}}) that the block's aggregation
#' function \code{\link{fn}} reads. Supports both free variables (with
#' coefficients estimated by the main model) and fixed variables (coefficients
#' held constant via \code{\link{fix}}).
#'
#' @param ... Unquoted variable names, combined with \code{+} (formula-style)
#'   or as separate arguments. Supports:
#'   \itemize{
#'     \item Simple variables: \code{vars(x + y)} or \code{vars(x, y)}
#'     \item Interactions: \code{vars(x * y)} (main effects + interaction) or
#'       \code{vars(x:y)} (interaction only; the member-paired interaction
#'       \eqn{H_{xz} = \sum_k w_k x_k z_k})
#'     \item Transformations: \code{vars(I(x^2))}
#'     \item Fixed coefficients: \code{vars(fix(x, 1.0) + y)}
#'   }
#'   Note: Numeric literals like \code{1}, \code{0}, or \code{-1} are ignored
#'   (no intercept support in mm blocks).
#'
#' @return A \code{bml_vars} object (list with \code{formula}, \code{free},
#'   \code{fixed}), or \code{NULL} if empty.
#'
#' @seealso \code{\link{fix}}, \code{\link{mm}}, \code{\link{fn}}
#'
#' @examples
#' \dontrun{
#' vars(income)
#' vars(income + education)
#' vars(income:education)      # member-paired interaction only
#' vars(x, z)                  # two attributes (e.g. for fn("cov"))
#' vars(fix(exposure, 1.0) + income)
#' }
#'
#' @export
vars <- function(...) {
  exprs <- as.list(substitute(list(...)))[-1]
  if (length(exprs) == 0) return(NULL)

  fixed_list <- list()
  free_terms <- character()

  extract_and_build <- function(e) {
    if (is.call(e)) {
      fn_name <- as.character(e[[1]])
      if (fn_name == "fix") {
        var_name <- as.character(e[[2]])
        value <- eval(e[[3]], parent.frame(2))
        fixed_list[[length(fixed_list) + 1]] <<- list(var = var_name, value = value)
        return(NULL)
      } else if (fn_name == "+") {
        left <- extract_and_build(e[[2]])
        right <- extract_and_build(e[[3]])
        return(c(left, right))
      } else {
        term_str <- paste(deparse(e, width.cutoff = 500), collapse = "")
        return(term_str)
      }
    } else if (is.name(e) && as.character(e) != "") {
      return(as.character(e))
    } else if (is.numeric(e)) {
      return(NULL) # numeric literals (1, 0, -1) ignored
    }
    NULL
  }

  for (e in exprs) {
    free_terms <- c(free_terms, extract_and_build(e))
  }
  free_terms <- unique(free_terms[!is.null(free_terms) & free_terms != ""])

  if (length(free_terms) == 0 && length(fixed_list) == 0) return(NULL)

  if (length(free_terms) > 0) {
    formula_str <- paste("~ 0 +", paste(free_terms, collapse = " + "))
    vars_formula <- as.formula(formula_str)
    base_vars <- all.vars(vars_formula)
  } else {
    vars_formula <- NULL
    base_vars <- character()
  }

  structure(
    list(
      formula = vars_formula,
      free = base_vars,
      fixed = if (length(fixed_list) > 0) fixed_list else NULL
    ),
    class = "bml_vars"
  )
}

#' Specify the membership weights of an mm() block
#'
#' @description
#' Defines the weights \eqn{w_{ikt}} — how much each member counts in the
#' block's aggregation. \code{w()} defines a within-group \emph{measure} over
#' members; the aggregation function \code{\link{fn}} chooses which function of
#' that measure (and of the \code{\link{vars}} attributes) to read off. Weights
#' can be fixed, taken from observed data, or estimated from member
#' characteristics:
#'
#' \preformatted{
#' w(~ 1/n)                          # equal weights
#' w(~ importance, scale = TRUE)     # observed, normalized
#' w(~ ilogit(b0 + b1 * q))          # estimated: b0, b1 are free parameters
#' }
#'
#' \strong{One parameter rule (shared with \code{\link{fn}}):} any symbol in the
#' formula that is not a data column or reserved word (\code{n}, the group
#' size) is a \emph{free parameter} with a default prior. The build always
#' messages the detected parameters (e.g. \code{Estimating parameters: b0, b1
#' (w)}), so a misspelled column name shows up immediately. Priors are set via
#' \code{prior(..., class = "w")}.
#'
#' \strong{Who counts vs. how contributions combine:} parameters that act on
#' \emph{external} member characteristics (importance, exposure time) belong in
#' \code{w()}; parameters that act on the \emph{effect attribute} (whatever is
#' in \code{vars()}) belong in \code{\link{fn}}. Building weights from the same
#' variable that is being aggregated (e.g. \code{w(~ exp(om * x))} with
#' \code{vars(x)}) is near-equivalent to a smooth-max in \code{fn()} and is
#' flagged with a warning.
#'
#' @param formula One-sided formula defining the raw weights. Group-level
#'   aggregation shortcuts (\code{min(x)}, \code{max(x)}, \code{mean(x)},
#'   \code{sum(x)}, \code{sd(x)}, \code{var(x)}, \code{median(x)},
#'   \code{mode(x)}, \code{range(x)}, \code{first(x)}, \code{last(x)},
#'   \code{quantile(x, p)}) are pre-computed in R. JAGS math functions
#'   (\code{exp}, \code{log}, \code{ilogit}, ...) pass through.
#' @param scale Logical; if \code{TRUE} (default), weights are normalized to
#'   sum to 1 within each group, making them a probability measure over the
#'   group's members. Emergent \code{fn()} types require \code{scale = TRUE}.
#'
#' @return A \code{bml_w} object.
#'
#' @details
#' Weight functions with estimated parameters must produce bounded, positive
#' values across plausible parameter values; unbounded weight functions can
#' crash the sampler. Use bounded forms such as
#' \code{w(~ ilogit(b0 + b1 * q), scale = TRUE)} or the generalized logistic
#' \code{w(~ 1 / (1 + (n - 1) * exp(-(b0 + b1 * q))), scale = TRUE)}, keep
#' \code{scale = TRUE}, standardize covariates, and use informative priors
#' (\code{prior(normal(0, 1), class = "w")}). Weight parameters are initialized
#' at 0 for stable starting values.
#'
#' @seealso \code{\link{mm}}, \code{\link{fn}}, \code{\link{bml}}
#'
#' @examples
#' \dontrun{
#' w(~ 1/n)                                # equal weights
#' w(~ importance, scale = TRUE)           # observed weights, normalized
#' w(~ ilogit(b0 + b1 * seniority))        # estimated weights
#' w(~ duration, scale = FALSE)            # unnormalized (aggregating sums)
#' }
#'
#' @references
#' Rosche, B. (2026). A Multilevel Model for Coalition Governments: Uncovering
#' Party-Level Dependencies Within and Between Governments. \emph{Political Analysis}.
#'
#' @export
w <- function(formula = ~ 1/n, scale = TRUE) {

  if (!inherits(formula, "formula")) {
    stop("w() takes a one-sided formula, e.g. w(~ 1/n) or w(~ importance).", call. = FALSE)
  }
  if (length(formula) == 3) {
    stop("w() takes a one-sided formula: the old 'w ~ ...' left-hand side is implied.\n",
         "  Old: fn(w ~ importance, c = TRUE)\n",
         "  New: w(~ importance, scale = TRUE)", call. = FALSE)
  }

  fn_string <- paste(deparse(formula[[2]], width.cutoff = 500), collapse = "")

  # ------------------------------------------------------------------------------------ #
  # Group-level aggregation shortcuts (pre-computed in R)
  # ------------------------------------------------------------------------------------ #

  supported_agg_simple <- c("min", "max", "mean", "sum", "sd", "var", "first", "last",
                            "median", "mode", "range")
  agg_pattern_simple <- "\\b(min|max|mean|sum|sd|var|first|last|median|mode|range)\\s*\\(\\s*([a-zA-Z_][a-zA-Z0-9_\\.]*)\\s*\\)"
  agg_pattern_quantile <- "\\bquantile\\s*\\(\\s*([a-zA-Z_][a-zA-Z0-9_\\.]*)\\s*,\\s*(-?[0-9]*\\.?[0-9]+)\\s*\\)"
  supported_agg <- c(supported_agg_simple, "quantile")

  jags_math_funcs <- c("exp", "log", "log10", "sqrt", "abs", "pow",
                       "sin", "cos", "tan", "asin", "acos", "atan",
                       "sinh", "cosh", "tanh",
                       "logit", "ilogit", "probit", "iprobit",
                       "cloglog", "icloglog",
                       "round", "trunc", "floor", "ceiling")
  all_allowed_funcs <- c(supported_agg, jags_math_funcs)

  nested_pattern <- "\\b(min|max|mean|sum|sd|var|first|last|median|mode|range|quantile)\\s*\\([^)]*\\b(min|max|mean|sum|sd|var|first|last|median|mode|range|quantile)\\s*\\("
  if (grepl(nested_pattern, fn_string)) {
    stop("Nested aggregation functions are not supported in weight formulas.\n",
         "  Example of invalid: w(~ min(max(x)))\n",
         "  Valid: w(~ min(x) + max(y))", call. = FALSE)
  }

  all_func_calls <- stringr::str_match_all(fn_string, "\\b([a-zA-Z_][a-zA-Z0-9_]*)\\s*\\(")[[1]]
  if (nrow(all_func_calls) > 0) {
    func_names <- all_func_calls[, 2]
    unsupported <- setdiff(func_names, all_allowed_funcs)
    if (length(unsupported) > 0) {
      stop("Unsupported function(s) in weight formula: ", paste(unsupported, collapse = ", "), "\n",
           "  Supported aggregation functions: ", paste(supported_agg, collapse = ", "), "\n",
           "  Supported JAGS math functions: ", paste(jags_math_funcs, collapse = ", "), call. = FALSE)
    }
  }

  agg_matches_simple <- stringr::str_match_all(fn_string, agg_pattern_simple)[[1]]
  agg_matches_quantile <- stringr::str_match_all(fn_string, agg_pattern_quantile)[[1]]

  agg_funcs <- NULL
  agg_vars <- NULL

  if (nrow(agg_matches_simple) > 0) {
    agg_funcs <- list()
    for (i in 1:nrow(agg_matches_simple)) {
      agg_funcs[[length(agg_funcs) + 1]] <- list(
        original = agg_matches_simple[i, 1],
        func     = agg_matches_simple[i, 2],
        var      = agg_matches_simple[i, 3],
        prob     = NULL,
        col_name = paste0(agg_matches_simple[i, 3], "_", agg_matches_simple[i, 2])
      )
    }
  }

  if (nrow(agg_matches_quantile) > 0) {
    if (is.null(agg_funcs)) agg_funcs <- list()
    for (i in 1:nrow(agg_matches_quantile)) {
      var_name <- agg_matches_quantile[i, 2]
      prob_val <- as.numeric(agg_matches_quantile[i, 3])
      if (prob_val < 0 || prob_val > 1) {
        stop("quantile() probability must be between 0 and 1. Got: ", prob_val, call. = FALSE)
      }
      prob_pct <- round(prob_val * 100)
      col_name <- paste0(var_name, "_q", prob_pct)
      agg_funcs[[length(agg_funcs) + 1]] <- list(
        original = agg_matches_quantile[i, 1],
        func     = "quantile",
        var      = var_name,
        prob     = prob_val,
        col_name = col_name
      )
    }
  }

  if (!is.null(agg_funcs)) {
    agg_vars <- unique(sapply(agg_funcs, function(x) x$var))
  }

  fn_string_transformed <- fn_string
  if (!is.null(agg_funcs)) {
    for (agg in agg_funcs) {
      escaped <- gsub("([\\(\\)\\.])", "\\\\\\1", agg$original)
      fn_string_transformed <- gsub(escaped, agg$col_name, fn_string_transformed)
    }
  }

  # ------------------------------------------------------------------------------------ #
  # Collect symbols; parameter/variable split is resolved against the data later
  # (one parameter rule: free symbols become parameters)
  # ------------------------------------------------------------------------------------ #

  all_symbols <- all.vars(as.formula(paste("~", fn_string_transformed)))

  structure(
    list(
      formula    = formula,
      string     = fn_string_transformed,
      symbols    = all_symbols,
      vars       = NULL,   # resolved in dissectFormula (data columns + agg cols + n)
      vars_p     = NULL,   # resolved in dissectFormula (b * x pattern, for labels)
      params     = NULL,   # resolved in dissectFormula (free symbols)
      scale      = isTRUE(scale),
      constraint = isTRUE(scale),  # internal alias used by the codegen
      agg_funcs  = agg_funcs,
      agg_vars   = agg_vars
    ),
    class = "bml_w"
  )
}

# Internal: defaults evaluated in the package namespace (the mm() arguments
# `w` and `fn` shadow the constructors inside mm()'s body)
.default_w <- function() w(~ 1 / n, scale = TRUE)
.default_fn <- function() fn("sum")

#' Define a multiple-membership block (the micro-to-macro aggregation)
#'
#' @description
#' Specifies a multiple-membership level where group-level units (e.g.,
#' occupations) are composed of member-level units (e.g., tasks). Each block
#' aggregates the weighted member records into one or more \emph{group-level
#' features} whose coefficients are estimated by the main model. The components
#' map onto the framework's notation \eqn{\theta^{micro,f}(M_{it})}:
#'
#' \itemize{
#'   \item \code{id = id(mmid, mainid)}: the membership structure \eqn{S_{it}}
#'   \item \code{vars = vars(x)}: the member attributes \eqn{x_{kt}}
#'   \item \code{w = w(~ ...)}: the weights \eqn{w_{ikt}} (who counts)
#'   \item \code{fn = fn("...")}: the aggregation function \eqn{f} (how
#'     contributions combine); \code{fn("sum")} is the additive case, other
#'     types are emergent features (variance, concentration, thresholds, ...)
#'   \item \code{RE = re(...)}: member random effects, aggregated by the weights
#'     (\eqn{\sum_k w_{ikt} u_{0k}}, plus random slopes via \code{re(1 + x)})
#'   \item \code{name}: block name for referencing the feature in main-formula
#'     interactions (auto-generated when omitted)
#' }
#'
#' @param id An \code{\link{id}} object: \code{id(mmid, mainid)}.
#' @param vars A \code{\link{vars}} object with member attributes, or
#'   \code{NULL} for weights-only blocks (RE-only, or \code{fn("hhi")}-family).
#' @param w A \code{\link{w}} object specifying the weights. Default:
#'   \code{w(~ 1/n, scale = TRUE)} (equal weights).
#' @param fn A \code{\link{fn}} object selecting the aggregation function.
#'   Default: \code{fn("sum")} (the additive weighted mean).
#' @param RE Random effects: \code{TRUE} (shorthand for \code{re(1)}), a
#'   \code{\link{re}} object, or \code{NULL}/\code{FALSE} for none. Only
#'   available with \code{fn("sum")}: member random effects do not compose with
#'   dispersion-type features.
#' @param FE Fixed effects: a \code{\link{fe}} object, or \code{NULL}. Mutually
#'   exclusive with \code{RE}. With many members this is weakly identified —
#'   prefer \code{RE} (partial pooling) at the member level.
#' @param name Optional block name (unquoted or string) used to reference the
#'   block's feature in main-formula interactions (e.g. \code{Ax:education}).
#'   Auto-generated from the feature when omitted (e.g. \code{A_x}, \code{V_x}).
#'   Must not collide with a data column.
#' @param ... Not used; catches removed arguments (\code{c =}, and \code{ar =},
#'   which moved into the effects grammar: \code{RE = re(1, ar = TRUE)}) with a
#'   migration message.
#'
#' @return A \code{bml_mm} object.
#'
#' @details
#' \strong{One mechanism.} Every block emits group-level feature(s) whose
#' coefficients are main-model coefficients (class \code{"b"}, labeled by the
#' feature name, e.g. \code{A_x}). \code{fn("sum")} is the additive case
#' \eqn{\theta = \beta A_x + \sum_k w_k u_{0k}}: the feature \eqn{A_x} plus the
#' optional \code{RE}. There is no separate "sum mode".
#'
#' \strong{Multiple blocks} can be combined with \code{+} to stack features
#' (mean + variance + concentration, ...). \code{RE} can be specified for one
#' block per member-id group.
#'
#' @seealso \code{\link{bml}}, \code{\link{id}}, \code{\link{vars}},
#'   \code{\link{w}}, \code{\link{fn}}, \code{\link{re}}, \code{\link{fe}},
#'   \code{\link{hm}}
#'
#' @examples
#' \dontrun{
#' # Additive aggregation (weighted mean effect)
#' mm(id = id(task, occ), vars = vars(x), w = w(~ importance, scale = TRUE))
#'
#' # With member random intercepts
#' mm(id = id(task, occ), vars = vars(x), w = w(~ 1/n), RE = TRUE)
#'
#' # Random intercept + slope (residual effect heterogeneity)
#' mm(id = id(task, occ), vars = vars(x), w = w(~ 1/n), RE = re(1 + x))
#'
#' # Emergent features
#' mm(id = id(task, occ), vars = vars(x), w = w(~ 1/n), fn = fn("var"))
#' mm(id = id(task, occ), w = w(~ importance), fn = fn("hhi"))
#' mm(id = id(task, occ), vars = vars(x), fn = fn("threshold", c = est(), kappa = 10))
#' mm(id = id(task, occ), vars = vars(x), fn = fn("smax", kappa = est()))
#'
#' # Named block, referenced in a cross-level interaction
#' bml(Y ~ education + Ax:education +
#'       mm(name = Ax, id = id(task, occ), vars = vars(x), w = w(~ 1/n)),
#'     data = dat)
#' }
#'
#' @references
#' Rosche, B. (2026). A Multilevel Model for Coalition Governments: Uncovering
#' Party-Level Dependencies Within and Between Governments. \emph{Political Analysis}.
#'
#' @export
mm <- function(id, vars = NULL, w = NULL, fn = NULL, RE = NULL, FE = NULL,
               name = NULL, ...) {

  dots <- list(...)
  if ("c" %in% names(dots)) {
    stop("The 'c' argument moved: normalization is now w(..., scale = TRUE/FALSE).", call. = FALSE)
  }
  if ("ar" %in% names(dots)) {
    stop("The 'ar' argument moved into the effects grammar:\n",
         "  mm(..., RE = TRUE, ar = TRUE)  ->  mm(..., RE = re(1, ar = TRUE))\n",
         "  For a calendar-time walk: RE = re(1, ar = <time column>).", call. = FALSE)
  }
  if (length(dots) > 0) {
    stop("Unknown mm() argument(s): ", paste(names(dots), collapse = ", "), call. = FALSE)
  }

  # Validate id
  if (!inherits(id, "bml_id")) {
    stop("'id' must be specified using id(mmid, mainid)", call. = FALSE)
  }
  if (length(id) != 2) {
    stop("'id' must have exactly 2 identifiers: id(mmid, mainid)", call. = FALSE)
  }

  # Weights: default equal weights; catch the old fn(w ~ ...) spelling
  if (is.null(w)) {
    w <- .default_w()
  }
  if (!inherits(w, "bml_w")) {
    stop("'w' must be specified using w(~ ...), e.g. w(~ 1/n) or w(~ importance, scale = TRUE).",
         call. = FALSE)
  }

  # Aggregation function: default sum
  if (is.null(fn)) {
    fn <- .default_fn()
  }
  if (!inherits(fn, "bml_fn")) {
    stop("'fn' must be specified using fn(), e.g. fn(\"sum\") or fn(\"var\").\n",
         "  Note: weights moved to w(); the old fn(w ~ ...) spelling was removed.", call. = FALSE)
  }

  # Effects grammar
  if (isTRUE(RE)) RE <- re(1)
  if (isFALSE(RE)) RE <- NULL
  if (!is.null(RE) && !inherits(RE, "bml_re")) {
    stop("'RE' must be TRUE, FALSE, or an re() object (e.g. RE = re(1 + x)).", call. = FALSE)
  }
  if (!is.null(FE) && !inherits(FE, "bml_fe")) {
    stop("'FE' must be a fe() object (e.g. FE = fe(1)).", call. = FALSE)
  }
  if (!is.null(RE) && !is.null(FE)) {
    stop("mm(): specify either RE = re(...) or FE = fe(...), not both.", call. = FALSE)
  }

  # Block name (unquoted symbol or string) - never force the promise directly:
  # bare symbols like name = Ax are taken as literal names
  name_expr <- substitute(name)
  block_name <- if (is.null(name_expr)) {
    NULL
  } else {
    val <- tryCatch(eval(name_expr, parent.frame()), error = function(e) NULL)
    if (is.character(val) && length(val) == 1) val else deparse(name_expr)
  }

  # fn-type validations that don't need the data
  weights_only <- fn$type %in% c("hhi", "effn", "entropy")
  if (weights_only && !is.null(vars)) {
    stop("fn(\"", fn$type, "\") is a function of the weights alone and takes no vars(). ",
         "Remove vars= or choose an attribute-based fn type.", call. = FALSE)
  }
  if (!weights_only && fn$type != "sum" && is.null(vars)) {
    stop("fn(\"", fn$type, "\") requires member attributes: specify vars().", call. = FALSE)
  }

  # RE/FE composition rules
  if (fn$type != "sum" && !is.null(RE)) {
    stop("RE is only available with fn(\"sum\"): member random effects do not compose with ",
         "dispersion/emergent features (the variance of x'b + u_k has no feature form). ",
         "Put the RE in a separate fn(\"sum\") block.", call. = FALSE)
  }
  if (fn$type != "sum" && !is.null(FE)) {
    stop("FE is only available with fn(\"sum\").", call. = FALSE)
  }
  if (!is.null(FE)) {
    warning("mm(FE = fe(...)): with many members, member-level fixed effects are weakly ",
            "identified (one parameter per member, collinear with the intercept and A_x). ",
            "Consider RE = re(...) (partial pooling) instead.", call. = FALSE)
  }

  # Emergent types require a normalized measure
  if (fn$type != "sum" && !isTRUE(w$scale)) {
    warning("fn(\"", fn$type, "\") is only interpretable for normalized weights; ",
            "forcing w(scale = TRUE).", call. = FALSE)
    w$scale <- TRUE
    w$constraint <- TRUE
  }

  # RE default: weights-only sum blocks are RE blocks (variance decomposition)
  if (is.null(vars) && fn$type == "sum" && is.null(FE)) {
    if (is.null(RE)) RE <- re(1)
  }

  structure(
    list(
      id   = as.character(id),
      vars = vars,
      w    = w,
      fn   = fn,
      RE   = RE,
      FE   = FE,
      name = block_name,
      # the ar spec lives on re(); carried here for the downstream pipeline
      ar   = if (!is.null(RE)) RE$ar else NULL
    ),
    class = "bml_mm"
  )
}

#' Define a hierarchical nesting structure
#'
#' @description
#' Specifies a hierarchical (nesting) level: each group-level unit belongs to
#' exactly one nesting-level unit (e.g., governments within countries).
#' \code{hm()} carries only the effect \emph{structure} — random effects
#' (\code{RE = re(...)}, partial pooling) or fixed effects (\code{FE = fe(...)},
#' no pooling). Fixed effects of unit-level \emph{covariates} belong in the main
#' formula, following the lme4/brms convention:
#'
#' \preformatted{
#' Y ~ gdp + hm(id = id(cid), RE = re(1 + gdp))   # like y ~ gdp + (1 + gdp | cid)
#' }
#'
#' @param id An \code{\link{id}} object: \code{id(hmid)}.
#' @param RE Random effects: \code{TRUE} (shorthand for \code{re(1)}), an
#'   \code{\link{re}} object, or \code{NULL}. Slope variables are main-level
#'   columns varying within the nesting units. Default when neither \code{RE}
#'   nor \code{FE} is given: \code{re(1)}.
#' @param FE Fixed effects: a \code{\link{fe}} object (unit dummies via
#'   \code{fe(1)}; unit-specific slopes via \code{fe(1 + x)}; report per-unit
#'   estimates with \code{fe(1, showFE = TRUE)}). Mutually exclusive with
#'   \code{RE}.
#' @param labels Unquoted variable name holding display labels for the nesting
#'   units (used when reporting per-unit fixed effects).
#' @param ... Not used; catches the removed \code{type =}, \code{vars =},
#'   \code{name =}, \code{showFE =}, and \code{ar =} arguments (the last moved
#'   into the effects grammar: \code{RE = re(1, ar = TRUE)}) with a migration
#'   message.
#'
#' @return A \code{bml_hm} object.
#'
#' @details
#' Cross-classified structures are modeled by including multiple \code{hm()}
#' blocks. The old \code{type = "RE"/"FE"} and \code{vars =} arguments were
#' removed: use \code{RE = re(...)} / \code{FE = fe(...)}, and move unit-level
#' covariates into the main formula (their coefficients are ordinary main-model
#' coefficients).
#'
#' @seealso \code{\link{bml}}, \code{\link{mm}}, \code{\link{id}},
#'   \code{\link{re}}, \code{\link{fe}}
#'
#' @examples
#' \dontrun{
#' hm(id = id(cid))                          # random intercepts (default)
#' hm(id = id(cid), RE = re(1 + x))          # random intercept + slope on x
#' hm(id = id(cid), FE = fe(1))              # country dummies (no pooling)
#' hm(id = id(cid), FE = fe(1, showFE = TRUE), labels = cname)
#' }
#'
#' @references
#' Goldstein, H. (2011). \emph{Multilevel Statistical Models} (4th ed.). Wiley.
#'
#' @export
hm <- function(id, RE = NULL, FE = NULL, labels = NULL, ...) {

  dots <- list(...)
  if (any(c("type", "vars", "name", "showFE") %in% names(dots))) {
    stop("hm() was redesigned:\n",
         "  - type = \"RE\"/\"FE\"  ->  RE = re(...) / FE = fe(...)\n",
         "  - vars = vars(gdp)   ->  put gdp in the main formula (Y ~ gdp + hm(...))\n",
         "  - showFE = TRUE      ->  FE = fe(1, showFE = TRUE)\n",
         "  - name = cname       ->  labels = cname", call. = FALSE)
  }
  if ("ar" %in% names(dots)) {
    stop("The 'ar' argument moved into the effects grammar:\n",
         "  hm(..., ar = TRUE)  ->  hm(..., RE = re(1, ar = TRUE))\n",
         "  For a calendar-time walk: RE = re(1, ar = <time column>).", call. = FALSE)
  }
  if (length(dots) > 0) {
    stop("Unknown hm() argument(s): ", paste(names(dots), collapse = ", "), call. = FALSE)
  }

  # Validate id
  if (!inherits(id, "bml_id")) {
    stop("'id' must be specified using id(hmid)", call. = FALSE)
  }
  if (length(id) != 1) {
    stop("'id' must have exactly 1 identifier for hm(): id(hmid)", call. = FALSE)
  }

  if (isTRUE(RE)) RE <- re(1)
  if (isFALSE(RE)) RE <- NULL
  if (!is.null(RE) && !inherits(RE, "bml_re")) {
    stop("'RE' must be TRUE, FALSE, or an re() object (e.g. RE = re(1 + x)).", call. = FALSE)
  }
  if (!is.null(FE) && !inherits(FE, "bml_fe")) {
    stop("'FE' must be a fe() object (e.g. FE = fe(1)).", call. = FALSE)
  }
  if (!is.null(RE) && !is.null(FE)) {
    stop("hm(): specify either RE = re(...) or FE = fe(...), not both.", call. = FALSE)
  }
  if (is.null(RE) && is.null(FE)) {
    RE <- re(1)
  }

  labels_expr <- substitute(labels)
  labels_char <- if (is.null(labels_expr)) {
    NULL
  } else {
    val <- tryCatch(eval(labels_expr, parent.frame()), error = function(e) NULL)
    if (is.character(val) && length(val) == 1) val else deparse(labels_expr)
  }

  structure(
    list(
      id     = as.character(id),
      RE     = RE,
      FE     = FE,
      labels = labels_char,
      # the ar spec lives on re(); carried here for the downstream pipeline
      ar     = if (!is.null(RE)) RE$ar else NULL
    ),
    class = "bml_hm"
  )
}
