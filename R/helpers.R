# ================================================================================================ #
# Helper functions for bml formula syntax
# ================================================================================================ #

#' Specify identifier variables for multiple-membership and hierarchical structures
#'
#' @description
#' Helper function used within \code{\link{mm}} and \code{\link{hm}} to specify
#' the identifier variables that define memberships and nesting structures. In
#' multiple-membership models, \code{id()} links member-level units (e.g., party
#' IDs) to group-level units (e.g., government IDs). In hierarchical models,
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
#' @seealso \code{\link{vars}}, \code{\link{mm}}, \code{\link{hm}}
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

#' Specify covariates for multiple-membership or hierarchical models
#'
#' @description
#' Helper function used within \code{\link{mm}} and \code{\link{hm}} to specify
#' which variables should be included at each level of the model. Supports both
#' free variables (with coefficients to be estimated) and fixed variables (with
#' coefficients held constant using \code{\link{fix}}).
#'
#' @param ... Unquoted variable names from your data, combined using \code{+}
#'   (formula-style). Supports:
#'   \itemize{
#'     \item Simple variables: \code{vars(x + y)}
#'     \item Interactions: \code{vars(x * y)} or \code{vars(x:y)}
#'     \item Transformations: \code{vars(I(x^2))} or \code{vars(I(x + y))}
#'     \item Fixed coefficients: \code{vars(fix(x, 1.0) + y)}
#'   }
#'   Note: Numeric literals like \code{1}, \code{0}, or \code{-1} are ignored
#'   (no intercept support in mm/hm blocks).
#'
#' @return A \code{bml_vars} object containing:
#'   \itemize{
#'     \item \code{formula}: Formula object for use with \code{model.matrix()}
#'     \item \code{free}: Character vector of base variable names
#'     \item \code{fixed}: List of variables with fixed coefficients (if any)
#'   }
#'   Returns \code{NULL} if no variables are specified.
#'
#' @seealso \code{\link{fix}}, \code{\link{mm}}, \code{\link{hm}}
#'
#' @examples
#' \dontrun{
#' # Simple variable specification (formula-style with +)
#' vars(income + education)
#'
#' # Single variable
#' vars(income)
#'
#' # Interactions
#' vars(income * education)  # expands to income + education + income:education
#' vars(income:education)    # interaction only
#'
#' # Transformations
#' vars(I(income^2))         # squared term
#' vars(income + I(income^2)) # linear and squared
#'
#' # Mix free and fixed variables
#' vars(fix(exposure, 1.0) + income + education)
#'
#' # Use in mm() specification
#' mm(
#'   id = id(pid, gid),
#'   vars = vars(rile + ipd),
#'   fn = fn(w ~ 1/n),
#'   RE = FALSE
#' )
#' }
#'
#' @export
vars <- function(...) {
  expr <- substitute(list(...))

  # Check for comma-separated arguments (should use formula-style with +)
  if (length(expr) > 2) {
    stop("Variables in vars() must be specified using formula-style with '+', not comma-separated.\n",
         "  Correct:   vars(x + y + z)\n",
         "  Incorrect: vars(x, y, z)")
  }

  fixed_list <- list()
  free_terms <- character()  # Store terms (not just var names) for formula

  # Recursive function to extract fix() calls and build formula terms
  extract_and_build <- function(e, depth = 0) {
    if (is.call(e)) {
      fn_name <- as.character(e[[1]])

      if (fn_name == "fix") {
        # It's a fix() call - extract var name and value
        var_name <- as.character(e[[2]])
        value <- eval(e[[3]], parent.frame(2))
        fixed_list[[length(fixed_list) + 1]] <<- list(var = var_name, value = value)
        return(NULL)  # Don't include in formula
      } else if (fn_name == "+") {
        # Addition - process both sides
        left <- extract_and_build(e[[2]], depth + 1)
        right <- extract_and_build(e[[3]], depth + 1)
        # Combine non-NULL results
        return(c(left, right))
      } else if (fn_name == "list" && depth == 0) {
        # Top-level list() wrapper
        if (length(e) >= 2) {
          return(extract_and_build(e[[2]], depth + 1))
        }
        return(NULL)
      } else {
        # Other calls (I(), :, *, etc.) - keep as term
        term_str <- deparse(e, width.cutoff = 500)
        term_str <- paste(term_str, collapse = "")
        return(term_str)
      }
    } else if (is.name(e) && as.character(e) != "") {
      var_name <- as.character(e)
      if (var_name != "list") {
        return(var_name)
      }
    } else if (is.numeric(e)) {
      # Numeric literals (1, 0, -1) - ignore for intercept
      return(NULL)
    }
    return(NULL)
  }

  free_terms <- extract_and_build(expr)
  free_terms <- free_terms[!is.null(free_terms) & free_terms != ""]
  free_terms <- unique(free_terms)

  if (length(free_terms) == 0 && length(fixed_list) == 0) return(NULL)

  # Create formula for model.matrix() (no intercept in mm/hm)
  if (length(free_terms) > 0) {
    formula_str <- paste("~ 0 +", paste(free_terms, collapse = " + "))
    vars_formula <- as.formula(formula_str)
    # Extract base variable names using all.vars()
    base_vars <- all.vars(vars_formula)
  } else {
    vars_formula <- NULL
    base_vars <- character()
  }

  # Return structure with formula
  structure(
    list(
      formula = vars_formula,
      free = base_vars,
      fixed = if (length(fixed_list) > 0) fixed_list else NULL
    ),
    class = "bml_vars"
  )
}

#' Specify a weight function for multiple-membership models
#'
#' @description
#' Defines how member-level contributions are weighted when aggregating to the
#' group level (the "micro-macro link"). The weight function can be a simple
#' formula (e.g., \code{1/n} for equal weights) or can include parameters to be
#' estimated from the data.
#'
#' @param w A two-sided formula specifying the weight function. The left-hand side
#'   must be \code{w}; the right-hand side defines the weighting scheme:
#'   \itemize{
#'     \item Simple: \code{w ~ 1/n} (equal weights based on group size)
#'     \item Parameterized: \code{w ~ b0 + b1 * tenure} (weights depend on
#'       member characteristics and estimated parameters)
#'     \item With group aggregates: \code{w ~ b1 * min(x) + (1-b1) * mean(x)}
#'       (weights based on group-level summaries; see Details)
#'   }
#'   Parameters must be named \code{b0}, \code{b1}, \code{b2}, etc.
#'
#' @param c Logical; if \code{TRUE} (default), weights are normalized to sum to 1
#'   within each group. Set to \code{FALSE} for unnormalized weights.
#'
#' @return A \code{bml_fn} object containing the parsed weight function specification.
#'
#' @details
#' \strong{Weight Function Components:}
#' \itemize{
#'   \item \strong{Variables} (e.g., \code{n}, \code{tenure}): Data from your dataset
#'   \item \strong{Parameters} (e.g., \code{b0}, \code{b1}): Estimated from the data
#'   \item \strong{Operations}: Standard R arithmetic (\code{+}, \code{-}, \code{*},
#'     \code{/}, \code{^}, etc.)
#' }
#'
#' \strong{Common Weight Functions:}
#' \itemize{
#'   \item Equal weights: \code{w ~ 1/n}
#'   \item Duration-based: \code{w ~ duration}
#'   \item Flexible parameterized: \code{w ~ b0 + b1 * seniority}
#'   \item Group aggregates: \code{w ~ b1 * min(x) + (1-b1) * mean(x)}
#' }
#'
#' When \code{c = TRUE}, the weights are constrained: \eqn{\sum_{k \in group} w_k = 1}.
#'
#' \strong{Group-Level Aggregation Functions:}
#'
#' The weight function supports aggregation functions that compute summaries
#' within each group (mainid). These are pre-computed in R before passing to JAGS.
#' Supported functions:
#' \itemize{
#'   \item \code{min(var)}, \code{max(var)}: Minimum/maximum value within the group
#'   \item \code{mean(var)}, \code{sum(var)}: Mean/sum of values within the group
#'   \item \code{median(var)}, \code{mode(var)}: Median/mode (most frequent) value within the group
#'   \item \code{sd(var)}, \code{var(var)}, \code{range(var)}: Standard deviation/variance/range (max-min) within the group
#'   \item \code{first(var)}, \code{last(var)}: First/last value (based on data order)
#'   \item \code{quantile(var, prob)}: Quantile at probability \code{prob} (0 to 1).
#'     For example, \code{quantile(x, 0.25)} computes the 25th percentile.
#' }
#'
#' Example: \code{fn(w ~ b1 * min(tenure) + (1-b1) * max(tenure))} creates weights
#' that blend the minimum and maximum tenure within each group, with the blend
#' controlled by the estimated parameter \code{b1}.
#'
#' Example with quantile: \code{fn(w ~ quantile(tenure, 0.75) / max(tenure))} uses
#' the 75th percentile relative to the maximum within each group.
#'
#' Note: Nested aggregation functions (e.g., \code{min(max(x))}) are not supported.
#'
#' \strong{JAGS Mathematical Functions:}
#'
#' The following mathematical functions are passed directly to JAGS and can be
#' used in weight formulas:
#' \itemize{
#'   \item \code{exp}, \code{log}, \code{log10}, \code{sqrt}, \code{abs}, \code{pow}
#'   \item \code{sin}, \code{cos}, \code{tan}, \code{asin}, \code{acos}, \code{atan}
#'   \item \code{sinh}, \code{cosh}, \code{tanh}
#'   \item \code{logit}, \code{ilogit}, \code{probit}, \code{iprobit}, \code{cloglog}, \code{icloglog}
#'   \item \code{round}, \code{trunc}, \code{floor}, \code{ceiling}
#' }
#'
#' Example: \code{fn(w ~ 1 / (1 + (n - 1) * exp(-(b1 * x))))} uses an exponential
#' decay function where weights depend on member characteristics.
#'
#' \strong{Ensuring Numerical Stability:}
#'
#' Weight functions with estimated parameters (\code{b0}, \code{b1}, ...) must
#' produce bounded, positive values across all plausible parameter values.
#' Unbounded weight functions can cause the MCMC sampler to crash
#' (e.g., \code{"Error in node w.1[...]: Invalid parent values"}).
#' During sampling, weight parameters can take on extreme values, and if the
#' weight function is not bounded, this will destabilize the likelihood.
#'
#' Recommendations:
#' \itemize{
#'   \item \strong{Use bounded weight functions.} Two options:
#'     \itemize{
#'       \item \code{ilogit()}: Bounds weights between 0 and 1 with a zero-point
#'         at 0.5: \code{fn(w ~ ilogit(b0 + b1 * x), c = TRUE)}
#'       \item \strong{Generalized logistic} (Rosche, 2026): Bounds weights between
#'         0 and 1 with a zero-point at \eqn{1/n} (equal weights), so deviations
#'         from equal weighting are estimated as a function of covariates:
#'         \code{fn(w ~ 1 / (1 + (n - 1) * exp(-(b0 + b1 * x))), c = TRUE)}
#'     }
#'   \item \strong{Use \code{c = TRUE}} (weight normalization) to prevent weights
#'     from growing without bound
#'   \item \strong{Standardize covariates} in the weight function. Variables with
#'     large ranges (e.g., income in thousands) can cause \code{b * x} to overflow
#'   \item \strong{Use informative priors} for weight parameters via the
#'     \code{priors} argument in \code{\link{bml}} (e.g.,
#'     \code{priors = list("b.w.1[1] ~ dnorm(0, 1)")})
#'   \item \strong{Avoid unbounded functions} like \code{exp(b * x)} without
#'     normalization (\code{c = TRUE}) or wrapping (e.g., inside \code{ilogit()})
#' }
#'
#' Weight parameters are initialized at 0 by default to ensure numerically stable
#' starting values. See \code{vignette("faq")} (Question 7) for detailed
#' troubleshooting of numerical issues.
#'
#' @seealso \code{\link{mm}}, \code{\link{bml}}, \code{vignette("model")} for
#'   the model structure, \code{vignette("faq")} for troubleshooting
#'
#' @examples
#' \dontrun{
#' # Equal weights (standard multiple-membership)
#' fn(w ~ 1/n, c = TRUE)
#'
#' # Tenure-based weights (proportional to time served)
#' fn(w ~ tenure, c = TRUE)
#'
#' # Flexible parameterized weights
#' fn(w ~ b0 + b1 * seniority, c = TRUE)
#'
#' # Unconstrained weights
#' fn(w ~ importance, c = FALSE)
#'
#' # Weights based on group aggregates
#' fn(w ~ b1 * min(tenure) + (1 - b1) * mean(tenure), c = TRUE)
#'
#' # Combining individual and aggregate measures
#' fn(w ~ b0 + b1 * (tenure / max(tenure)), c = TRUE)
#'
#' # Using median for robust central tendency
#' fn(w ~ tenure / median(tenure), c = TRUE)
#'
#' # Using quantiles for percentile-based weights
#' fn(w ~ quantile(tenure, 0.75) - quantile(tenure, 0.25), c = TRUE)
#' }
#'
#' @references
#' Rosche, B. (2026). A Multilevel Model for Theorizing and Estimating the
#' Micro-Macro Link. \emph{Political Analysis}.
#'
#' Browne, W. J., Goldstein, H., & Rasbash, J. (2001). Multiple membership
#' multiple classification (MMMC) models. \emph{Statistical Modelling}, 1(2), 103-124.
#'
#' @export
fn <- function(w = w ~ 1/n, c = TRUE) {

  # Convert to formula if needed
  if (!inherits(w, "formula")) {
    stop("Weight function 'w' must be a formula (e.g., w ~ 1/n)")
  }

  # Extract components from formula
  fn_string <- deparse(w[[3]], width.cutoff = 500)
  fn_string <- paste(fn_string, collapse = "")

  # Insert intercept variable for b0
  fn_string <- gsub("\\bb0\\b", "b0 * X0", fn_string)

  # ========================================================================================== #
  # Detect aggregation functions: min, max, mean, sum, sd, var, first, last, median, mode, range, quantile
  # ========================================================================================== #

  # Single-argument aggregation functions (computed in R before passing to JAGS)
  supported_agg_simple <- c("min", "max", "mean", "sum", "sd", "var", "first", "last", "median", "mode", "range")
  agg_pattern_simple <- "\\b(min|max|mean|sum|sd|var|first|last|median|mode|range)\\s*\\(\\s*([a-zA-Z_][a-zA-Z0-9_\\.]*)\\s*\\)"

  # Two-argument aggregation function: quantile(var, prob)
  agg_pattern_quantile <- "\\bquantile\\s*\\(\\s*([a-zA-Z_][a-zA-Z0-9_\\.]*)\\s*,\\s*(-?[0-9]*\\.?[0-9]+)\\s*\\)"

  # All supported aggregation functions
  supported_agg <- c(supported_agg_simple, "quantile")

  # JAGS mathematical functions (passed through directly to JAGS)
  jags_math_funcs <- c("exp", "log", "log10", "sqrt", "abs", "pow",
                       "sin", "cos", "tan", "asin", "acos", "atan",
                       "sinh", "cosh", "tanh",
                       "logit", "ilogit", "probit", "iprobit",
                       "cloglog", "icloglog",
                       "round", "trunc", "floor", "ceiling")

  # All allowed functions (aggregation + JAGS math)
  all_allowed_funcs <- c(supported_agg, jags_math_funcs)

  # Check for nested aggregation functions (not supported)
  nested_pattern <- "\\b(min|max|mean|sum|sd|var|first|last|median|mode|range|quantile)\\s*\\([^)]*\\b(min|max|mean|sum|sd|var|first|last|median|mode|range|quantile)\\s*\\("
  if (grepl(nested_pattern, fn_string)) {
    stop("Nested aggregation functions are not supported in weight formulas.\n",
         "  Example of invalid: fn(w ~ min(max(x)))\n",
         "  Valid: fn(w ~ min(x) + max(y))", call. = FALSE)
  }

  # Check for unsupported functions (any function call not in allowed list)
  all_func_calls <- stringr::str_match_all(fn_string, "\\b([a-zA-Z_][a-zA-Z0-9_]*)\\s*\\(")[[1]]
  if (nrow(all_func_calls) > 0) {
    func_names <- all_func_calls[, 2]
    unsupported <- setdiff(func_names, all_allowed_funcs)
    if (length(unsupported) > 0) {
      stop("Unsupported function(s) in weight formula: ", paste(unsupported, collapse = ", "), "\n",
           "  Supported aggregation functions: ", paste(supported_agg, collapse = ", "), "\n",
           "  Supported JAGS math functions: ", paste(jags_math_funcs, collapse = ", "), "\n",
           "  Note: Use variables directly (e.g., tenure) or with parameters (e.g., b1 * tenure)", call. = FALSE)
    }
  }

  # Extract simple aggregation function calls (single argument)
  agg_matches_simple <- stringr::str_match_all(fn_string, agg_pattern_simple)[[1]]

  # Extract quantile function calls (two arguments)
  agg_matches_quantile <- stringr::str_match_all(fn_string, agg_pattern_quantile)[[1]]

  agg_funcs <- NULL
  agg_vars <- NULL

  # Process simple aggregation functions
  if (nrow(agg_matches_simple) > 0) {
    agg_funcs <- list()
    for (i in 1:nrow(agg_matches_simple)) {
      agg_funcs[[length(agg_funcs) + 1]] <- list(
        original = agg_matches_simple[i, 1],  # e.g., "min(fdp)"
        func     = agg_matches_simple[i, 2],  # e.g., "min"
        var      = agg_matches_simple[i, 3],  # e.g., "fdp"
        prob     = NULL,                       # Not applicable for simple functions
        col_name = paste0(agg_matches_simple[i, 3], "_", agg_matches_simple[i, 2])  # e.g., "fdp_min"
      )
    }
  }

  # Process quantile function calls
  if (nrow(agg_matches_quantile) > 0) {
    if (is.null(agg_funcs)) agg_funcs <- list()
    for (i in 1:nrow(agg_matches_quantile)) {
      var_name <- agg_matches_quantile[i, 2]
      prob_val <- as.numeric(agg_matches_quantile[i, 3])

      # Validate probability is between 0 and 1
      if (prob_val < 0 || prob_val > 1) {
        stop("quantile() probability must be between 0 and 1. Got: ", prob_val, call. = FALSE)
      }

      # Create column name: var_q25 for quantile(var, 0.25)
      prob_pct <- round(prob_val * 100)
      col_name <- paste0(var_name, "_q", prob_pct)

      agg_funcs[[length(agg_funcs) + 1]] <- list(
        original = agg_matches_quantile[i, 1],  # e.g., "quantile(fdp, 0.25)"
        func     = "quantile",
        var      = var_name,
        prob     = prob_val,                     # Store the probability
        col_name = col_name                      # e.g., "fdp_q25"
      )
    }
  }

  # Collect base variables needed for computing aggregates
  if (!is.null(agg_funcs)) {
    agg_vars <- unique(sapply(agg_funcs, function(x) x$var))
  }

  # Transform formula string: replace aggregation calls with column names
  fn_string_transformed <- fn_string
  if (!is.null(agg_funcs)) {
    for (agg in agg_funcs) {
      # Escape special regex characters in the original match
      escaped <- gsub("([\\(\\)\\.])", "\\\\\\1", agg$original)
      fn_string_transformed <- gsub(escaped, agg$col_name, fn_string_transformed)
    }
  }

  # ========================================================================================== #
  # Extract variables and parameters from the transformed string
  # ========================================================================================== #

  # Extract parameters (b0, b1, b2, ...)
  params <- unique(stringr::str_extract_all(fn_string_transformed, "\\bb\\d+\\b")[[1]])

  # Get all variables from the transformed formula string
  # Parse the transformed string to extract variable names
  all_v_transformed <- all.vars(as.formula(paste("w ~", fn_string_transformed)))
  wvars <- setdiff(all_v_transformed, c("w", params))

  # Variables that have parameters (b * X pattern)
  vars_p <- stringr::str_match_all(fn_string_transformed, "b\\d+\\s*\\*\\s*([a-zA-Z_][a-zA-Z0-9_\\.]*)")[[1]]
  vars_p <- if (nrow(vars_p) > 0) vars_p[, 2] else character(0)

  structure(
    list(
      formula    = w,
      string     = fn_string_transformed,
      vars       = wvars,
      vars_p     = vars_p,
      params     = params,
      constraint = c,
      agg_funcs  = agg_funcs,
      agg_vars   = agg_vars
    ),
    class = "bml_fn"
  )
}

#' Define a multiple-membership structure
#'
#' @description
#' Specifies a multiple-membership level in the model where group-level units
#' (e.g., governments) are composed of multiple member-level units (e.g., political
#' parties). Unlike pure hierarchical nesting, members can belong to multiple groups,
#' and their contributions are aggregated using a user-specified weight function.
#'
#' @param id An \code{\link{id}} object specifying the member-level and group-level
#'   identifiers: \code{id(mmid, mainid)} where \code{mmid} identifies members
#'   and \code{mainid} identifies groups.
#'
#' @param vars A \code{\link{vars}} object specifying member-level covariates to
#'   aggregate, or \code{NULL} for random effects only. Supports interactions
#'   (\code{*}, \code{:}) and transformations (\code{I()}). Variables are weighted
#'   according to the function specified in \code{fn}.
#'
#' @param fn A \code{\link{fn}} object specifying the weight function (default:
#'   \code{fn(w ~ 1/n, c = TRUE)} for equal weights). Note: Weight functions do
#'   NOT support interactions or \code{I()} - pre-create any needed transformed
#'   variables in your data. See \code{\link{fn}} for details.
#'
#' @param RE Logical; if \code{TRUE}, include random effects for member-level units.
#'   Automatically set to \code{TRUE} if \code{vars = NULL} (random effects only).
#'
#' @param ar Logical; if \code{TRUE}, random effects evolve autoregressively across
#'   participations. Requires members to have sequential participation indicators
#'   in the data. Default: \code{FALSE}.
#'
#' @return A \code{bml_mm} object containing the multiple-membership specification.
#'
#' @details
#' \strong{Multiple-Membership Models:}
#'
#' In standard hierarchical models, each observation belongs to exactly one group.
#' Multiple-membership models relax this assumption, allowing groups to be composed
#' of multiple members, with flexible weighting of member contributions.
#'
#' \strong{Model Structure:}
#'
#' The contribution from mm block \eqn{k} to group \eqn{j} is:
#'
#' \deqn{\mathrm{mm}_{kj} = \sum_{i \in group_j} w_{ki} (x_{ki}' \beta_k + \alpha_{ki})}
#'
#' where:
#' \itemize{
#'   \item \eqn{w_{ki}}: Weight for member \eqn{i} in group \eqn{j} (from \code{fn})
#'   \item \eqn{x_{ki}}: Member-level covariates (from \code{vars})
#'   \item \eqn{\beta_k}: Regression coefficients (estimated)
#'   \item \eqn{\alpha_{ki}}: Member-level random effect (if \code{RE = TRUE})
#' }
#'
#' \strong{Multiple mm() Blocks:}
#'
#' You can specify multiple \code{mm()} blocks with different weight functions,
#' variables, or random effect specifications. This allows modeling different
#' aggregation mechanisms simultaneously.
#'
#' @seealso \code{\link{bml}}, \code{\link{id}}, \code{\link{vars}}, \code{\link{fn}},
#'   \code{\link{hm}}
#'
#' @examples
#' \dontrun{
#' # Equal weights with variables
#' mm(
#'   id = id(pid, gid),
#'   vars = vars(rile + ipd),
#'   fn = fn(w ~ 1/n, c = TRUE),
#'   RE = FALSE
#' )
#'
#' # Random effects only (no variables)
#' mm(
#'   id = id(pid, gid),
#'   vars = NULL,
#'   fn = fn(w ~ 1/n, c = TRUE),
#'   RE = TRUE  # Automatically TRUE when vars = NULL
#' )
#'
#' # Flexible weights with parameter
#' mm(
#'   id = id(pid, gid),
#'   vars = vars(org_structure),
#'   fn = fn(w ~ ilogit(b0 + b1 * pseat), c = TRUE),
#'   RE = TRUE
#' )
#'
#' # Autoregressive random effects
#' mm(
#'   id = id(pid, gid),
#'   vars = NULL,
#'   fn = fn(w ~ 1/n, c = TRUE),
#'   RE = TRUE,
#'   ar = TRUE  # Random effects evolve over participations
#' )
#'
#' # Interactions and transformations in vars
#' mm(
#'   id = id(pid, gid),
#'   vars = vars(rile * ipd),  # Main effects plus interaction
#'   fn = fn(w ~ 1/n, c = TRUE),
#'   RE = FALSE
#' )
#'
#' mm(
#'   id = id(pid, gid),
#'   vars = vars(rile + I(rile^2)),  # Quadratic term
#'   fn = fn(w ~ 1/n, c = TRUE),
#'   RE = FALSE
#' )
#' }
#'
#' @references
#' Rosche, B. (2026). A Multilevel Model for Theorizing and Estimating the
#' Micro-Macro Link. \emph{Political Analysis}.
#'
#' Browne, W. J., Goldstein, H., & Rasbash, J. (2001). Multiple membership
#' multiple classification (MMMC) models. \emph{Statistical Modelling}, 1(2), 103-124.
#'
#' Fielding, A., & Goldstein, H. (2006). \emph{Cross-classified and multiple
#' membership structures in multilevel models: An introduction and review}.
#' Research Report RR791, Department for Education and Skills.
#'
#' @export
mm <- function(id, vars = NULL, fn = NULL, RE = NULL, ar = FALSE) {


  # Validate id
  if (!inherits(id, "bml_id")) {
    stop("'id' must be specified using id(mmid, mainid)")
  }
  if (length(id) != 2) {
    stop("'id' must have exactly 2 identifiers: id(mmid, mainid)")
  }

  # Validate fn
  if (is.null(fn) || !inherits(fn, "bml_fn")) {
    stop("'fn' must be specified using fn(w ~ ..., c = TRUE/FALSE)")
  }

  # RE logic: if vars is NULL, RE must be TRUE
  if (is.null(vars)) {
    if (!is.null(RE) && isFALSE(RE)) {
      stop("RE must be TRUE when vars is NULL (otherwise mm() block has no effect)")
    }
    RE <- TRUE
  } else if (is.null(RE)) {
    RE <- FALSE
  }

  # Keep vars as-is (can be character vector or list structure)
  structure(
    list(
      id   = as.character(id),
      vars = vars,
      fn   = fn,
      RE   = RE,
      ar   = ar
    ),
    class = "bml_mm"
  )
}

#' Define a hierarchical nesting structure
#'
#' @description
#' Specifies a hierarchical (nesting) level in the model where group-level units
#' are nested within higher-level entities. Unlike multiple-membership structures,
#' each group belongs to exactly one nesting-level unit. Can model either random
#' effects or fixed effects at the nesting level.
#'
#' @param id An \code{\link{id}} object specifying the nesting-level identifier:
#'   \code{id(hmid)} where \code{hmid} identifies the higher-level units (e.g.,
#'   countries, regions).
#'
#' @param vars A \code{\link{vars}} object specifying nesting-level covariates,
#'   or \code{NULL} for intercept-only effects. Supports interactions (\code{*}, \code{:})
#'   and transformations (\code{I()}).
#'
#' @param name Unquoted variable name for nesting-level labels (optional). If
#'   provided, these labels will be displayed in model output for fixed effects.
#'
#' @param type Character; either \code{"RE"} for random effects (default) or
#'   \code{"FE"} for fixed effects at the nesting level.
#'
#' @param showFE Logical; if \code{TRUE} and \code{type = "FE"}, fixed effect
#'   estimates for each nesting-level unit are included in output. Default: \code{FALSE}.
#'
#' @param ar Logical; if \code{TRUE}, random effects evolve autoregressively across
#'   participations at the nesting level. Requires sequential participation
#'   indicators in the data. Default: \code{FALSE}.
#'
#' @return A \code{bml_hm} object containing the hierarchical specification.
#'
#' @details
#' \strong{Hierarchical vs. Multiple-Membership:}
#'
#' Hierarchical structures (\code{hm}) model strict nesting: each group belongs to
#' exactly one higher-level unit. Use \code{\link{mm}} when groups can have
#' memberships in multiple units.
#'
#' \strong{Random vs. Fixed Effects:}
#' \itemize{
#'   \item \strong{Random effects} (\code{type = "RE"}): Nesting-level units are
#'     treated as a random sample from a population. Best when you have many units
#'     and want to generalize.
#'   \item \strong{Fixed effects} (\code{type = "FE"}): Each unit gets its own
#'     parameter. Best when you have few units or want to estimate unit-specific effects.
#' }
#'
#' \strong{Cross-Classification:}
#'
#' Multiple \code{hm()} blocks create cross-classified models where groups are
#' simultaneously nested within multiple non-nested hierarchies (e.g., schools
#' within both neighborhoods and districts).
#'
#' @seealso \code{\link{bml}}, \code{\link{mm}}, \code{\link{id}}, \code{\link{vars}}
#'
#' @examples
#' \dontrun{
#' # Random effects with covariates
#' hm(
#'   id = id(cid),
#'   vars = vars(gdp + democracy),
#'   name = cname,
#'   type = "RE"
#' )
#'
#' # Random intercepts only
#' hm(
#'   id = id(cid),
#'   vars = NULL,
#'   type = "RE"
#' )
#'
#' # Fixed effects
#' hm(
#'   id = id(cid),
#'   vars = NULL,
#'   name = cname,
#'   type = "FE",
#'   showFE = TRUE  # Show estimates for each country
#' )
#'
#' # Autoregressive random effects
#' hm(
#'   id = id(cid),
#'   vars = NULL,
#'   type = "RE",
#'   ar = TRUE  # Effects evolve over time
#' )
#' }
#'
#' @references
#' Goldstein, H. (2011). \emph{Multilevel Statistical Models} (4th ed.). Wiley.
#'
#' Rabe-Hesketh, S., & Skrondal, A. (2012). \emph{Multilevel and Longitudinal
#' Modeling Using Stata} (3rd ed.). Stata Press.
#'
#' @export
hm <- function(id, vars = NULL, name = NULL, type = "RE", showFE = FALSE, ar = FALSE) {

  # Validate id
  if (!inherits(id, "bml_id")) {
    stop("'id' must be specified using id(hmid)")
  }
  if (length(id) != 1) {
    stop("'id' must have exactly 1 identifier for hm(): id(hmid)")
  }

  # Validate type
  type <- match.arg(type, c("RE", "FE"))

  # Handle name - capture as string
  name_char <- if (!missing(name) && !is.null(substitute(name))) {
    deparse(substitute(name))
  } else {
    NULL
  }

  # Keep vars as-is (can be character vector or list structure)
  structure(
    list(
      id     = as.character(id),
      vars   = vars,
      name   = name_char,
      type   = type,
      showFE = showFE,
      ar     = ar
    ),
    class = "bml_hm"
  )
}
