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
#' # Multiple-membership: parties (pid) within governments (gid)
#' id(pid, gid)
#'
#' # Hierarchical: governments within countries
#' id(cid)
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
#' # Fix a coefficient to 1.0 (standard offset)
#' fix(exposure, 1.0)
#'
#' # Use within vars() for multiple-membership models
#' vars(fix(population, 0.5) + income + education)
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
#' @param ... Unquoted variable names from your data. Variables can be separated
#'   by \code{+} or commas. To fix a coefficient, wrap the variable in
#'   \code{\link{fix}}.
#'
#' @return A \code{bml_vars} object containing:
#'   \itemize{
#'     \item \code{free}: Character vector of variables with coefficients to estimate
#'     \item \code{fixed}: List of variables with fixed coefficients (if any)
#'   }
#'   Returns \code{NULL} if no variables are specified.
#'
#' @seealso \code{\link{fix}}, \code{\link{mm}}, \code{\link{hm}}
#'
#' @examples
#' # Simple variable specification
#' vars(income + education)
#'
#' # Equivalent using comma separation
#' vars(income, education)
#'
#' # Mix free and fixed variables
#' vars(fix(exposure, 1.0) + income + education)
#'
#' # Use in mm() specification
#' \dontrun{
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

  fixed_list <- list()
  free_vars <- character()

  # Recursive function to extract variables and fix() calls
  extract_vars <- function(e) {
    if (is.call(e)) {
      if (length(e) > 0 && as.character(e[[1]]) == "fix") {
        # It's a fix() call - extract var name and value
        var_name <- as.character(e[[2]])
        value <- eval(e[[3]], parent.frame(2))
        fixed_list[[length(fixed_list) + 1]] <<- list(var = var_name, value = value)
      } else {
        # Recurse through other calls (e.g., +)
        for (i in seq_along(e)[-1]) {
          extract_vars(e[[i]])
        }
      }
    } else if (is.name(e) && as.character(e) != "") {
      # It's a free variable name
      var_name <- as.character(e)
      if (var_name != "list") {  # Filter out the list() wrapper
        free_vars <<- c(free_vars, var_name)
      }
    }
  }

  extract_vars(expr)

  if (length(free_vars) == 0 && length(fixed_list) == 0) return(NULL)

  # If only free vars (no fixed), return character vector for backward compatibility
  if (length(fixed_list) == 0) {
    return(structure(unique(free_vars), class = "bml_vars"))
  }

  # Otherwise return list structure
  structure(
    list(
      free = unique(free_vars),
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
#' }
#'
#' When \code{c = TRUE}, the weights are constrained: \eqn{\sum_{k \in group} w_k = 1}.
#'
#' @seealso \code{\link{mm}}, \code{\link{bml}}
#'
#' @examples
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
#' @references
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

  # Extract variables and parameters
  all_v <- all.vars(w)
  params <- unique(stringr::str_extract_all(fn_string, "\\bb\\d+\\b")[[1]])
  wvars <- setdiff(all_v, c("w", params))

  # Variables that have parameters (b * X pattern)
  vars_p <- stringr::str_match_all(fn_string, "b\\d+\\s*\\*\\s*([a-zA-Z_][a-zA-Z0-9_\\.]*)")[[1]]
  vars_p <- if (nrow(vars_p) > 0) vars_p[, 2] else character(0)

  structure(
    list(
      formula    = w,
      string     = fn_string,
      vars       = wvars,
      vars_p     = vars_p,
      params     = params,
      constraint = c
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
#'   aggregate, or \code{NULL} for random effects only. Variables are weighted
#'   according to the function specified in \code{fn}.
#'
#' @param fn A \code{\link{fn}} object specifying the weight function (default:
#'   \code{fn(w ~ 1/n, c = TRUE)} for equal weights). See \code{\link{fn}} for
#'   details on weight function specification.
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
#' \deqn{\text{mm}_{kj} = \sum_{i \in group_j} w_{ki} (\mathbf{x}_{ki}'\boldsymbol{\beta}_k + \alpha_{ki})}
#'
#' where:
#' \itemize{
#'   \item \eqn{w_{ki}}: Weight for member \eqn{i} in group \eqn{j} (from \code{fn})
#'   \item \eqn{\mathbf{x}_{ki}}: Member-level covariates (from \code{vars})
#'   \item \eqn{\boldsymbol{\beta}_k}: Regression coefficients (estimated)
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
#'   vars = vars(party_strength),
#'   fn = fn(w ~ b0 + b1 * tenure, c = TRUE),
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
#' @references
#' Browne, W. J., Goldstein, H., & Rasbash, J. (2001). Multiple membership
#' multiple classification (MMMC) models. \emph{Statistical Modelling}, 1(2), 103-124.
#'
#' Fielding, A., & Goldstein, H. (2006). \emph{Cross-classified and multiple
#' membership structures in multilevel models: An introduction and review}.
#' Research Report RR791, Department for Education and Skills.
#'
#' @export
mm <- function(id, vars = NULL, fn = fn(), RE = NULL, ar = FALSE) {


  # Validate id
  if (!inherits(id, "bml_id")) {
    stop("'id' must be specified using id(mmid, mainid)")
  }
  if (length(id) != 2) {
    stop("'id' must have exactly 2 identifiers: id(mmid, mainid)")
  }

  # Validate fn
  if (!inherits(fn, "bml_fn")) {
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
#'   or \code{NULL} for intercept-only effects.
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
