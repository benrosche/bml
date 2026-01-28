# ================================================================================================ #
# Helper functions for bml formula syntax
# ================================================================================================ #

#' @title Specify ID variables
#' @description Helper function to specify level identifiers in mm() and hm() objects.
#' @param ... Unquoted variable names (e.g., id(mmid, mainid) for mm(), id(hmid) for hm())
#' @return A bml_id object containing the variable names as characters
#' @export
id <- function(...) {
  ids <- as.character(match.call(expand.dots = FALSE)$...)
  structure(ids, class = "bml_id")
}

#' @title Fix covariate coefficient
#' @description Specify a covariate with a fixed coefficient value (not estimated)
#' @param var Unquoted variable name
#' @param value Numeric value for the fixed coefficient
#' @return A bml_fix object
#' @export
fix <- function(var, value) {
  var_name <- deparse(substitute(var))
  structure(list(var = var_name, value = as.numeric(value)), class = "bml_fix")
}

#' @title Specify variables
#' @description Helper function to specify variables in mm() and hm() objects.
#' @param ... Unquoted variable names, can use + syntax (e.g., vars(x1 + x2) or vars(x1, x2))
#'            Can also include fix() calls for fixed coefficients (e.g., vars(fix(x1, 1.0) + x2))
#' @return A bml_vars object containing free and fixed variable specifications, or NULL if empty
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

#' @title Specify weight function
#' @description Helper function to specify the weight function in mm() objects.
#' @param w A formula specifying the weight function (default: w ~ 1/n)
#' @param c Logical; if TRUE (default), weights are constrained to sum to 1
#' @return A bml_fn object containing the weight function specification
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
  params <- stringr::str_extract_all(fn_string, "\\bb\\d+\\b")[[1]] %>% unique()
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

#' @title Multiple membership specification
#' @description Specifies a multiple membership structure including variables, weight function, and random effects.
#' @param id An id() object specifying member-level and group-level identifiers: id(mmid, mainid)
#' @param vars A vars() object specifying member-level variables to aggregate, or NULL for RE-only
#' @param fn A fn() object specifying the weight function (default: fn(w ~ 1/n, c = TRUE))
#' @param RE Logical; if TRUE, include random effects. Automatically TRUE if vars is NULL.
#' @param ar Logical; if TRUE, random effects evolve autoregressively across participations (default: FALSE)
#' @return A bml_mm object containing the multiple membership specification
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

#' @title Hierarchical membership specification
#' @description Specifies a hierarchical (nesting-level) structure including variables and effect type.
#' @param id An id() object specifying the nesting-level identifier: id(hmid)
#' @param vars A vars() object specifying nesting-level variables, or NULL
#' @param name Unquoted variable name for nesting-level labels (optional)
#' @param type Character; "RE" for random effects (default) or "FE" for fixed effects
#' @param showFE Logical; if TRUE and type = "FE", report the fixed effects (default: FALSE)
#' @param ar Logical; if TRUE, random effects evolve autoregressively across participations (default: FALSE)
#' @return A bml_hm object containing the hierarchical membership specification
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
