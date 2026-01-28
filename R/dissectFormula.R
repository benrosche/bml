# ================================================================================================ #
# Function dissectFormula
# ================================================================================================ #

dissectFormula <- function(formula, family, data) {

 # This function takes the formula object and turns it into variables for further processing.
  # Arguments:
  # - formula (formula) formula object
  # - family (str) model family
  # - data (df) data.frame to check whether specified variables are in the dataset
  # Return:
  # - Returns a list with four elements: lhs, mainvars, mm, hm

  # --------------------------------------------------------------------------------------------- #
  # Extract left- and right-hand side from formula
  # --------------------------------------------------------------------------------------------- #

  lhs <- formula[[2]] %>% all.vars()
  rhs_terms <- attr(terms(formula), "term.labels")

  # Validate LHS based on family
  if (family == "Gaussian" & length(lhs) > 1) {
    stop("family=\"Gaussian\" takes only one variable on the left-hand side of the formula.")
  }
  if (family == "Weibull" & length(lhs) != 2) {
    stop("family=\"Weibull\" takes two variables on the left-hand side: 'Surv(survtime, event)'")
  }

  # --------------------------------------------------------------------------------------------- #
  # Separate mm(), hm(), fix(), and regular terms
  # --------------------------------------------------------------------------------------------- #

  mm_terms <- rhs_terms[startsWith(rhs_terms, "mm(")]
  hm_terms <- rhs_terms[startsWith(rhs_terms, "hm(")]
  fix_terms <- rhs_terms[startsWith(rhs_terms, "fix(")]
  regular_terms <- rhs_terms[!startsWith(rhs_terms, "mm(") & !startsWith(rhs_terms, "hm(") & !startsWith(rhs_terms, "fix(")]

  # Parse fix() terms at main level
  eval_env_fix <- list2env(list(fix = fix), parent = baseenv())
  mainvars_fixed <- list()

  if (length(fix_terms) > 0) {
    for (term in fix_terms) {
      parsed <- tryCatch(
        eval(parse(text = term), envir = eval_env_fix),
        error = function(e) {
          stop("Error parsing fix() term: ", term, "\n", e$message)
        }
      )
      if (inherits(parsed, "bml_fix")) {
        mainvars_fixed[[length(mainvars_fixed) + 1]] <- parsed
      }
    }
  }

  mainvars <- c(
    if (attr(terms(formula), "intercept") == 1) "X0",
    regular_terms
  )

  # --------------------------------------------------------------------------------------------- #
  # Parse mm() blocks
  # --------------------------------------------------------------------------------------------- #

  # Create evaluation environment with helper functions
  eval_env <- list2env(list(
    id   = id,
    vars = vars,
    fn   = fn,
    mm   = mm
  ), parent = baseenv())

  mm_list <- list()

  if (length(mm_terms) > 0) {

    # Evaluate each mm() term
    mm_list <- lapply(mm_terms, function(term) {
      tryCatch(
        eval(parse(text = term), envir = eval_env),
        error = function(e) {
          stop("Error parsing mm() block: ", term, "\n", e$message)
        }
      )
    })

    # Validate: all mm() blocks must have the same id
    first_id <- mm_list[[1]]$id
    for (i in seq_along(mm_list)) {
      if (!identical(mm_list[[i]]$id, first_id)) {
        stop("All mm() blocks must have the same id(mmid, mainid). ",
             "Found: id(", paste(first_id, collapse = ", "), ") and ",
             "id(", paste(mm_list[[i]]$id, collapse = ", "), ")")
      }
    }

    # Validate: id variables exist in data
    if (!all(first_id %in% names(data))) {
      missing_ids <- first_id[!first_id %in% names(data)]
      stop("ID variable(s) not found in data: ", paste(missing_ids, collapse = ", "))
    }

    # Validate: mm vars exist in data
    all_mm_vars <- unlist(lapply(mm_list, function(m) {
      if (is.list(m$vars) && !is.null(m$vars$free)) {
        # New structure: list with free and fixed
        c(m$vars$free, sapply(m$vars$fixed, function(x) x$var))
      } else {
        # Old structure: character vector
        m$vars
      }
    }))
    if (length(all_mm_vars) > 0 && !all(all_mm_vars %in% names(data))) {
      missing_vars <- all_mm_vars[!all_mm_vars %in% names(data)]
      stop("mm() variable(s) not found in data: ", paste(missing_vars, collapse = ", "))
    }

    # Validate: weight function variables exist in data (excluding w, n, X0, and parameters)
    all_fn_vars <- unlist(lapply(mm_list, function(m) m$fn$vars))
    all_fn_vars <- setdiff(all_fn_vars, c("X0", "n"))
    if (length(all_fn_vars) > 0 && !all(all_fn_vars %in% names(data))) {
      missing_vars <- all_fn_vars[!all_fn_vars %in% names(data)]
      stop("Weight function variable(s) not found in data: ", paste(missing_vars, collapse = ", "))
    }

    # Check: at least one mm() block must have RE = TRUE or vars specified
    has_re <- any(sapply(mm_list, function(m) m$RE))
    has_vars <- any(sapply(mm_list, function(m) {
      if (is.list(m$vars)) {
        !is.null(m$vars$free) || !is.null(m$vars$fixed)
      } else {
        !is.null(m$vars)
      }
    }))
    if (!has_re && !has_vars) {
      stop("At least one mm() block must have RE = TRUE or specify vars")
    }

    # Check: at most one mm() block can have RE = TRUE
    n_re <- sum(sapply(mm_list, function(m) m$RE))
    if (n_re > 1) {
      stop("RE = TRUE cannot be specified for more than one mm() block.")
    }
  }

  # --------------------------------------------------------------------------------------------- #
  # Parse hm() blocks
  # --------------------------------------------------------------------------------------------- #

  # Create evaluation environment with helper functions
  eval_env_hm <- list2env(list(
    id   = id,
    vars = vars,
    hm   = hm
  ), parent = baseenv())

  hm_list <- list()

  if (length(hm_terms) > 0) {

    # Evaluate each hm() term
    hm_list <- lapply(hm_terms, function(term) {
      tryCatch(
        eval(parse(text = term), envir = eval_env_hm),
        error = function(e) {
          stop("Error parsing hm() block: ", term, "\n", e$message)
        }
      )
    })

    # Validate: id variables exist in data
    for (i in seq_along(hm_list)) {
      if (!hm_list[[i]]$id %in% names(data)) {
        stop("hm() id variable not found in data: ", hm_list[[i]]$id)
      }
    }

    # Validate: hm vars exist in data
    all_hm_vars <- unlist(lapply(hm_list, function(h) {
      if (is.list(h$vars) && !is.null(h$vars$free)) {
        # New structure: list with free and fixed
        c(h$vars$free, sapply(h$vars$fixed, function(x) x$var))
      } else {
        # Old structure: character vector
        h$vars
      }
    }))
    if (length(all_hm_vars) > 0 && !all(all_hm_vars %in% names(data))) {
      missing_vars <- all_hm_vars[!all_hm_vars %in% names(data)]
      stop("hm() variable(s) not found in data: ", paste(missing_vars, collapse = ", "))
    }

    # Validate: name variable exists in data
    for (i in seq_along(hm_list)) {
      if (!is.null(hm_list[[i]]$name) && !hm_list[[i]]$name %in% names(data)) {
        stop("hm() name variable not found in data: ", hm_list[[i]]$name)
      }
    }

    # Identify which mainvars are actually hm-level vars (don't vary within hm groups)
    # This moves appropriate variables from mainvars to the hm blocks
    if (length(mainvars) > 1) {  # more than just intercept
      hmid <- hm_list[[1]]$id
      candidate_vars <- mainvars[mainvars != "X0"]

      if (length(candidate_vars) > 0) {
        hmvars_detected <- data %>%
          dplyr::select(all_of(hmid), all_of(candidate_vars)) %>%
          dplyr::group_by(across(all_of(hmid))) %>%
          dplyr::summarise(across(everything(), ~var(., na.rm = TRUE)), .groups = "drop") %>%
          dplyr::select(-all_of(hmid)) %>%
          dplyr::select_if(~all(. == 0 | is.na(.))) %>%
          names()

        # Add detected hm vars to first hm block if not already specified
        if (length(hmvars_detected) > 0) {
          existing_hm_vars <- hm_list[[1]]$vars
          new_vars <- union(existing_hm_vars, hmvars_detected)
          hm_list[[1]]$vars <- new_vars

          # Remove from mainvars
          mainvars <- mainvars[!mainvars %in% hmvars_detected]
        }
      }
    }
  }

  # --------------------------------------------------------------------------------------------- #
  # Validate mainvars exist in data
  # --------------------------------------------------------------------------------------------- #

  mainvars_check <- mainvars[mainvars != "X0"]
  if (length(mainvars_check) > 0 && !all(mainvars_check %in% names(data))) {
    missing_vars <- mainvars_check[!mainvars_check %in% names(data)]
    stop("Main-level variable(s) not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Validate fixed main vars exist in data
  if (length(mainvars_fixed) > 0) {
    fixed_var_names <- sapply(mainvars_fixed, function(x) x$var)
    # Allow X0 to be fixed
    fixed_var_names_check <- fixed_var_names[fixed_var_names != "X0"]
    if (length(fixed_var_names_check) > 0 && !all(fixed_var_names_check %in% names(data))) {
      missing_vars <- fixed_var_names_check[!fixed_var_names_check %in% names(data)]
      stop("Fixed main-level variable(s) not found in data: ", paste(missing_vars, collapse = ", "))
    }
  }

  # --------------------------------------------------------------------------------------------- #
  # Return
  # --------------------------------------------------------------------------------------------- #

  return(
    list(
      lhs            = lhs,
      mainvars       = mainvars,
      mainvars_fixed = if (length(mainvars_fixed) > 0) mainvars_fixed else NULL,
      mm             = mm_list,
      hm             = hm_list
    )
  )
}
