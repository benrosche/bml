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


  # Validate family
  valid_families <- c("Gaussian", "Binomial", "Weibull", "Cox")
  if (!family %in% valid_families) {
    stop("Invalid family: '", family, "'. Must be one of: ", paste(valid_families, collapse = ", "))
  }

  # Validate LHS based on family
  if (family == "Gaussian" && length(lhs) > 1) {
    stop("family=\"Gaussian\" takes only one variable on the left-hand side of the formula.")
  }
  if (family == "Binomial" && length(lhs) > 1) {
    stop("family=\"Binomial\" takes only one variable on the left-hand side of the formula.")
  }
  if (family == "Weibull" && length(lhs) != 2) {
    stop("family=\"Weibull\" takes two variables on the left-hand side: 'Surv(survtime, event)'")
  }
  if (family == "Cox" && length(lhs) != 2) {
    stop("family=\"Cox\" takes two variables on the left-hand side: 'Surv(survtime, event)'")
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

  # Create a main-level formula for model.matrix() - supports interactions and I()
  # This formula only includes regular terms (no mm(), hm(), fix())
  has_intercept <- attr(terms(formula), "intercept") == 1
  if (length(regular_terms) > 0) {
    main_formula_str <- paste(regular_terms, collapse = " + ")
    if (has_intercept) {
      main_formula <- as.formula(paste("~", main_formula_str))
    } else {
      main_formula <- as.formula(paste("~ 0 +", main_formula_str))
    }
  } else {
    # No regular terms, just intercept or nothing
    if (has_intercept) {
      main_formula <- as.formula("~ 1")
    } else {
      main_formula <- NULL
    }
  }

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
  mmid_groups <- NULL

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

    # Validate: all mm() blocks must have the same mainid (but can have different mmid)
    first_mainid <- mm_list[[1]]$id[2]
    for (i in seq_along(mm_list)) {
      if (mm_list[[i]]$id[2] != first_mainid) {
        stop("All mm() blocks must share the same mainid (second element of id()). ",
             "Found: mainid='", first_mainid, "' and mainid='", mm_list[[i]]$id[2], "'")
      }
    }

    # Group blocks by mmid for tracking
    mmid_groups <- list()
    for (i in seq_along(mm_list)) {
      mmid_name <- mm_list[[i]]$id[1]
      mmid_groups[[mmid_name]] <- c(mmid_groups[[mmid_name]], i)
    }

    # Collect all unique id variable names
    all_id_vars <- unique(c(
      sapply(mm_list, function(m) m$id[1]),  # all mmid names
      first_mainid                            # mainid name
    ))

    # Validate: id variables exist in data
    if (!all(all_id_vars %in% names(data))) {
      missing_ids <- all_id_vars[!all_id_vars %in% names(data)]
      stop("ID variable(s) not found in data: ", paste(missing_ids, collapse = ", "))
    }

    # Validate: mm vars exist in data
    # Use all.vars() on formula to get base variable names (handles interactions, I())
    all_mm_vars <- unlist(lapply(mm_list, function(m) {
      if (is.null(m$vars)) return(NULL)
      if (is.list(m$vars)) {
        # New structure with formula
        base_vars <- if (!is.null(m$vars$formula)) all.vars(m$vars$formula) else m$vars$free
        fixed_vars <- if (!is.null(m$vars$fixed)) sapply(m$vars$fixed, function(x) x$var) else NULL
        c(base_vars, fixed_vars)
      } else {
        # Legacy: character vector
        as.character(m$vars)
      }
    }))
    if (length(all_mm_vars) > 0 && !all(all_mm_vars %in% names(data))) {
      missing_vars <- all_mm_vars[!all_mm_vars %in% names(data)]
      stop("mm() variable(s) not found in data: ", paste(missing_vars, collapse = ", "))
    }

    # Validate: weight function variables exist in data (excluding w, n, X0, and parameters)
    # Note: For aggregation functions like min(fdp), we need to check the base variable (fdp),
    # not the derived column name (min_fdp) which will be created in createData()
    all_fn_vars <- unlist(lapply(mm_list, function(m) m$fn$vars))
    all_fn_vars <- setdiff(all_fn_vars, c("X0", "n"))

    # Get base variables needed for aggregation functions
    all_agg_vars <- unlist(lapply(mm_list, function(m) m$fn$agg_vars))

    # Get aggregated column names (these will be created, not checked against data)
    all_agg_cols <- unlist(lapply(mm_list, function(m) {
      if (!is.null(m$fn$agg_funcs)) sapply(m$fn$agg_funcs, function(a) a$col_name)
      else NULL
    }))

    # Check regular weight function variables (excluding aggregated columns)
    regular_fn_vars <- setdiff(all_fn_vars, all_agg_cols)
    if (length(regular_fn_vars) > 0 && !all(regular_fn_vars %in% names(data))) {
      missing_vars <- regular_fn_vars[!regular_fn_vars %in% names(data)]
      stop("Weight function variable(s) not found in data: ", paste(missing_vars, collapse = ", "))
    }

    # Check base variables used in aggregation functions
    if (length(all_agg_vars) > 0 && !all(all_agg_vars %in% names(data))) {
      missing_vars <- all_agg_vars[!all_agg_vars %in% names(data)]
      stop("Variable(s) used in aggregation functions not found in data: ", paste(missing_vars, collapse = ", "))
    }

    # Check: at least one mm() block must have RE = TRUE or vars specified
    has_re <- any(sapply(mm_list, function(m) m$RE))
    has_vars <- any(sapply(mm_list, function(m) {
      if (is.null(m$vars)) return(FALSE)
      if (is.list(m$vars)) {
        !is.null(m$vars$formula) || !is.null(m$vars$free) || !is.null(m$vars$fixed)
      } else {
        length(m$vars) > 0
      }
    }))
    if (!has_re && !has_vars) {
      stop("At least one mm() block must have RE = TRUE or specify vars")
    }

    # Check: RE = TRUE can be specified for multiple blocks only if they have different mmid
    # Within blocks sharing the same mmid, only one can have RE = TRUE
    for (mmid_name in names(mmid_groups)) {
      block_indices <- mmid_groups[[mmid_name]]
      n_re_in_group <- sum(sapply(block_indices, function(i) mm_list[[i]]$RE))
      if (n_re_in_group > 1) {
        stop("RE = TRUE can only be specified for one mm() block per mmid. ",
             "Found ", n_re_in_group, " blocks with RE = TRUE sharing mmid='", mmid_name, "'")
      }
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
    # Use all.vars() on formula to get base variable names (handles interactions, I())
    all_hm_vars <- unlist(lapply(hm_list, function(h) {
      if (is.null(h$vars)) return(NULL)
      if (is.list(h$vars)) {
        # New structure with formula
        base_vars <- if (!is.null(h$vars$formula)) all.vars(h$vars$formula) else h$vars$free
        fixed_vars <- if (!is.null(h$vars$fixed)) sapply(h$vars$fixed, function(x) x$var) else NULL
        c(base_vars, fixed_vars)
      } else {
        # Legacy: character vector
        as.character(h$vars)
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
      # Filter to only variables that exist in data (excludes interaction terms like a:b)
      candidate_vars <- candidate_vars[candidate_vars %in% names(data)]

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

  # Use all.vars() to extract base variable names - handles interactions (a:b, a*b) and I()
  if (!is.null(main_formula)) {
    mainvars_base <- all.vars(main_formula)
  } else {
    mainvars_base <- c()
  }

  if (length(mainvars_base) > 0 && !all(mainvars_base %in% names(data))) {
    missing_vars <- mainvars_base[!mainvars_base %in% names(data)]
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
      main_formula   = main_formula,
      mm             = mm_list,
      mm_groups      = if (length(mm_list) > 0) mmid_groups else NULL,
      hm             = hm_list
    )
  )
}
