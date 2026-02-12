# ================================================================================================ #
# Function createData
# ================================================================================================ #

createData <- function(data, formula_parts) {

  # This function takes the parsed formula and prepares data structures for JAGS.
  # Arguments:
  # - data (df) data.frame
  # - formula_parts (list) output from dissectFormula containing lhs, mainvars, mm, hm
  # Return:
  # - Returns a list with: data, mm_blocks, main, hm_blocks

  # Ensure data is ungrouped to avoid dplyr "Adding missing grouping variables" messages
  data <- dplyr::ungroup(data)

  # Unpack formula parts -------------------------------------------------------------------- #

  lhs            <- formula_parts$lhs
  mainvars       <- formula_parts$mainvars
  mainvars_fixed <- formula_parts$mainvars_fixed
  main_formula   <- formula_parts$main_formula
  mm             <- formula_parts$mm
  hm             <- formula_parts$hm

  has_mm <- length(mm) > 0
  has_hm <- length(hm) > 0

  # ========================================================================================= #
  # Rename and regroup IDs
  # ========================================================================================= #

  # Get ID variable names
  if (has_mm) {
    mainid_name <- mm[[1]]$id[2]
    # Collect all unique mmid variable names and group blocks by mmid
    all_mmid_names <- unique(sapply(mm, function(m) m$id[1]))
    mmid_to_blocks <- list()
    for (i in seq_along(mm)) {
      mmid_nm <- mm[[i]]$id[1]
      mmid_to_blocks[[mmid_nm]] <- c(mmid_to_blocks[[mmid_nm]], i)
    }
  } else {
    mainid_name <- "mainid"
    all_mmid_names <- "mmid"
    mmid_to_blocks <- list(mmid = integer(0))
    data <- data %>% dplyr::mutate(mmid = 1, mainid = 1)
  }

  if (has_hm) {
    hmid_name <- hm[[1]]$id
  } else {
    hmid_name <- "hmid"
    data <- data %>% dplyr::mutate(hmid = 1)
  }

  # Rename mainid and hmid to standard names
  data <- data %>%
    dplyr::rename(mainid = all_of(mainid_name), hmid = all_of(hmid_name)) %>%
    dplyr::group_by(mainid) %>% dplyr::mutate(mainid = dplyr::cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::group_by(hmid) %>% dplyr::mutate(hmid = dplyr::cur_group_id()) %>% dplyr::ungroup()

  # Rename each unique mmid variable to mmid.g format and create sequential IDs
  # NAs are preserved (rows where a member doesn't belong to this mmid group)
  for (g in seq_along(all_mmid_names)) {
    mmid_orig_name <- all_mmid_names[g]
    mmid_new_name <- paste0("mmid.", g)
    data <- data %>%
      dplyr::rename(!!mmid_new_name := all_of(mmid_orig_name))
    non_na <- !is.na(data[[mmid_new_name]])
    if (any(non_na)) {
      orig_vals <- data[[mmid_new_name]]
      data[[mmid_new_name]][non_na] <- dplyr::dense_rank(orig_vals[non_na])
    }
  }

  # For backward compatibility, create 'mmid' as alias for 'mmid.1'
  if (has_mm && "mmid.1" %in% colnames(data)) {
    data$mmid <- data$mmid.1
  }

  # Sort data by mainid
  data <- data %>% dplyr::arrange(mainid)

  # ========================================================================================= #
  # MM level: Process mm() blocks
  # ========================================================================================= #

  mm_blocks <- list()

  if (has_mm) {

    # Collect all base mm variables across all mm blocks (for subsetting data)
    all_mmvars_base <- unique(unlist(lapply(mm, function(m) {
      if (is.null(m$vars)) return(NULL)
      if (is.list(m$vars)) {
        # Use all.vars() on formula to get base variable names
        base_vars <- if (!is.null(m$vars$formula)) all.vars(m$vars$formula) else m$vars$free
        fixed_vars <- if (!is.null(m$vars$fixed)) sapply(m$vars$fixed, function(x) x$var) else NULL
        c(base_vars, fixed_vars)
      } else {
        as.character(m$vars)
      }
    })))

    # Create base mm data with all base variables (include all mmid columns)
    mmid_cols <- paste0("mmid.", seq_along(all_mmid_names))
    mmdat_base <- data %>%
      dplyr::arrange(mainid) %>%
      dplyr::select(all_of(mmid_cols), mainid, any_of(all_mmvars_base))

    # Process each mm block
    for (i in seq_along(mm)) {
      block <- mm[[i]]

      # Determine which mmid group this block belongs to
      block_mmid_name <- block$id[1]
      mmid_group <- which(all_mmid_names == block_mmid_name)
      mmid_col <- paste0("mmid.", mmid_group)

      # Filter to rows belonging to this mmid group (drop NAs from other groups)
      mmdat_block <- mmdat_base %>% dplyr::filter(!is.na(.data[[mmid_col]]))

      # Extract formula and fixed vars from new structure
      vars_formula <- NULL
      vars_fixed <- NULL

      if (!is.null(block$vars) && is.list(block$vars)) {
        vars_formula <- block$vars$formula
        vars_fixed <- block$vars$fixed
      } else if (!is.null(block$vars)) {
        # Legacy: character vector - create formula
        vars_formula <- as.formula(paste("~ 0 +", paste(block$vars, collapse = " + ")))
      }

      # Build design matrix using model.matrix() - handles interactions and I()
      if (!is.null(vars_formula)) {
        X_mm <- model.matrix(vars_formula, data = mmdat_block)
        X_mm_df <- as.data.frame(X_mm)
        dat_free <- dplyr::bind_cols(
          mmdat_block %>% dplyr::select(all_of(mmid_col), mainid) %>% dplyr::rename(mmid = all_of(mmid_col)),
          X_mm_df
        )
        # Get actual column names from model.matrix (these are the "vars")
        vars_free <- colnames(X_mm_df)
      } else {
        dat_free <- NULL
        vars_free <- NULL
      }

      # Handle fixed variables
      dat_fixed <- if (!is.null(vars_fixed)) {
        var_names <- sapply(vars_fixed, function(x) x$var)
        mmdat_block %>% dplyr::select(all_of(mmid_col), mainid, all_of(var_names)) %>% dplyr::rename(mmid = all_of(mmid_col))
      } else NULL

      fix_values <- if (!is.null(vars_fixed)) {
        sapply(vars_fixed, function(x) x$value)
      } else NULL

      # Weight function data for this block
      wvars <- block$fn$vars
      agg_funcs <- block$fn$agg_funcs
      agg_vars <- block$fn$agg_vars

      wdat <- data %>%
        dplyr::filter(!is.na(.data[[mmid_col]])) %>%
        dplyr::add_count(mainid, name = "n") %>%
        dplyr::mutate(X0 = 1)

      # Compute group-level aggregates if any aggregation functions were used
      if (!is.null(agg_funcs) && length(agg_funcs) > 0) {

        # Helper function for statistical mode (most frequent value)
        stat_mode <- function(x) {
          x <- x[!is.na(x)]
          if (length(x) == 0) return(NA)
          ux <- unique(x)
          ux[which.max(tabulate(match(x, ux)))]
        }

        for (agg in agg_funcs) {
          col_name <- agg$col_name
          var_name <- agg$var
          func_name <- agg$func
          prob_val <- agg$prob  # Only used for quantile

          # Skip if column already exists (same aggregation may appear multiple times)
          if (col_name %in% names(wdat)) next

          # Compute group-level aggregate
          agg_result <- wdat %>%
            dplyr::group_by(mainid) %>%
            dplyr::summarise(
              !!col_name := switch(func_name,
                "min"      = min(.data[[var_name]], na.rm = TRUE),
                "max"      = max(.data[[var_name]], na.rm = TRUE),
                "mean"     = mean(.data[[var_name]], na.rm = TRUE),
                "sum"      = sum(.data[[var_name]], na.rm = TRUE),
                "sd"       = sd(.data[[var_name]], na.rm = TRUE),
                "var"      = var(.data[[var_name]], na.rm = TRUE),
                "first"    = dplyr::first(.data[[var_name]]),
                "last"     = dplyr::last(.data[[var_name]]),
                "median"   = median(.data[[var_name]], na.rm = TRUE),
                "mode"     = stat_mode(.data[[var_name]]),
                "range"    = max(.data[[var_name]], na.rm = TRUE) - min(.data[[var_name]], na.rm = TRUE),
                "quantile" = quantile(.data[[var_name]], probs = prob_val, na.rm = TRUE)
              ),
              .groups = "drop"
            )

          # Join back to replicate aggregate value to all members in the group
          wdat <- wdat %>%
            dplyr::left_join(agg_result, by = "mainid")
        }
      }

      wdat <- wdat %>%
        dplyr::select(all_of(mmid_col), mainid, all_of(wvars)) %>%
        dplyr::rename(mmid = all_of(mmid_col))

      mm_blocks[[i]] <- list(
        vars       = vars_free,
        vars_fixed = vars_fixed,
        dat        = dat_free,
        dat_fixed  = dat_fixed,
        fix_values = fix_values,
        wdat       = wdat,
        fn         = block$fn,
        RE         = block$RE,
        ar         = block$ar,
        mmid_group = mmid_group
      )
    }

    # Store common info
    attr(mm_blocks, "all_mmid_names") <- all_mmid_names
    attr(mm_blocks, "mmid_to_blocks") <- mmid_to_blocks
    attr(mm_blocks, "mainid_name")    <- mainid_name
    attr(mm_blocks, "has_RE")      <- any(sapply(mm, function(m) m$RE))
    attr(mm_blocks, "has_vars")    <- any(sapply(mm_blocks, function(m) {
      !is.null(m$vars) || !is.null(m$vars_fixed)
    }))

  } else {
    mm_blocks <- NULL
  }

  # ========================================================================================= #
  # HM level: Process hm() blocks
  # ========================================================================================= #

  hm_blocks <- list()

  if (has_hm) {

    for (i in seq_along(hm)) {
      block <- hm[[i]]

      # Extract formula and fixed vars from new structure
      vars_formula <- NULL
      vars_fixed <- NULL

      if (!is.null(block$vars) && is.list(block$vars)) {
        vars_formula <- block$vars$formula
        vars_fixed <- block$vars$fixed
      } else if (!is.null(block$vars)) {
        # Legacy: character vector - create formula
        vars_formula <- as.formula(paste("~ 0 +", paste(block$vars, collapse = " + ")))
      }

      hmtype <- block$type
      hmname <- block$name
      showFE <- block$showFE
      hmar   <- block$ar

      if (hmtype == "FE") {
        # Fixed effects: create dummy variables
        # Note: fix() and formula vars are not supported for FE type

        hmdat <- data %>%
          dplyr::select(hmid, any_of(hmname)) %>%
          dplyr::group_by(hmid) %>%
          dplyr::filter(row_number() == 1) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(rn = row_number(), val = 1) %>%
          tidyr::pivot_wider(names_from = hmid, names_prefix = "hmid", values_from = val, values_fill = list(val = 0)) %>%
          dplyr::rename(hmid = rn)

        if (!is.null(hmname)) {
          hmdat <- hmdat %>% dplyr::rename(hmname = all_of(hmname))
        }

        # FE variables (excluding reference category)
        hmvars_fe <- paste0("hmid", 2:nrow(hmdat))

        hm_blocks[[i]] <- list(
          id         = block$id,
          vars       = hmvars_fe,
          vars_fixed = NULL,
          dat        = hmdat,
          dat_fixed  = NULL,
          fix_values = NULL,
          name       = hmname,
          type       = hmtype,
          showFE     = showFE,
          ar         = hmar
        )

      } else {
        # Random effects - use model.matrix() for interactions and I()

        # Get base variable names for subsetting data
        all_hmvars_base <- c(
          if (!is.null(vars_formula)) all.vars(vars_formula) else NULL,
          if (!is.null(vars_fixed)) sapply(vars_fixed, function(x) x$var) else NULL
        )

        hmdat_base <- data %>%
          dplyr::select(hmid, any_of(all_hmvars_base)) %>%
          dplyr::group_by(hmid) %>%
          dplyr::filter(row_number() == 1) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(hmid)

        # Build design matrix using model.matrix() - handles interactions and I()
        if (!is.null(vars_formula)) {
          X_hm <- model.matrix(vars_formula, data = hmdat_base)
          X_hm_df <- as.data.frame(X_hm)
          dat_free <- dplyr::bind_cols(
            hmdat_base %>% dplyr::select(hmid),
            X_hm_df
          )
          # Get actual column names from model.matrix (these are the "vars")
          vars_free <- colnames(X_hm_df)
        } else {
          dat_free <- hmdat_base %>% dplyr::select(hmid)
          vars_free <- NULL
        }

        # Handle fixed variables
        dat_fixed <- if (!is.null(vars_fixed)) {
          var_names <- sapply(vars_fixed, function(x) x$var)
          hmdat_base %>% dplyr::select(hmid, all_of(var_names))
        } else NULL

        fix_values <- if (!is.null(vars_fixed)) {
          sapply(vars_fixed, function(x) x$value)
        } else NULL

        hm_blocks[[i]] <- list(
          id         = block$id,
          vars       = vars_free,
          vars_fixed = vars_fixed,
          dat        = dat_free,
          dat_fixed  = dat_fixed,
          fix_values = fix_values,
          name       = hmname,
          type       = hmtype,
          showFE     = showFE,
          ar         = hmar
        )
      }
    }

    attr(hm_blocks, "hmid_name") <- hmid_name

  } else {
    hm_blocks <- NULL
  }

  # ========================================================================================= #
  # Main level
  # ========================================================================================= #

  if (!has_mm) {
    data <- data %>% dplyr::mutate(mainid = row_number())
  }

  # Create base maindat with index variables and LHS
  # When multiple mmid groups exist (stacked data with NAs), count members from
  # the first mmid group only for the legacy mmn/mmi1/mmi2 values.
  # Per-group values are recomputed in createJagsVars.R.
  mmn_data <- if (has_mm && "mmid.1" %in% colnames(data)) {
    data %>% dplyr::filter(!is.na(.data[["mmid.1"]]))
  } else {
    data
  }
  maindat_base <- mmn_data %>%
    dplyr::arrange(mainid) %>%
    dplyr::group_by(mainid) %>%
    dplyr::add_count(mainid, name = "mmn") %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(mainid) %>%
    dplyr::mutate(mmi2 = cumsum(mmn), mmi1 = mmi2 - mmn + 1) %>%
    dplyr::select(mainid, hmid, mmi1, mmi2, mmn, all_of(lhs))

  # Create design matrix using model.matrix() - handles interactions and I()
  if (!is.null(main_formula)) {
    # Get one row per mainid from original data (which has all variables)
    main_data_for_mm <- data %>%
      dplyr::arrange(mainid) %>%
      dplyr::group_by(mainid) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(mainid)

    # model.matrix() handles interactions (a*b, a:b) and I() transformations
    X_main <- model.matrix(main_formula, data = main_data_for_mm)
    X_main_df <- as.data.frame(X_main)

    # Rename intercept column to X0 for consistency
    if ("(Intercept)" %in% colnames(X_main_df)) {
      colnames(X_main_df)[colnames(X_main_df) == "(Intercept)"] <- "X0"
    }

    # Update mainvars to reflect actual column names from model.matrix
    mainvars <- colnames(X_main_df)

    maindat <- dplyr::bind_cols(maindat_base, X_main_df)
  } else {
    # No main-level predictors
    maindat <- maindat_base
    mainvars <- c()
  }

  # Build separate data for fixed variables (from fix() syntax)
  # These are extracted from the original data, not the model.matrix output
  maindat_fixed <- if (!is.null(mainvars_fixed)) {
    fixed_var_names <- sapply(mainvars_fixed, function(x) x$var)
    # Handle X0 specially - add it to maindat if not already there
    if ("X0" %in% fixed_var_names && !"X0" %in% colnames(maindat)) {
      maindat <- maindat %>% dplyr::mutate(X0 = 1)
    }
    # For non-X0 fixed vars, ensure they're in maindat (from original data)
    other_fixed_vars <- fixed_var_names[fixed_var_names != "X0"]
    for (fv in other_fixed_vars) {
      if (!fv %in% colnames(maindat) && fv %in% colnames(data)) {
        # Get the value from original data, one per mainid
        fv_data <- data %>%
          dplyr::arrange(mainid) %>%
          dplyr::group_by(mainid) %>%
          dplyr::filter(row_number() == 1) %>%
          dplyr::ungroup() %>%
          dplyr::pull(!!rlang::sym(fv))
        maindat[[fv]] <- fv_data
      }
    }
    # Now select the fixed variable columns
    if ("X0" %in% fixed_var_names) {
      fixed_var_names_data <- c("X0", other_fixed_vars)
    } else {
      fixed_var_names_data <- other_fixed_vars
    }
    maindat %>% dplyr::select(all_of(fixed_var_names_data))
  } else NULL

  main_fix_values <- if (!is.null(mainvars_fixed)) {
    sapply(mainvars_fixed, function(x) x$value)
  } else NULL

  main <- list(
    dat        = maindat,
    vars       = mainvars,
    vars_fixed = mainvars_fixed,
    dat_fixed  = maindat_fixed,
    fix_values = main_fix_values,
    lhs        = lhs,
    formula    = main_formula
  )

  # ========================================================================================= #
  # Return
  # ========================================================================================= #

  return(
    list(
      data      = data,
      mm_blocks = mm_blocks,
      main      = main,
      hm_blocks = hm_blocks
    )
  )
}
