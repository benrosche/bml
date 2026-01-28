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

  # Unpack formula parts -------------------------------------------------------------------- #

  lhs            <- formula_parts$lhs
  mainvars       <- formula_parts$mainvars
  mainvars_fixed <- formula_parts$mainvars_fixed
  mm             <- formula_parts$mm
  hm             <- formula_parts$hm

  has_mm <- length(mm) > 0
  has_hm <- length(hm) > 0

  # ========================================================================================= #
  # Rename and regroup IDs
  # ========================================================================================= #

  # Get ID variable names
  if (has_mm) {
    mmid_name   <- mm[[1]]$id[1]
    mainid_name <- mm[[1]]$id[2]
  } else {
    mmid_name   <- "mmid"
    mainid_name <- "mainid"
    data <- data %>% dplyr::mutate(mmid = 1, mainid = 1)
  }

  if (has_hm) {
    hmid_name <- hm[[1]]$id
  } else {
    hmid_name <- "hmid"
    data <- data %>% dplyr::mutate(hmid = 1)
  }

  # Rename to standard names and create sequential IDs
  data <- data %>%
    dplyr::rename(mmid = all_of(mmid_name), mainid = all_of(mainid_name), hmid = all_of(hmid_name)) %>%
    dplyr::group_by(mmid) %>% dplyr::mutate(mmid = dplyr::cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::group_by(mainid) %>% dplyr::mutate(mainid = dplyr::cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::group_by(hmid) %>% dplyr::mutate(hmid = dplyr::cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::arrange(mainid, mmid)

  # ========================================================================================= #
  # MM level: Process mm() blocks
  # ========================================================================================= #

  mm_blocks <- list()

  if (has_mm) {

    # Collect all unique mm variables across all mm blocks (both free and fixed)
    all_mmvars <- unique(unlist(lapply(mm, function(m) {
      if (is.list(m$vars) && !is.null(m$vars$free)) {
        c(m$vars$free, sapply(m$vars$fixed, function(x) x$var))
      } else {
        m$vars
      }
    })))

    # Create base mm data with all variables
    mmdat_base <- data %>%
      dplyr::arrange(mainid, mmid) %>%
      dplyr::select(mmid, mainid, all_of(all_mmvars))

    # Process each mm block
    for (i in seq_along(mm)) {
      block <- mm[[i]]

      # Separate free and fixed variables
      if (is.list(block$vars) && !is.null(block$vars$free)) {
        # New structure: list with free and fixed
        vars_free <- block$vars$free
        vars_fixed <- block$vars$fixed
      } else {
        # Old structure: character vector (backward compatibility)
        vars_free <- as.character(block$vars)
        vars_fixed <- NULL
      }

      # Build design matrices
      dat_free <- if (!is.null(vars_free) && length(vars_free) > 0) {
        mmdat_base %>% dplyr::select(mmid, mainid, all_of(vars_free))
      } else NULL

      dat_fixed <- if (!is.null(vars_fixed)) {
        var_names <- sapply(vars_fixed, function(x) x$var)
        mmdat_base %>% dplyr::select(mmid, mainid, all_of(var_names))
      } else NULL

      fix_values <- if (!is.null(vars_fixed)) {
        sapply(vars_fixed, function(x) x$value)
      } else NULL

      # Weight function data for this block
      wvars <- block$fn$vars

      wdat <- data %>%
        dplyr::add_count(mainid, name = "n") %>%
        dplyr::mutate(X0 = 1) %>%
        dplyr::select(mmid, mainid, all_of(wvars))

      mm_blocks[[i]] <- list(
        vars       = vars_free,
        vars_fixed = vars_fixed,
        dat        = dat_free,
        dat_fixed  = dat_fixed,
        fix_values = fix_values,
        wdat       = wdat,
        fn         = block$fn,
        RE         = block$RE,
        ar         = block$ar
      )
    }

    # Store common info
    attr(mm_blocks, "mmid_name")   <- mmid_name
    attr(mm_blocks, "mainid_name") <- mainid_name
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

      # Separate free and fixed variables
      if (is.list(block$vars) && !is.null(block$vars$free)) {
        # New structure: list with free and fixed
        vars_free <- block$vars$free
        vars_fixed <- block$vars$fixed
      } else {
        # Old structure: character vector (backward compatibility)
        vars_free <- as.character(block$vars)
        vars_fixed <- NULL
      }

      hmtype <- block$type
      hmname <- block$name
      showFE <- block$showFE
      hmar   <- block$ar

      if (hmtype == "FE") {
        # Fixed effects: create dummy variables
        # Note: fix() is not supported for FE type (would be redundant)

        all_hmvars <- c(vars_free, if (!is.null(vars_fixed)) sapply(vars_fixed, function(x) x$var) else NULL)

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
        # Random effects

        all_hmvars <- c(vars_free, if (!is.null(vars_fixed)) sapply(vars_fixed, function(x) x$var) else NULL)

        hmdat <- data %>%
          dplyr::select(hmid, all_of(all_hmvars)) %>%
          dplyr::group_by(hmid) %>%
          dplyr::filter(row_number() == 1) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(hmid)

        # Split into free and fixed data
        dat_free <- if (!is.null(vars_free) && length(vars_free) > 0) {
          hmdat %>% dplyr::select(hmid, all_of(vars_free))
        } else {
          hmdat %>% dplyr::select(hmid)
        }

        dat_fixed <- if (!is.null(vars_fixed)) {
          var_names <- sapply(vars_fixed, function(x) x$var)
          hmdat %>% dplyr::select(hmid, all_of(var_names))
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

  # Collect all main-level variable names (free and fixed)
  all_mainvars <- mainvars
  if (!is.null(mainvars_fixed)) {
    fixed_var_names <- sapply(mainvars_fixed, function(x) x$var)
    # Don't include X0 if it's fixed (will be added separately)
    fixed_var_names_data <- fixed_var_names[fixed_var_names != "X0"]
    all_mainvars <- c(all_mainvars, fixed_var_names_data)
  }

  maindat <- data %>%
    dplyr::arrange(mainid) %>%
    dplyr::group_by(mainid) %>%
    dplyr::add_count(mainid, name = "mmn") %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(mainid) %>%
    dplyr::mutate(mmi2 = cumsum(mmn), mmi1 = lag(mmi2) + 1) %>%
    dplyr::mutate(mmi1 = ifelse(row_number() == 1, 1, mmi1)) %>%
    dplyr::mutate(X0 = 1) %>%
    dplyr::select(mainid, hmid, mmi1, mmi2, mmn, all_of(lhs), all_of(all_mainvars))

  # Build separate data for fixed variables
  maindat_fixed <- if (!is.null(mainvars_fixed)) {
    fixed_var_names <- sapply(mainvars_fixed, function(x) x$var)
    # Handle X0 specially
    if ("X0" %in% fixed_var_names) {
      fixed_var_names_data <- c("X0", fixed_var_names[fixed_var_names != "X0"])
    } else {
      fixed_var_names_data <- fixed_var_names
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
    lhs        = lhs
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
