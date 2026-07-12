# ================================================================================================ #
# Function createData
# ================================================================================================ #

createData <- function(data, formula_parts) {

  # This function takes the parsed formula and prepares data structures for JAGS.
  # Arguments:
  # - data (df) data.frame
  # - formula_parts (list) output from dissectFormula
  # Return:
  # - Returns a list with: data, mm_blocks, main, hm_blocks

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

  if (has_mm) {
    mainid_name <- mm[[1]]$id[2]
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
    data[["mmid"]] <- 1
    data[["mainid"]] <- 1
  }

  if (has_hm) {
    hmid_name <- hm[[1]]$id
  } else {
    hmid_name <- "hmid"
    data[["hmid"]] <- 1
  }

  data <- data %>%
    dplyr::rename(mainid = all_of(mainid_name), hmid = all_of(hmid_name)) %>%
    dplyr::group_by(mainid) %>% dplyr::mutate(mainid = dplyr::cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::group_by(hmid) %>% dplyr::mutate(hmid = dplyr::cur_group_id()) %>% dplyr::ungroup()

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

  if (has_mm && "mmid.1" %in% colnames(data)) {
    data$mmid <- data$mmid.1
  }

  data <- data %>% dplyr::arrange(mainid)

  # ========================================================================================= #
  # MM level: Process mm() blocks
  # ========================================================================================= #

  mm_blocks <- list()

  if (has_mm) {

    for (i in seq_along(mm)) {
      block <- mm[[i]]

      mmid_group <- which(all_mmid_names == block$id[1])
      mmid_col <- paste0("mmid.", mmid_group)

      # Canonical member-row ordering: (mainid, mmid). All member-level structures
      # (design matrices, weights, slope data, and the mmid index vectors built in
      # createJagsVars) must share this ordering.
      mmdat_block <- data %>%
        dplyr::filter(!is.na(.data[[mmid_col]])) %>%
        dplyr::arrange(mainid, .data[[mmid_col]])

      # ----- member attribute design matrix ------------------------------------------------- #
      vars_formula <- if (!is.null(block$vars)) block$vars$formula else NULL
      vars_fixed <- if (!is.null(block$vars)) block$vars$fixed else NULL

      if (!is.null(vars_formula)) {
        X_mm <- model.matrix(vars_formula, data = mmdat_block)
        X_mm_df <- as.data.frame(X_mm)
        dat_free <- dplyr::bind_cols(
          mmdat_block %>% dplyr::select(all_of(mmid_col), mainid) %>%
            dplyr::rename(mmid = all_of(mmid_col)),
          X_mm_df
        )
        vars_free <- colnames(X_mm_df)
      } else {
        dat_free <- NULL
        vars_free <- NULL
      }

      dat_fixed <- if (!is.null(vars_fixed)) {
        var_names <- sapply(vars_fixed, function(x) x$var)
        mmdat_block %>% dplyr::select(all_of(mmid_col), mainid, all_of(var_names)) %>%
          dplyr::rename(mmid = all_of(mmid_col))
      } else NULL

      fix_values <- if (!is.null(vars_fixed)) {
        sapply(vars_fixed, function(x) x$value)
      } else NULL

      # ----- weight data --------------------------------------------------------------------- #
      wobj <- block$w
      wvars <- wobj$vars
      agg_funcs <- wobj$agg_funcs

      wdat <- data %>%
        dplyr::filter(!is.na(.data[[mmid_col]])) %>%
        dplyr::arrange(mainid, .data[[mmid_col]]) %>%
        dplyr::add_count(mainid, name = "n") %>%
        dplyr::mutate(X0 = 1)

      if (!is.null(agg_funcs) && length(agg_funcs) > 0) {

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
          prob_val <- agg$prob

          if (col_name %in% names(wdat)) next

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

          wdat <- wdat %>%
            dplyr::left_join(agg_result, by = "mainid")
        }
      }

      wdat <- wdat %>%
        dplyr::select(all_of(mmid_col), mainid, all_of(wvars)) %>%
        dplyr::rename(mmid = all_of(mmid_col))

      # ----- random/fixed effect slope data (member level) ------------------------------------ #
      eff <- block$RE %||% block$FE
      redat <- NULL
      if (!is.null(eff) && length(eff$slopes) > 0) {
        redat <- mmdat_block %>%
          dplyr::select(all_of(mmid_col), mainid, all_of(eff$slopes)) %>%
          dplyr::rename(mmid = all_of(mmid_col))
      }

      # ----- time values for a time-indexed AR walk (member rows, canonical order) ------------ #
      ardat <- NULL
      if (!is.null(block$RE) && !is.null(block$RE$ar) && !is.null(block$RE$ar$time)) {
        # id columns were renamed above; map the time variable accordingly
        tcol <- block$RE$ar$time
        if (tcol == mainid_name) tcol <- "mainid"
        if (tcol == hmid_name) tcol <- "hmid"
        ardat <- mmdat_block %>%
          dplyr::select(all_of(mmid_col), mainid, artime = all_of(tcol)) %>%
          dplyr::rename(mmid = all_of(mmid_col))
      }

      mm_blocks[[i]] <- list(
        vars           = vars_free,
        vars_fixed     = vars_fixed,
        dat            = dat_free,
        dat_fixed      = dat_fixed,
        fix_values     = fix_values,
        wdat           = wdat,
        redat          = redat,
        ardat          = ardat,
        w              = wobj,
        fn             = block$fn,
        RE             = block$RE,
        FE             = block$FE,
        ar             = block$ar,
        name           = block$name,
        feature_labels = block$feature_labels,
        attr_cols      = block$attr_cols,
        mmid_group     = mmid_group
      )
    }

    attr(mm_blocks, "all_mmid_names") <- all_mmid_names
    attr(mm_blocks, "mmid_to_blocks") <- mmid_to_blocks
    attr(mm_blocks, "mainid_name")    <- mainid_name
    attr(mm_blocks, "has_RE")         <- any(sapply(mm, function(m) !is.null(m$RE)))
    attr(mm_blocks, "has_vars")       <- any(sapply(mm_blocks, function(m) {
      !is.null(m$vars) || !is.null(m$vars_fixed)
    }))

  } else {
    mm_blocks <- NULL
  }

  # ========================================================================================= #
  # HM level: Process hm() blocks (structure only: RE/FE, labels)
  # ========================================================================================= #

  hm_blocks <- list()

  if (has_hm) {

    for (i in seq_along(hm)) {
      block <- hm[[i]]

      hmdat <- data %>%
        dplyr::select(hmid, any_of(block$labels)) %>%
        dplyr::group_by(hmid) %>%
        dplyr::filter(row_number() == 1) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(hmid)

      if (!is.null(block$labels)) {
        hmdat <- hmdat %>% dplyr::rename(hmlabel = all_of(block$labels))
      }

      hm_blocks[[i]] <- list(
        id     = block$id,
        RE     = block$RE,
        FE     = block$FE,
        labels = block$labels,
        dat    = hmdat,
        ar     = block$ar
      )
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

  # One row per mainid from original data (for main-level design + hm slopes)
  main_data_for_mm <- data %>%
    dplyr::arrange(mainid) %>%
    dplyr::group_by(mainid) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(mainid)

  if (!is.null(main_formula)) {
    X_main <- model.matrix(main_formula, data = main_data_for_mm)
    X_main_df <- as.data.frame(X_main)

    if ("(Intercept)" %in% colnames(X_main_df)) {
      colnames(X_main_df)[colnames(X_main_df) == "(Intercept)"] <- "X0"
    }

    mainvars <- colnames(X_main_df)
    maindat <- dplyr::bind_cols(maindat_base, X_main_df)
  } else {
    maindat <- maindat_base
    mainvars <- c()
  }

  # hm slope data (main-level values of the slope variables, per hm block)
  if (has_hm) {
    for (i in seq_along(hm_blocks)) {
      eff <- hm_blocks[[i]]$RE %||% hm_blocks[[i]]$FE
      if (!is.null(eff) && length(eff$slopes) > 0) {
        hm_blocks[[i]]$slopedat <- main_data_for_mm %>%
          dplyr::select(mainid, all_of(eff$slopes))
      } else {
        hm_blocks[[i]]$slopedat <- NULL
      }

      # time values for a time-indexed AR walk (one per main-level observation);
      # id columns were renamed, so map the time variable accordingly
      arspec <- hm_blocks[[i]]$ar
      if (!is.null(arspec) && !is.null(arspec$time)) {
        tcol <- arspec$time
        if (has_mm && tcol == mainid_name) tcol <- "mainid"
        if (tcol == hmid_name) tcol <- "hmid"
        hm_blocks[[i]]$artime <- main_data_for_mm %>% dplyr::pull(!!rlang::sym(tcol))
      } else {
        hm_blocks[[i]]$artime <- NULL
      }
    }
  }

  # Fixed main-level variables (fix() syntax)
  maindat_fixed <- if (!is.null(mainvars_fixed)) {
    fixed_var_names <- sapply(mainvars_fixed, function(x) x$var)
    if ("X0" %in% fixed_var_names && !"X0" %in% colnames(maindat)) {
      maindat <- maindat %>% dplyr::mutate(X0 = 1)
    }
    other_fixed_vars <- fixed_var_names[fixed_var_names != "X0"]
    for (fv in other_fixed_vars) {
      if (!fv %in% colnames(maindat) && fv %in% colnames(data)) {
        fv_data <- main_data_for_mm %>% dplyr::pull(!!rlang::sym(fv))
        maindat[[fv]] <- fv_data
      }
    }
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
      data         = data,
      mm_blocks    = mm_blocks,
      main         = main,
      hm_blocks    = hm_blocks,
      interactions = formula_parts$interactions
    )
  )
}
