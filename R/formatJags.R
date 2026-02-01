# ================================================================================================ #
# Function formatJags
# ================================================================================================ #

formatJags <- function(jags.out, monitor, Ns, mm_blocks, main, hm_blocks, mm, hm, family, cox_intervals = NULL) {

 # ========================================================================================== #
  # Flags and setup
  # ========================================================================================== #

  has_mm <- !is.null(mm_blocks) && length(mm_blocks) > 0
  has_hm <- !is.null(hm_blocks) && length(hm_blocks) > 0
  has_mm_RE <- has_mm && attr(mm_blocks, "has_RE")

  n.umm      <- Ns$n.umm
  mmn        <- Ns$mmn
  n.main     <- Ns$n.main
  n.hm       <- Ns$n.hm
  n.GPN      <- Ns$n.GPN
  n.HMN      <- Ns$n.HMN
  n.mmblocks <- Ns$n.mmblocks

  # Per-mmid-group info
  all_mmid_names <- Ns$all_mmid_names
  mmid_to_blocks <- Ns$mmid_to_blocks
  n.umm_list     <- Ns$n.umm_list
  n.GPN_list     <- Ns$n.GPN_list

  mainvars <- main$vars
  lhs      <- main$lhs

  # Check if any mm block uses AR
  any_ar <- has_mm && any(sapply(mm_blocks, function(b) b$ar))

  # ========================================================================================== #
  # Create reg.table from JAGS output
  # ========================================================================================== #

  reg.table <-
    tibble::as_tibble(jags.out$BUGSoutput$summary[, c(1, 2, 3, 7)], rownames = "name") %>%
    dplyr::rename(mean = 2, sd = 3, lb = 4, ub = 5)

  # ========================================================================================== #
  # Add fixed variables to reg.table
  # ========================================================================================== #

  fixed_rows <- list()

  # Main-level fixed variables
  if (!is.null(main$vars_fixed)) {
    for (i in seq_along(main$vars_fixed)) {
      var_info <- main$vars_fixed[[i]]
      fixed_rows[[length(fixed_rows) + 1]] <- data.frame(
        name = paste0("fix.main[", i, "]"),
        mean = var_info$value,
        sd = NA_real_,
        lb = NA_real_,
        ub = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  # MM-level fixed variables (per block)
  if (has_mm) {
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]
      if (!is.null(block$vars_fixed)) {
        for (i in seq_along(block$vars_fixed)) {
          var_info <- block$vars_fixed[[i]]
          fixed_rows[[length(fixed_rows) + 1]] <- data.frame(
            name = paste0("fix.mm.", k, "[", i, "]"),
            mean = var_info$value,
            sd = NA_real_,
            lb = NA_real_,
            ub = NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  # HM-level fixed variables (per block)
  if (has_hm) {
    for (k in seq_along(hm_blocks)) {
      block <- hm_blocks[[k]]
      if (!is.null(block$vars_fixed)) {
        for (i in seq_along(block$vars_fixed)) {
          var_info <- block$vars_fixed[[i]]
          fixed_rows[[length(fixed_rows) + 1]] <- data.frame(
            name = paste0("fix.hm.", k, "[", i, "]"),
            mean = var_info$value,
            sd = NA_real_,
            lb = NA_real_,
            ub = NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  # Bind fixed rows to reg.table
  if (length(fixed_rows) > 0) {
    fixed_df <- dplyr::bind_rows(fixed_rows)
    reg.table <- dplyr::bind_rows(reg.table, fixed_df)
  }

  # ========================================================================================== #
  # Extract and organize outputs
  # ========================================================================================== #

  re.mm <- list()
  re.hm <- list()
  w <- list()
  pred <- c()

  if (monitor) {

    # MM-level Random Effects - per mmid group -------------------------------------------- #

    if (has_mm && !is.null(all_mmid_names)) {
      for (g in seq_along(all_mmid_names)) {
        block_indices <- mmid_to_blocks[[all_mmid_names[g]]]
        has_re_in_group <- any(sapply(block_indices, function(i) mm_blocks[[i]]$RE))
        any_ar_in_group <- any(sapply(block_indices, function(i) mm_blocks[[i]]$ar))

        if (has_re_in_group) {
          re.mm_raw <- reg.table %>%
            dplyr::filter(startsWith(name, paste0("re.mm.", g, "["))) %>%
            dplyr::select(-sd, -lb, -ub)

          if (any_ar_in_group) {
            # Autoregressive structure
            re.mm_df <- re.mm_raw %>%
              tidyr::separate(name, c("i", "j"), ",", remove = FALSE) %>%
              dplyr::mutate(
                i = as.numeric(stringr::str_remove(i, paste0("re.mm.", g, "\\["))),
                j = as.numeric(stringr::str_remove(j, "]"))
              ) %>%
              dplyr::arrange(i, j)

            n_umm_g <- n.umm_list[[g]]
            n_GPN_g <- n.GPN_list[[g]]
            remat <- matrix(NA, nrow = n_umm_g, ncol = n_GPN_g)
            rownames(remat) <- paste0("MM unit ", seq_len(n_umm_g))
            colnames(remat) <- paste0("Random walk ", seq_len(n_GPN_g))

            for (r in seq_len(nrow(re.mm_df))) {
              remat[re.mm_df$i[r], re.mm_df$j[r]] <- re.mm_df$mean[r]
            }

            re.mm[[g]] <- remat
          } else {
            re.mm[[g]] <- re.mm_raw %>% dplyr::pull(mean)
          }

          reg.table <- reg.table %>% dplyr::filter(!startsWith(name, paste0("re.mm.", g, "[")))
        }
      }
    }

    # HM-level Random Effects --------------------------------------------------------------- #

    if (has_hm) {
      for (k in seq_along(hm_blocks)) {
        block <- hm_blocks[[k]]
        if (block$type == "RE") {
          re.hm_raw <- reg.table %>%
            dplyr::filter(startsWith(name, paste0("re.hm.", k, "["))) %>%
            dplyr::select(-sd, -lb, -ub)

          if (block$ar) {
            # Autoregressive structure
            re.hm_df <- re.hm_raw %>%
              tidyr::separate(name, c("i", "j"), ",", remove = FALSE) %>%
              dplyr::mutate(
                i = as.numeric(stringr::str_remove(i, paste0("re.hm.", k, "\\["))),
                j = as.numeric(stringr::str_remove(j, "]"))
              ) %>%
              dplyr::arrange(i, j)

            remat <- matrix(NA, nrow = n.hm, ncol = n.HMN)
            rownames(remat) <- paste0("HM unit ", seq_len(n.hm))
            colnames(remat) <- paste0("Random walk ", seq_len(n.HMN))

            for (r in seq_len(nrow(re.hm_df))) {
              remat[re.hm_df$i[r], re.hm_df$j[r]] <- re.hm_df$mean[r]
            }

            re.hm[[k]] <- remat
          } else {
            # Non-AR structure
            re.hm[[k]] <- re.hm_raw %>% dplyr::pull(mean)
          }
        } else {
          re.hm[[k]] <- c()
        }
      }
      reg.table <- reg.table %>% dplyr::filter(!stringr::str_detect(name, "^re\\.hm\\.\\d+\\["))
    }

    # Weights for each mm block ------------------------------------------------------------ #

    if (has_mm) {
      for (k in seq_along(mm_blocks)) {
        w_raw <- reg.table %>%
          dplyr::filter(startsWith(name, paste0("w.", k, "["))) %>%
          dplyr::pull(mean)

        if (length(w_raw) > 0) {
          id1 <- cumsum(mmn) - mmn + 1
          id2 <- cumsum(mmn)

          wmat <- matrix(NA, nrow = n.main, ncol = max(mmn))
          rownames(wmat) <- paste0("Main unit ", seq_len(n.main))
          colnames(wmat) <- paste0("W", seq_len(max(mmn)))

          for (i in seq_len(n.main)) {
            wmat[i, seq_len(mmn[i])] <- w_raw[id1[i]:id2[i]]
          }

          w[[k]] <- wmat
        } else {
          w[[k]] <- NULL
        }
      }
      reg.table <- reg.table %>% dplyr::filter(!stringr::str_detect(name, "^w\\.\\d+\\["))
    }

    # Predicted values --------------------------------------------------------------------- #

    pred <- reg.table %>%
      dplyr::filter(startsWith(name, "pred")) %>%
      dplyr::select(-sd, -lb, -ub) %>%
      dplyr::pull(mean)

    reg.table <- reg.table %>% dplyr::filter(!startsWith(name, "pred"))
  }

  # ========================================================================================== #
  # Rename parameters to meaningful names
  # ========================================================================================== #

  newnames <- reg.table %>% dplyr::pull(name)

  # Main-level variables (estimated)
  main_indices <- stringr::str_detect(newnames, "^b\\[")
  if (any(main_indices)) {
    main_nums <- as.numeric(stringr::str_extract(newnames[main_indices], "(?<=\\[)\\d+(?=\\])"))
    newnames[main_indices] <- ifelse(mainvars[main_nums] == "X0", "Intercept", mainvars[main_nums])
  }

  # Main-level variables (fixed)
  if (!is.null(main$vars_fixed)) {
    fix_main_indices <- stringr::str_detect(newnames, "^fix\\.main\\[")
    if (any(fix_main_indices)) {
      fix_nums <- as.numeric(stringr::str_extract(newnames[fix_main_indices], "(?<=\\[)\\d+(?=\\])"))
      var_names <- sapply(fix_nums, function(i) main$vars_fixed[[i]]$var)
      newnames[fix_main_indices] <- ifelse(var_names == "X0", "Intercept (fixed)", paste0(var_names, " (fixed)"))
    }
  }

  # MM-level variables (per block)
  if (has_mm) {
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]

      # Estimated variables
      if (!is.null(block$vars) && length(block$vars) > 0) {
        # Match both b.mm.k[x] (array) and b.mm.k (scalar, when single variable)
        pattern <- paste0("^b\\.mm\\.", k, "($|\\[)")
        mm_indices <- stringr::str_detect(newnames, pattern)
        if (any(mm_indices)) {
          matched <- newnames[mm_indices]
          has_bracket <- stringr::str_detect(matched, "\\[")
          mm_nums <- as.numeric(ifelse(
            has_bracket,
            stringr::str_extract(matched, "(?<=\\[)\\d+(?=\\])"),
            "1"
          ))
          newnames[mm_indices] <- paste0(block$vars[mm_nums], " (mm.", k, ")")
        }
      }

      # Fixed variables
      if (!is.null(block$vars_fixed)) {
        pattern <- paste0("^fix\\.mm\\.", k, "\\[")
        fix_mm_indices <- stringr::str_detect(newnames, pattern)
        if (any(fix_mm_indices)) {
          fix_nums <- as.numeric(stringr::str_extract(newnames[fix_mm_indices], "(?<=\\[)\\d+(?=\\])"))
          var_names <- sapply(fix_nums, function(i) block$vars_fixed[[i]]$var)
          newnames[fix_mm_indices] <- paste0(var_names, " (mm.", k, ", fixed)")
        }
      }
    }
  }

  # HM-level variables (per block)
  if (has_hm) {
    for (k in seq_along(hm_blocks)) {
      block <- hm_blocks[[k]]

      # Estimated variables
      if (!is.null(block$vars) && length(block$vars) > 0) {
        # Match both b.hm.k[x] (array) and b.hm.k (scalar, when single variable)
        pattern <- paste0("^b\\.hm\\.", k, "($|\\[)")
        hm_indices <- stringr::str_detect(newnames, pattern)
        if (any(hm_indices)) {
          matched <- newnames[hm_indices]
          has_bracket <- stringr::str_detect(matched, "\\[")
          hm_nums <- as.numeric(ifelse(
            has_bracket,
            stringr::str_extract(matched, "(?<=\\[)\\d+(?=\\])"),
            "1"
          ))
          # Handle FE case where vars are hmid2, hmid3, etc.
          if (block$type == "FE" && !is.null(block$name)) {
            # Use the hmname labels if available
            hm_labels <- hm_blocks[[k]]$dat %>% dplyr::pull(hmname)
            newnames[hm_indices] <- hm_labels[hm_nums + 1]  # +1 because reference is excluded
          } else {
            newnames[hm_indices] <- paste0(block$vars[hm_nums], " (hm.", k, ")")
          }
        }
      }

      # Fixed variables
      if (!is.null(block$vars_fixed)) {
        pattern <- paste0("^fix\\.hm\\.", k, "\\[")
        fix_hm_indices <- stringr::str_detect(newnames, pattern)
        if (any(fix_hm_indices)) {
          fix_nums <- as.numeric(stringr::str_extract(newnames[fix_hm_indices], "(?<=\\[)\\d+(?=\\])"))
          var_names <- sapply(fix_nums, function(i) block$vars_fixed[[i]]$var)
          newnames[fix_hm_indices] <- paste0(var_names, " (hm.", k, ", fixed)")
        }
      }
    }
  }

  # Weight parameters (per block)
  if (has_mm) {
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]
      if (length(block$fn$params) > 0) {
        # Match both b.w.k[x] (array) and b.w.k (scalar, when single parameter)
        pattern <- paste0("^b\\.w\\.", k, "($|\\[)")
        w_indices <- stringr::str_detect(newnames, pattern)
        if (any(w_indices)) {
          matched <- newnames[w_indices]
          has_bracket <- stringr::str_detect(matched, "\\[")
          w_nums <- as.numeric(ifelse(
            has_bracket,
            stringr::str_extract(matched, "(?<=\\[)\\d+(?=\\])"),
            "1"
          ))
          # Use vars_p if available
          if (length(block$fn$vars_p) > 0) {
            newnames[w_indices] <- paste0(block$fn$vars_p[w_nums], " (w.", k, ")")
          } else {
            newnames[w_indices] <- paste0(block$fn$params[w_nums], " (w.", k, ")")
          }
        }
      }
    }
  }

  # Variance parameters: annotate sigma.mm.g with the mm blocks that have RE = TRUE
  if (has_mm && !is.null(all_mmid_names)) {
    for (g in seq_along(all_mmid_names)) {
      block_indices <- mmid_to_blocks[[all_mmid_names[g]]]
      re_blocks_in_group <- block_indices[sapply(block_indices, function(i) mm_blocks[[i]]$RE)]

      if (length(re_blocks_in_group) == 1) {
        mm_tag <- paste0(" (mm.", re_blocks_in_group, ")")
        sigma_mm_idx <- which(newnames == paste0("sigma.mm.", g))
        if (length(sigma_mm_idx) > 0) {
          newnames[sigma_mm_idx] <- paste0("sigma.mm.", g, mm_tag)
        }
      }
    }
  }

  # Variance parameters: annotate sigma.hm.k with the hm block index
  if (has_hm) {
    for (k in seq_along(hm_blocks)) {
      block <- hm_blocks[[k]]
      if (block$type == "RE") {
        sigma_hm_idx <- which(newnames == paste0("sigma.hm.", k))
        if (length(sigma_hm_idx) > 0) {
          newnames[sigma_hm_idx] <- paste0("sigma.hm (hm.", k, ")")
        }
      }
    }
  }

  # ========================================================================================== #
  # Standardize parameter names for consistent bracket notation
  # ========================================================================================== #

  # Standardize parameter names: drop block index when only one block exists
  single_mm <- n.mmblocks == 1
  n.hmblocks <- if (has_hm) length(hm_blocks) else 0
  single_hm <- n.hmblocks == 1
  n.main.params <- length(mainvars)

  reg.table <- reg.table %>%
    dplyr::mutate(name = dplyr::case_when(
      # b.mm: single mm block -> b.mm / b.mm[x]
      single_mm & stringr::str_detect(name, "^b\\.mm\\.\\d+\\[") ~
        stringr::str_replace(name, "^b\\.mm\\.\\d+\\[(\\d+)\\]", "b.mm[\\1]"),
      single_mm & stringr::str_detect(name, "^b\\.mm\\.\\d+$") ~
        stringr::str_replace(name, "^b\\.mm\\.\\d+$", "b.mm"),
      # b.mm: multiple mm blocks -> b.mm[k] / b.mm[k,x]
      stringr::str_detect(name, "^b\\.mm\\.\\d+\\[") ~
        stringr::str_replace(name, "^b\\.mm\\.(\\d+)\\[(\\d+)\\]", "b.mm[\\1,\\2]"),
      stringr::str_detect(name, "^b\\.mm\\.\\d+$") ~
        stringr::str_replace(name, "^b\\.mm\\.(\\d+)$", "b.mm[\\1]"),
      # fix.mm: single mm block -> fix.mm[x]
      single_mm & stringr::str_detect(name, "^fix\\.mm\\.\\d+\\[") ~
        stringr::str_replace(name, "^fix\\.mm\\.\\d+\\[(\\d+)\\]", "fix.mm[\\1]"),
      # fix.mm: multiple mm blocks -> fix.mm[k,x]
      stringr::str_detect(name, "^fix\\.mm\\.\\d+\\[") ~
        stringr::str_replace(name, "^fix\\.mm\\.(\\d+)\\[(\\d+)\\]", "fix.mm[\\1,\\2]"),
      # b.w: single mm block -> b.w / b.w[x]
      single_mm & stringr::str_detect(name, "^b\\.w\\.\\d+\\[") ~
        stringr::str_replace(name, "^b\\.w\\.\\d+\\[(\\d+)\\]", "b.w[\\1]"),
      single_mm & stringr::str_detect(name, "^b\\.w\\.\\d+$") ~
        stringr::str_replace(name, "^b\\.w\\.\\d+$", "b.w"),
      # b.w: multiple mm blocks -> b.w[k] / b.w[k,x]
      stringr::str_detect(name, "^b\\.w\\.\\d+\\[") ~
        stringr::str_replace(name, "^b\\.w\\.(\\d+)\\[(\\d+)\\]", "b.w[\\1,\\2]"),
      stringr::str_detect(name, "^b\\.w\\.\\d+$") ~
        stringr::str_replace(name, "^b\\.w\\.(\\d+)$", "b.w[\\1]"),
      # b.hm: single hm block -> b.hm / b.hm[x]
      single_hm & stringr::str_detect(name, "^b\\.hm\\.\\d+\\[") ~
        stringr::str_replace(name, "^b\\.hm\\.\\d+\\[(\\d+)\\]", "b.hm[\\1]"),
      single_hm & stringr::str_detect(name, "^b\\.hm\\.\\d+$") ~
        stringr::str_replace(name, "^b\\.hm\\.\\d+$", "b.hm"),
      # b.hm: multiple hm blocks -> b.hm[k] / b.hm[k,x]
      stringr::str_detect(name, "^b\\.hm\\.\\d+\\[") ~
        stringr::str_replace(name, "^b\\.hm\\.(\\d+)\\[(\\d+)\\]", "b.hm[\\1,\\2]"),
      stringr::str_detect(name, "^b\\.hm\\.\\d+$") ~
        stringr::str_replace(name, "^b\\.hm\\.(\\d+)$", "b.hm[\\1]"),
      # fix.hm: single hm block -> fix.hm[x]
      single_hm & stringr::str_detect(name, "^fix\\.hm\\.\\d+\\[") ~
        stringr::str_replace(name, "^fix\\.hm\\.\\d+\\[(\\d+)\\]", "fix.hm[\\1]"),
      # fix.hm: multiple hm blocks -> fix.hm[k,x]
      stringr::str_detect(name, "^fix\\.hm\\.\\d+\\[") ~
        stringr::str_replace(name, "^fix\\.hm\\.(\\d+)\\[(\\d+)\\]", "fix.hm[\\1,\\2]"),
      # sigma.hm: single hm block -> sigma.hm
      single_hm & stringr::str_detect(name, "^sigma\\.hm\\.\\d+$") ~
        stringr::str_replace(name, "^sigma\\.hm\\.\\d+$", "sigma.hm"),
      # sigma.hm: multiple hm blocks -> sigma.hm[k]
      stringr::str_detect(name, "^sigma\\.hm\\.\\d+$") ~
        stringr::str_replace(name, "^sigma\\.hm\\.(\\d+)$", "sigma.hm[\\1]"),
      # b: single parameter -> drop brackets
      n.main.params == 1 & stringr::str_detect(name, "^b\\[\\d+\\]$") ~
        stringr::str_replace(name, "^b\\[\\d+\\]$", "b"),
      TRUE ~ name
    ))

  # ========================================================================================== #
  # Finalize reg.table
  # ========================================================================================== #
  #
  # Final structure of reg.table:
  #   - Rownames: Original JAGS parameter names (e.g., "b[1]", "sigma.mm")
  #   - Columns:
  #       * Parameter: Cleaned/labeled parameter names (accessible via $Parameter)
  #       * mean: Posterior means (accessible via $mean)
  #       * sd: Posterior standard deviations (accessible via $sd)
  #       * lb: Lower 95% credible interval bound (accessible via $lb)
  #       * ub: Upper 95% credible interval bound (accessible via $ub)
  #
  # ========================================================================================== #

  reg.table <- reg.table %>%
    dplyr::mutate(Parameter = newnames) %>%
    dplyr::relocate(Parameter, .before = mean) %>%
    dplyr::filter(Parameter != "deviance") %>%
    tibble::column_to_rownames(var = "name")

  # Add metadata about posterior estimates
  attr(reg.table, "estimate_type") <- "Posterior mean (MCMC)"
  attr(reg.table, "credible_interval") <- "95% equal-tailed credible intervals [2.5%, 97.5%]"
  attr(reg.table, "DIC") <- as.numeric(jags.out$BUGSoutput$DIC)

  # Add outcome family and link information
  outcome_desc <- switch(family,
    "Gaussian" = "Gaussian (identity link)",
    "Binomial" = "Binomial (logit link)",
    "Weibull"  = {
      if (length(lhs) >= 2) {
        paste0("Weibull survival (duration: ", lhs[1], ", event: ", lhs[2], ")")
      } else {
        "Weibull survival (log link)"
      }
    },
    "Cox" = {
      base_desc <- if (length(lhs) >= 2) {
        paste0("Cox proportional hazards (duration: ", lhs[1], ", event: ", lhs[2])
      } else {
        "Cox proportional hazards (log link"
      }

      # Add interval information if using piecewise constant baseline hazard
      if (!is.null(cox_intervals) && is.numeric(cox_intervals) && cox_intervals > 0) {
        paste0(base_desc, ", ", cox_intervals, " baseline hazard intervals)")
      } else {
        paste0(base_desc, ")")
      }
    },
    paste0(family, " (unknown link)")  # fallback
  )
  attr(reg.table, "outcome_family") <- outcome_desc

  # Build level specification
  level_spec_lines <- c()

  # HM blocks
  if (has_hm) {
    for (k in seq_along(hm_blocks)) {
      block <- hm_blocks[[k]]
      id_var <- block$id

      # Effect type
      effect_type <- if (block$type == "RE") "RE" else "FE"

      # Temporal structure (only for RE)
      temporal <- if (block$type == "RE" && block$ar) "AR" else if (block$type == "RE") "indep" else "indep"

      # Format: hm.k: id (type, temporal)
      level_spec_lines <- c(level_spec_lines,
                           paste0("  hm.", k, ": ", id_var, " (", effect_type, ", ", temporal, ")"))
    }
  }

  # MM blocks
  if (has_mm) {
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]
      id_vars <- paste(mm[[k]]$id, collapse = ", ")

      if (block$RE) {
        # Show RE specification with temporal structure
        temporal <- if (block$ar) "AR" else "indep"
        level_spec_lines <- c(level_spec_lines,
                             paste0("  mm.", k, ": ", id_vars, " (RE, ", temporal, ")"))
      } else {
        # No parentheses when only vars, no RE
        level_spec_lines <- c(level_spec_lines,
                             paste0("  mm.", k, ": ", id_vars))
      }
    }
  }

  # Store as attribute
  if (length(level_spec_lines) > 0) {
    attr(reg.table, "level_spec") <- paste(level_spec_lines, collapse = "\n")
  } else {
    attr(reg.table, "level_spec") <- NULL
  }

  # ========================================================================================== #
  # Return
  # ========================================================================================== #

  if (monitor) {
    return(list(
      reg.table = reg.table,
      w         = w,
      re.mm     = re.mm,
      re.hm     = re.hm,
      pred      = pred
    ))
  } else {
    return(list(
      reg.table = reg.table,
      w         = list(),
      re.mm     = list(),
      re.hm     = list(),
      pred      = c()
    ))
  }

}
