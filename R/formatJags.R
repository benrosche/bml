# ================================================================================================ #
# Function formatJags
# ================================================================================================ #

formatJags <- function(jags.out, monitor_spec, Ns, Xs, mm_blocks, main, hm_blocks, interactions,
                       family, cox_intervals = NULL) {

  # ========================================================================================== #
  # Flags and setup
  # ========================================================================================== #

  has_mm <- !is.null(mm_blocks) && length(mm_blocks) > 0
  has_hm <- !is.null(hm_blocks) && length(hm_blocks) > 0
  has_int <- length(interactions %||% list()) > 0
  extract_re <- monitor_has(monitor_spec, "random_effects")
  extract_weights <- monitor_has(monitor_spec, "weights") ||
    monitor_has(monitor_spec, "parameters")
  extract_pred <- monitor_has(monitor_spec, "predictive")
  extract_fitted <- monitor_has(monitor_spec, "fitted")

  n.main     <- Ns$n.main
  n.hm       <- Ns$n.hm
  n.HMN      <- Ns$n.HMN

  all_mmid_names <- Ns$all_mmid_names
  mmid_to_blocks <- Ns$mmid_to_blocks
  n.umm_list     <- Ns$n.umm_list
  n.GPN_list     <- Ns$n.GPN_list

  mainvars <- main$vars
  lhs      <- main$lhs

  # ========================================================================================== #
  # Create reg.table from JAGS output
  # ========================================================================================== #

  reg.table <-
    tibble::as_tibble(jags.out$BUGSoutput$summary[, c(1, 2, 3, 7)], rownames = "name") %>%
    dplyr::rename(mean = 2, sd = 3, lb = 4, ub = 5)

  # ========================================================================================== #
  # Add fixed (fix()) variables to reg.table
  # ========================================================================================== #

  fixed_rows <- list()

  if (!is.null(main$vars_fixed)) {
    for (i in seq_along(main$vars_fixed)) {
      var_info <- main$vars_fixed[[i]]
      fixed_rows[[length(fixed_rows) + 1]] <- data.frame(
        name = paste0("fix.main[", i, "]"),
        mean = var_info$value, sd = NA_real_, lb = NA_real_, ub = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  if (has_mm) {
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]
      if (!is.null(block$vars_fixed)) {
        for (i in seq_along(block$vars_fixed)) {
          var_info <- block$vars_fixed[[i]]
          fixed_rows[[length(fixed_rows) + 1]] <- data.frame(
            name = paste0("fix.mm.", k, "[", i, "]"),
            mean = var_info$value, sd = NA_real_, lb = NA_real_, ub = NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(fixed_rows) > 0) {
    fixed_df <- dplyr::bind_rows(fixed_rows)
    reg.table <- dplyr::bind_rows(reg.table, fixed_df)
  }

  # ========================================================================================== #
  # Extract and organize monitored quantities
  # ========================================================================================== #

  re.mm <- list()
  re.hm <- list()
  w <- list()
  pred <- c()
  fitted_mu <- c()

  # MM-level random effects (per mmid group; intercepts + slopes) ------------------------- #

    if (extract_re && has_mm && !is.null(all_mmid_names)) {
      for (g in seq_along(all_mmid_names)) {
        block_indices <- mmid_to_blocks[[all_mmid_names[g]]]
        re_idx <- block_indices[sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$RE))]
        if (length(re_idx) == 0) next
        re_spec <- mm_blocks[[re_idx[1]]]$RE
        any_ar_in_group <- any(sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$ar)))

        re.mm_raw <- reg.table %>%
          dplyr::filter(startsWith(name, paste0("re.mm.", g, "["))) %>%
          dplyr::select(-sd, -lb, -ub)

        if (any_ar_in_group && nrow(re.mm_raw) > 0) {
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
        } else if (nrow(re.mm_raw) > 0) {
          out_g <- list(Intercept = re.mm_raw %>% dplyr::pull(mean))
          for (s in seq_along(re_spec$slopes)) {
            sl_raw <- reg.table %>%
              dplyr::filter(startsWith(name, paste0("re.mm.", g, ".s", s, "["))) %>%
              dplyr::pull(mean)
            out_g[[re_spec$slopes[s]]] <- sl_raw
          }
          re.mm[[g]] <- if (length(out_g) == 1) out_g$Intercept else out_g
        }

        reg.table <- reg.table %>%
          dplyr::filter(!stringr::str_detect(name, paste0("^re\\.mm\\.", g, "(\\.s\\d+)?\\[")))
      }
    }

    # HM-level random effects ---------------------------------------------------------------- #

    if (extract_re && has_hm) {
      for (k in seq_along(hm_blocks)) {
        block <- hm_blocks[[k]]
        if (!is.null(block$RE)) {
          re.hm_raw <- reg.table %>%
            dplyr::filter(startsWith(name, paste0("re.hm.", k, "["))) %>%
            dplyr::select(-sd, -lb, -ub)

          if (!is.null(block$ar) && nrow(re.hm_raw) > 0) {
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
          } else if (nrow(re.hm_raw) > 0) {
            out_k <- list(Intercept = re.hm_raw %>% dplyr::pull(mean))
            for (s in seq_along(block$RE$slopes)) {
              sl_raw <- reg.table %>%
                dplyr::filter(startsWith(name, paste0("re.hm.", k, ".s", s, "["))) %>%
                dplyr::pull(mean)
              out_k[[block$RE$slopes[s]]] <- sl_raw
            }
            re.hm[[k]] <- if (length(out_k) == 1) out_k$Intercept else out_k
          }
        } else {
          re.hm[[k]] <- c()
        }
      }
      reg.table <- reg.table %>%
        dplyr::filter(!stringr::str_detect(name, "^re\\.hm\\.\\d+(\\.s\\d+)?\\["))
    }

    # Weights for each mm block --------------------------------------------------------------- #

    if (has_mm) {
      for (k in seq_along(mm_blocks)) {
        g <- mm_blocks[[k]]$mmid_group

        if (isTRUE(Xs$w.is.precomp[[k]])) {
          w_raw <- Xs$w.precomp[[k]]
        } else if (extract_weights) {
          w_raw <- reg.table %>%
            dplyr::filter(startsWith(name, paste0("w.", k, "["))) %>%
            dplyr::pull(mean)
        } else {
          w_raw <- numeric()
        }

        if (length(w_raw) > 0) {
          mmn_g <- Ns$mmn_list[[g]]  # member counts per main unit for this mmid group
          id2 <- cumsum(mmn_g)
          id1 <- id2 - mmn_g + 1

          wmat <- matrix(NA, nrow = n.main, ncol = max(mmn_g))
          rownames(wmat) <- paste0("Main unit ", seq_len(n.main))
          colnames(wmat) <- paste0("W", seq_len(max(mmn_g)))

          for (i in seq_len(n.main)) {
            if (mmn_g[i] > 0) {
              wmat[i, seq_len(mmn_g[i])] <- w_raw[id1[i]:id2[i]]
            }
          }

          w[[k]] <- wmat
        } else {
          w[[k]] <- NULL
        }
      }
      reg.table <- reg.table %>% dplyr::filter(!stringr::str_detect(name, "^w\\.\\d+\\["))
    }

    # Predicted values and fitted mu ----------------------------------------------------------- #

    if (extract_pred) {
      pred <- reg.table %>%
        dplyr::filter(stringr::str_detect(name, "^pred\\[")) %>%
        dplyr::pull(mean)
    }

    if (extract_fitted) {
      fitted_mu <- reg.table %>%
        dplyr::filter(stringr::str_detect(name, "^mu\\[")) %>%
        dplyr::pull(mean)
    }

    reg.table <- reg.table %>%
      dplyr::filter(!stringr::str_detect(
        name,
        "^(pred|mu|log_lik|w\\.\\d+)\\[|^re\\.(mm|hm)\\.\\d+(\\.s\\d+)?\\["
      ))

  # ========================================================================================== #
  # Rename parameters to user-facing labels
  # ========================================================================================== #

  newnames <- reg.table %>% dplyr::pull(name)
  component <- rep("other", length(newnames))

  # Main-level coefficients: b[x] (or scalar "b") -> term names
  main_indices <- stringr::str_detect(newnames, "^b(\\[|$)")
  if (any(main_indices)) {
    matched <- newnames[main_indices]
    has_bracket <- stringr::str_detect(matched, "\\[")
    main_nums <- as.numeric(ifelse(has_bracket,
                                   stringr::str_extract(matched, "(?<=\\[)\\d+(?=\\])"), "1"))
    newnames[main_indices] <- ifelse(mainvars[main_nums] == "X0", "(Intercept)",
                                     mainvars[main_nums])
    component[main_indices] <- "fixed"
  }

  # Fixed main-level variables
  if (!is.null(main$vars_fixed)) {
    fix_main_indices <- stringr::str_detect(newnames, "^fix\\.main\\[")
    if (any(fix_main_indices)) {
      fix_nums <- as.numeric(stringr::str_extract(newnames[fix_main_indices], "(?<=\\[)\\d+(?=\\])"))
      var_names <- sapply(fix_nums, function(i) main$vars_fixed[[i]]$var)
      newnames[fix_main_indices] <- ifelse(var_names == "X0", "(Intercept) (fixed)",
                                           paste0(var_names, " (fixed)"))
      component[fix_main_indices] <- "fixed"
    }
  }

  if (has_mm) {
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]
      g <- block$mmid_group

      # Feature coefficients: b.fn.k[x] -> feature labels
      if (length(block$feature_labels) > 0) {
        pattern <- paste0("^b\\.fn\\.", k, "($|\\[)")
        idxs <- stringr::str_detect(newnames, pattern)
        if (any(idxs)) {
          matched <- newnames[idxs]
          has_bracket <- stringr::str_detect(matched, "\\[")
          nums <- as.numeric(ifelse(has_bracket,
                                    stringr::str_extract(matched, "(?<=\\[)\\d+(?=\\])"), "1"))
          newnames[idxs] <- block$feature_labels[nums]
          component[idxs] <- "fixed"
        }
      }

      # fn shape parameters: fn.<p>.k -> fn[p] (mm.k)
      for (pn in names(block$fn$est_params)) {
        idxs <- which(newnames == paste0("fn.", pn, ".", k))
        if (length(idxs) > 0) {
          newnames[idxs] <- paste0("fn[", pn, "] (mm.", k, ")")
          component[idxs] <- "fn"
        }
      }

      # Weight parameters: b.w.k[p]
      if (length(block$w$params) > 0) {
        pattern <- paste0("^b\\.w\\.", k, "($|\\[)")
        idxs <- stringr::str_detect(newnames, pattern)
        if (any(idxs)) {
          matched <- newnames[idxs]
          has_bracket <- stringr::str_detect(matched, "\\[")
          nums <- as.numeric(ifelse(has_bracket,
                                    stringr::str_extract(matched, "(?<=\\[)\\d+(?=\\])"), "1"))
          if (length(block$w$vars_p) >= max(nums) && length(block$w$vars_p) > 0) {
            newnames[idxs] <- paste0("w[", block$w$vars_p[nums], "] (mm.", k, ")")
          } else {
            newnames[idxs] <- paste0("w[", block$w$params[nums], "] (mm.", k, ")")
          }
          component[idxs] <- "weights"
        }
      }

      # Fixed mm variables
      if (!is.null(block$vars_fixed)) {
        pattern <- paste0("^fix\\.mm\\.", k, "\\[")
        idxs <- stringr::str_detect(newnames, pattern)
        if (any(idxs)) {
          fix_nums <- as.numeric(stringr::str_extract(newnames[idxs], "(?<=\\[)\\d+(?=\\])"))
          var_names <- sapply(fix_nums, function(i) block$vars_fixed[[i]]$var)
          newnames[idxs] <- paste0(var_names, " (mm.", k, ", fixed)")
          component[idxs] <- "fixed"
        }
      }

      # Member fixed effects
      if (!is.null(block$FE) && block$FE$showFE) {
        pattern <- paste0("^alpha\\.mm\\.", k, "\\[")
        idxs <- stringr::str_detect(newnames, pattern)
        if (any(idxs)) {
          nums <- as.numeric(stringr::str_extract(newnames[idxs], "(?<=\\[)\\d+(?=\\])"))
          newnames[idxs] <- paste0("member ", nums, " (mm.", k, ", FE)")
          component[idxs] <- "fixed"
        }
      }
    }

    # RE sd / cor labels per mmid group
    for (g in seq_along(all_mmid_names)) {
      block_indices <- mmid_to_blocks[[all_mmid_names[g]]]
      re_idx <- block_indices[sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$RE))]
      if (length(re_idx) == 0) next
      re_spec <- mm_blocks[[re_idx[1]]]$RE

      idx <- which(newnames == paste0("sigma.mm.", g))
      if (length(idx) > 0) {
        newnames[idx] <- paste0("sd(Intercept) (mm.", re_idx[1], ")")
        component[idx] <- "random"
      }
      for (s in seq_along(re_spec$slopes)) {
        idx <- which(newnames == paste0("sigma.mm.", g, ".s", s))
        if (length(idx) > 0) {
          newnames[idx] <- paste0("sd(", re_spec$slopes[s], ") (mm.", re_idx[1], ")")
          component[idx] <- "random"
        }
      }
      idx <- which(newnames == paste0("rho.mm.", g))
      if (length(idx) > 0) {
        newnames[idx] <- paste0("cor(Intercept,", re_spec$slopes[1], ") (mm.", re_idx[1], ")")
        component[idx] <- "random"
      }
    }
  }

  # Interaction coefficients (scalar "b.int" when there is a single interaction)
  if (has_int) {
    idxs <- which(stringr::str_detect(newnames, "^b\\.int($|\\[)"))
    if (length(idxs) > 0) {
      matched <- newnames[idxs]
      has_bracket <- stringr::str_detect(matched, "\\[")
      nums <- as.numeric(ifelse(has_bracket,
                                stringr::str_extract(matched, "(?<=\\[)\\d+(?=\\])"), "1"))
      newnames[idxs] <- sapply(nums, function(t) interactions[[t]]$label)
      component[idxs] <- "fixed"
    }
  }

  # HM labels
  if (has_hm) {
    for (k in seq_along(hm_blocks)) {
      block <- hm_blocks[[k]]

      if (!is.null(block$RE)) {
        idx <- which(newnames == paste0("sigma.hm.", k))
        if (length(idx) > 0) {
          newnames[idx] <- paste0("sd(Intercept) (hm.", k, ")")
          component[idx] <- "random"
        }
        for (s in seq_along(block$RE$slopes)) {
          idx <- which(newnames == paste0("sigma.hm.", k, ".s", s))
          if (length(idx) > 0) {
            newnames[idx] <- paste0("sd(", block$RE$slopes[s], ") (hm.", k, ")")
            component[idx] <- "random"
          }
        }
        idx <- which(newnames == paste0("rho.hm.", k))
        if (length(idx) > 0) {
          newnames[idx] <- paste0("cor(Intercept,", block$RE$slopes[1], ") (hm.", k, ")")
          component[idx] <- "random"
        }
      }

      if (!is.null(block$FE) && block$FE$showFE) {
        pattern <- paste0("^alpha\\.hm\\.", k, "\\[")
        idxs <- stringr::str_detect(newnames, pattern)
        if (any(idxs)) {
          nums <- as.numeric(stringr::str_extract(newnames[idxs], "(?<=\\[)\\d+(?=\\])"))
          if ("hmlabel" %in% colnames(block$dat)) {
            unit_labels <- block$dat %>% dplyr::pull(hmlabel)
            newnames[idxs] <- paste0(unit_labels[nums], " (hm.", k, ", FE)")
          } else {
            newnames[idxs] <- paste0("unit ", nums, " (hm.", k, ", FE)")
          }
          component[idxs] <- "fixed"
        }
        for (s in seq_along(block$FE$slopes)) {
          pattern <- paste0("^alphas\\.hm\\.", k, "\\.s", s, "\\[")
          idxs <- stringr::str_detect(newnames, pattern)
          if (any(idxs)) {
            nums <- as.numeric(stringr::str_extract(newnames[idxs], "(?<=\\[)\\d+(?=\\])"))
            newnames[idxs] <- paste0("unit ", nums, ":", block$FE$slopes[s], " (hm.", k, ", FE)")
            component[idxs] <- "fixed"
          }
        }
      }
    }
  }

  # Family parameters
  idx <- which(newnames == "sigma")
  if (length(idx) > 0) { component[idx] <- "random" }
  idx <- which(newnames == "shape")
  if (length(idx) > 0) { component[idx] <- "random" }

  # ========================================================================================== #
  # Finalize reg.table
  # ========================================================================================== #

  reg.table <- reg.table %>%
    dplyr::mutate(Parameter = newnames, component = component) %>%
    dplyr::relocate(Parameter, .before = mean) %>%
    dplyr::filter(Parameter != "deviance") %>%
    tibble::column_to_rownames(var = "name")

  attr(reg.table, "estimate_type") <- "Posterior mean (MCMC)"
  attr(reg.table, "credible_interval") <- "95% equal-tailed credible intervals [2.5%, 97.5%]"
  attr(reg.table, "DIC") <- as.numeric(jags.out$BUGSoutput$DIC)

  outcome_desc <- switch(family,
    "Gaussian" = "gaussian (identity link)",
    "Binomial" = "bernoulli (logit link)",
    "Weibull"  = {
      if (length(lhs) >= 2) {
        paste0("weibull survival (duration: ", lhs[1], ", event: ", lhs[2], ")")
      } else {
        "weibull survival (log link)"
      }
    },
    "Cox" = {
      base_desc <- if (length(lhs) >= 2) {
        paste0("cox proportional hazards (duration: ", lhs[1], ", event: ", lhs[2])
      } else {
        "cox proportional hazards (log link"
      }
      if (!is.null(cox_intervals) && is.numeric(cox_intervals) && cox_intervals > 0) {
        paste0(base_desc, ", ", cox_intervals, " baseline hazard intervals)")
      } else {
        paste0(base_desc, ")")
      }
    },
    paste0(family, " (unknown link)")
  )
  attr(reg.table, "outcome_family") <- outcome_desc

  # Level specification
  level_spec_lines <- c()

  if (has_hm) {
    for (k in seq_along(hm_blocks)) {
      block <- hm_blocks[[k]]
      eff_desc <- if (!is.null(block$RE)) {
        paste0("RE = re(", if (block$RE$intercept) "1" else "0",
               if (length(block$RE$slopes) > 0) paste0(" + ", paste(block$RE$slopes, collapse = " + ")),
               if (isTRUE(block$RE$cor)) ", cor = TRUE", ")",
               if (!is.null(block$ar)) paste0(", ar", if (!is.null(block$ar$time)) paste0(" = ", block$ar$time)) else "")
      } else if (!is.null(block$FE)) {
        paste0("FE = fe(", if (block$FE$intercept) "1" else "0",
               if (length(block$FE$slopes) > 0) paste0(" + ", paste(block$FE$slopes, collapse = " + ")),
               ")")
      } else ""
      level_spec_lines <- c(level_spec_lines,
                            paste0("  hm.", k, ": ", block$id, " (", eff_desc, ")"))
    }
  }

  if (has_mm) {
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]
      parts <- paste0("fn(\"", block$fn$type, "\")")
      if (!is.null(block$RE)) {
        parts <- paste0(parts, ", RE = re(", if (block$RE$intercept) "1" else "0",
                        if (length(block$RE$slopes) > 0)
                          paste0(" + ", paste(block$RE$slopes, collapse = " + ")),
                        if (isTRUE(block$RE$cor)) ", cor = TRUE", ")")
      }
      if (!is.null(block$FE)) parts <- paste0(parts, ", FE")
      if (!is.null(block$ar)) parts <- paste0(parts, ", ar", if (!is.null(block$ar$time)) paste0(" = ", block$ar$time))
      level_spec_lines <- c(level_spec_lines,
                            paste0("  mm.", k, " [", block$name, "]: ", parts))
    }
  }

  if (length(level_spec_lines) > 0) {
    attr(reg.table, "level_spec") <- paste(level_spec_lines, collapse = "\n")
  } else {
    attr(reg.table, "level_spec") <- NULL
  }

  # ========================================================================================== #
  # Return
  # ========================================================================================== #

  list(
    reg.table = reg.table,
    w         = w,
    re.mm     = re.mm,
    re.hm     = re.hm,
    pred      = pred,
    fitted    = fitted_mu
  )
}
