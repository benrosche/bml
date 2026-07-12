# ================================================================================================ #
# Function createJagsVars
# ================================================================================================ #

createJagsVars <- function(data, family, mm_blocks, main, hm_blocks, interactions,
                           monitor, chains, inits, cox_intervals = NULL) {

  # Unpack main ------------------------------------------------------------------------------ #

  lhs      <- main$lhs
  mainvars <- main$vars
  maindat  <- main$dat

  # Flags -------------------------------------------------------------------------------------- #

  has_mm <- !is.null(mm_blocks) && length(mm_blocks) > 0
  has_hm <- !is.null(hm_blocks) && length(hm_blocks) > 0
  has_int <- length(interactions %||% list()) > 0

  # R-side equivalents of the JAGS math functions (for precomputation)
  r_math_env <- list(
    ilogit = stats::plogis, logit = stats::qlogis,
    probit = stats::qnorm, iprobit = stats::pnorm,
    cloglog = function(x) log(-log(1 - x)), icloglog = function(x) 1 - exp(-exp(x)),
    pow = function(a, b) a^b, log10 = log10
  )

  # ========================================================================================== #
  # Create IDs
  # ========================================================================================== #

  mainid <- maindat %>% dplyr::pull(mainid)
  hmid   <- if (has_hm) maindat %>% dplyr::pull(hmid) else c()
  n.main <- length(mainid)

  all_mmid_names <- if (has_mm) attr(mm_blocks, "all_mmid_names") else c()
  mmid_to_blocks <- if (has_mm) attr(mm_blocks, "mmid_to_blocks") else list()
  n_mmid_groups <- length(all_mmid_names)

  mmid_list   <- list()
  mmi1_list   <- list()
  mmi2_list   <- list()
  mmn_list    <- list()
  grp.mm_list <- list()
  n.mm_list   <- list()
  n.umm_list  <- list()

  n.GPN_list  <- list()
  n.GPNi_list <- list()
  n.GPn_list  <- list()
  # per group: normalized time-gap matrix for time-indexed AR walks (NULL when not used)
  gap.mm_list <- vector("list", n_mmid_groups)

  if (has_mm) {
    for (g in seq_along(all_mmid_names)) {
      block_indices <- mmid_to_blocks[[all_mmid_names[g]]]
      first_block <- mm_blocks[[block_indices[1]]]
      # wdat is already in the canonical (mainid, mmid) order from createData
      block_data <- first_block$wdat

      mmid_list[[g]] <- block_data %>% dplyr::pull(mmid)
      n.mm_list[[g]] <- length(mmid_list[[g]])
      n.umm_list[[g]] <- length(unique(mmid_list[[g]]))

      group_counts <- block_data %>%
        dplyr::group_by(mainid) %>%
        dplyr::summarise(mmn = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(mainid)

      mmn_full <- rep(0, n.main)
      mmn_full[group_counts$mainid] <- group_counts$mmn
      mmn_list[[g]] <- mmn_full

      mmi2_list[[g]] <- cumsum(mmn_list[[g]])
      mmi1_list[[g]] <- mmi2_list[[g]] - mmn_list[[g]] + 1

      grp.mm_list[[g]] <- rep(1:n.main, times = mmn_list[[g]])

      n.GPN_list[[g]]  <- block_data %>% dplyr::count(mmid) %>% dplyr::pull(n) %>% max()
      n.GPNi_list[[g]] <- block_data %>% dplyr::count(mmid) %>% dplyr::pull(n)

      # participation index: by time within member when the walk is time-indexed,
      # otherwise by row order (canonical (mainid, mmid) order)
      re_idx_g <- block_indices[sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$RE))]
      ar_spec_g <- if (length(re_idx_g) > 0) mm_blocks[[re_idx_g[1]]]$ar else NULL

      if (!is.null(ar_spec_g) && !is.null(ar_spec_g$time)) {
        ardat_g <- mm_blocks[[re_idx_g[1]]]$ardat  # canonical order, aligned with block_data
        n.GPn_list[[g]] <- stats::ave(ardat_g$artime, ardat_g$mmid,
                                      FUN = function(t) rank(t, ties.method = "first"))

        # normalized gap matrix [member, walk step]; step variance = sigma^2 * gap
        n_umm_g <- n.umm_list[[g]]
        n_GPN_g <- n.GPN_list[[g]]
        gap_mat <- matrix(1, nrow = n_umm_g, ncol = max(n_GPN_g, 2))
        times_by_member <- split(ardat_g$artime, ardat_g$mmid)
        all_gaps <- unlist(lapply(times_by_member, function(t) diff(sort(t))))
        mean_gap <- if (length(all_gaps) > 0) mean(all_gaps) else 1
        if (mean_gap <= 0) mean_gap <- 1
        for (m in names(times_by_member)) {
          t_sorted <- sort(times_by_member[[m]])
          if (length(t_sorted) >= 2) {
            gap_mat[as.integer(m), 2:length(t_sorted)] <- diff(t_sorted) / mean_gap
          }
        }
        gap.mm_list[[g]] <- gap_mat
      } else {
        n.GPn_list[[g]] <- block_data %>% dplyr::group_by(mmid) %>%
          dplyr::mutate(n = dplyr::row_number()) %>% dplyr::pull(n)
      }
    }
  }

  # Backward-compatible aliases for the first group
  mmid   <- if (has_mm && n_mmid_groups >= 1) mmid_list[[1]] else c()
  mmn    <- if (has_mm && n_mmid_groups >= 1) mmn_list[[1]] else c()
  mmi1   <- if (has_mm && n_mmid_groups >= 1) mmi1_list[[1]] else c()
  mmi2   <- if (has_mm && n_mmid_groups >= 1) mmi2_list[[1]] else c()

  # ========================================================================================== #
  # Create Ns
  # ========================================================================================== #

  n.mm   <- if (has_mm && n_mmid_groups >= 1) n.mm_list[[1]] else 0
  n.umm  <- if (has_mm && n_mmid_groups >= 1) n.umm_list[[1]] else 0
  n.hm   <- if (has_hm) length(unique(hmid)) else 0

  n.GPN  <- if (has_mm && n_mmid_groups >= 1) n.GPN_list[[1]] else c()
  n.GPNi <- if (has_mm && n_mmid_groups >= 1) n.GPNi_list[[1]] else c()
  n.GPn  <- if (has_mm && n_mmid_groups >= 1) n.GPn_list[[1]] else c()

  n.HMN  <- if (has_hm) maindat %>% dplyr::count(hmid) %>% dplyr::pull(n) %>% max() else c()
  n.HMNi <- if (has_hm) maindat %>% dplyr::count(hmid) %>% dplyr::pull(n) else c()
  n.HMn  <- if (has_hm) maindat %>% dplyr::group_by(hmid) %>%
    dplyr::mutate(n = row_number()) %>% dplyr::pull(n) else c()

  # Time-indexed hm AR walk: order each unit's observations by time and build gaps.
  # The walk position vector n.HMn is shared across hm blocks, so at most one hm block
  # can define a time-indexed ordering.
  gap.hm_list <- vector("list", if (has_hm) length(hm_blocks) else 0)
  if (has_hm) {
    ar_time_blocks <- which(sapply(hm_blocks, function(b) {
      !is.null(b$ar) && !is.null(b$ar$time)
    }))
    if (length(ar_time_blocks) > 1) {
      stop("Only one hm() block can carry a time-indexed AR walk (the walk ordering ",
           "is shared across hm blocks).", call. = FALSE)
    }
    if (length(ar_time_blocks) == 1) {
      k_ar <- ar_time_blocks[1]
      artime <- hm_blocks[[k_ar]]$artime
      dup <- duplicated(data.frame(hmid = hmid, t = artime))
      if (any(dup)) {
        stop("re(ar = ", hm_blocks[[k_ar]]$ar$time, "): duplicated time values within a ",
             "nesting unit (a zero time gap makes the random-walk step variance zero).",
             call. = FALSE)
      }
      n.HMn <- stats::ave(artime, hmid, FUN = function(t) rank(t, ties.method = "first"))

      gap_mat <- matrix(1, nrow = n.hm, ncol = max(n.HMN, 2))
      times_by_unit <- split(artime, hmid)
      all_gaps <- unlist(lapply(times_by_unit, function(t) diff(sort(t))))
      mean_gap <- if (length(all_gaps) > 0) mean(all_gaps) else 1
      if (mean_gap <= 0) mean_gap <- 1
      for (u in names(times_by_unit)) {
        t_sorted <- sort(times_by_unit[[u]])
        if (length(t_sorted) >= 2) {
          gap_mat[as.integer(u), 2:length(t_sorted)] <- diff(t_sorted) / mean_gap
        }
      }
      gap.hm_list[[k_ar]] <- gap_mat
    }
  }

  n.mmblocks <- if (has_mm) length(mm_blocks) else 0

  # ========================================================================================== #
  # Per-block computation modes and matrices
  # ========================================================================================== #

  n_blocks <- if (has_mm) length(mm_blocks) else 0

  X.mm       <- vector("list", n_blocks)
  offset.mm  <- vector("list", n_blocks)
  X.w        <- vector("list", n_blocks)
  X.re       <- vector("list", n_mmid_groups)  # per mmid group: RE slope covariates
  X.fe       <- vector("list", n_blocks)       # per block: FE slope covariates
  w.precomp  <- vector("list", n_blocks)
  w.is.precomp <- as.list(rep(FALSE, n_blocks))
  F.precomp  <- vector("list", n_blocks)
  F.is.precomp <- as.list(rep(FALSE, n_blocks))
  n.Xmm      <- as.list(rep(0, n_blocks))
  n.Xw       <- as.list(rep(0, n_blocks))

  w_in_jags <- function(block) length(block$w$params) > 0
  f_in_jags <- function(block) {
    (length(block$feature_labels) > 0) &&
      (w_in_jags(block) || length(block$fn$est_params) > 0)
  }

  if (has_mm) {
    for (i in seq_along(mm_blocks)) {
      block <- mm_blocks[[i]]
      g <- block$mmid_group

      # Member attribute design matrix
      if (!is.null(block$vars) && length(block$vars) > 0) {
        X.mm[[i]] <- block$dat %>% dplyr::select(all_of(block$vars)) %>% as.matrix()
        n.Xmm[[i]] <- ncol(X.mm[[i]])
      }

      # Fixed-coefficient offset
      if (!is.null(block$dat_fixed)) {
        var_names <- sapply(block$vars_fixed, function(x) x$var)
        X_fix <- block$dat_fixed %>% dplyr::select(all_of(var_names)) %>% as.matrix()
        offset.mm[[i]] <- as.vector(X_fix %*% block$fix_values)
      }

      # Weight covariate matrix
      wvars <- block$w$vars
      if (length(wvars) > 0) {
        X.w[[i]] <- block$wdat %>% dplyr::select(all_of(wvars)) %>% as.matrix()
        n.Xw[[i]] <- ncol(X.w[[i]])
      }

      # ---------------- Weight precomputation (no free parameters) ---------------- #
      if (!w_in_jags(block)) {
        wdat_env <- c(as.list(block$wdat), r_math_env)
        uw <- eval(parse(text = block$w$string), envir = wdat_env, enclos = baseenv())
        if (length(uw) == 1) uw <- rep(uw, nrow(block$wdat))
        uw <- as.numeric(uw)
        if (block$w$constraint) {
          sum_uw <- stats::ave(uw, block$wdat$mainid, FUN = sum)
          w.precomp[[i]] <- uw / sum_uw
        } else {
          w.precomp[[i]] <- uw
        }
        w.is.precomp[[i]] <- TRUE
      }

      # ---------------- Feature precomputation (fast path) ---------------- #
      if (length(block$feature_labels) > 0 && !f_in_jags(block)) {

        wv <- w.precomp[[i]]
        grp <- block$wdat$mainid
        n_feat <- length(block$feature_labels)
        fnobj <- block$fn
        type <- fnobj$type

        group_sum <- function(v) {
          out <- rep(0, n.main)
          s <- tapply(v, grp, sum)
          out[as.numeric(names(s))] <- as.numeric(s)
          out
        }

        if (type == "sum") {
          Fmat <- sapply(seq_len(ncol(X.mm[[i]])), function(x) group_sum(wv * X.mm[[i]][, x]))
          F.precomp[[i]] <- if (n_feat == 1) as.vector(Fmat) else Fmat
        } else if (type == "var") {
          m <- fnobj$moment %||% 2
          A <- group_sum(wv * X.mm[[i]][, 1])
          F.precomp[[i]] <- group_sum(wv * (X.mm[[i]][, 1] - A[grp])^m)
        } else if (type == "hhi") {
          F.precomp[[i]] <- group_sum(wv^2)
        } else if (type == "effn") {
          F.precomp[[i]] <- 1 / group_sum(wv^2)
        } else if (type == "entropy") {
          F.precomp[[i]] <- -group_sum(wv * log(wv))
        } else if (type == "threshold") {
          cv <- fnobj$fixed_params$c
          kv <- fnobj$fixed_params$kappa
          F.precomp[[i]] <- group_sum(wv * stats::plogis(kv * (X.mm[[i]][, 1] - cv)))
        } else if (type == "smax") {
          kv <- fnobj$fixed_params$kappa
          F.precomp[[i]] <- (1 / kv) * log(group_sum(wv * exp(kv * X.mm[[i]][, 1])))
        } else if (type == "gmean") {
          pv <- fnobj$fixed_params$p
          F.precomp[[i]] <- group_sum(wv * X.mm[[i]][, 1]^pv)^(1 / pv)
        } else if (type == "cov") {
          A1 <- group_sum(wv * X.mm[[i]][, 1])
          A2 <- group_sum(wv * X.mm[[i]][, 2])
          F.precomp[[i]] <- group_sum(wv * (X.mm[[i]][, 1] - A1[grp]) * (X.mm[[i]][, 2] - A2[grp]))
        } else if (type == "expr") {
          fixedp <- lapply(fnobj$fixed_params %||% list(), identity)
          Fv <- rep(NA_real_, n.main)
          for (j in unique(grp)) {
            rows <- which(grp == j)
            member_env <- list(w = wv[rows])
            for (a in seq_along(block$attr_cols)) {
              col <- match(block$attr_cols[a], block$vars)
              member_env[[block$attr_cols[a]]] <- X.mm[[i]][rows, col]
            }
            Fv[j] <- dsl_eval_group(fnobj$graph, member_env, fixedp)
          }
          Fv[is.na(Fv)] <- 0
          F.precomp[[i]] <- Fv
        }
        F.is.precomp[[i]] <- TRUE
      }

      # ---------------- FE slope covariates ---------------- #
      if (!is.null(block$FE) && length(block$FE$slopes) > 0) {
        X.fe[[i]] <- block$redat %>% dplyr::select(all_of(block$FE$slopes)) %>% as.matrix()
      }
    }

    # RE slope covariates: one matrix per mmid group (from the RE-carrying block)
    for (g in seq_along(all_mmid_names)) {
      block_indices <- mmid_to_blocks[[all_mmid_names[g]]]
      re_idx <- block_indices[sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$RE))]
      if (length(re_idx) > 0) {
        rb <- mm_blocks[[re_idx[1]]]
        if (length(rb$RE$slopes) > 0) {
          X.re[[g]] <- rb$redat %>% dplyr::select(all_of(rb$RE$slopes)) %>% as.matrix()
        }
      }
    }
  }

  # ========================================================================================== #
  # Main and HM matrices
  # ========================================================================================== #

  if (length(mainvars) > 0) {
    X.main <- maindat %>% dplyr::select(all_of(mainvars)) %>% as.matrix()
    n.Xmain <- ncol(X.main)
  } else {
    X.main <- NULL
    n.Xmain <- 0
  }

  if (!is.null(main$dat_fixed)) {
    var_names <- sapply(main$vars_fixed, function(x) x$var)
    var_names_data <- if ("X0" %in% var_names) c("X0", var_names[var_names != "X0"]) else var_names
    X.main.fix <- main$dat_fixed %>% dplyr::select(all_of(var_names_data)) %>% as.matrix()
    offset.main <- as.vector(X.main.fix %*% main$fix_values)
    n.Xmain.fix <- length(main$fix_values)
  } else {
    offset.main <- NULL
    n.Xmain.fix <- 0
  }

  # HM slope covariates (main-level values)
  X.re.hm <- vector("list", if (has_hm) length(hm_blocks) else 0)
  if (has_hm) {
    for (i in seq_along(hm_blocks)) {
      block <- hm_blocks[[i]]
      eff <- block$RE %||% block$FE
      if (!is.null(eff) && length(eff$slopes) > 0) {
        X.re.hm[[i]] <- block$slopedat %>%
          dplyr::arrange(mainid) %>%
          dplyr::select(all_of(eff$slopes)) %>% as.matrix()
      }
    }
  }

  # ========================================================================================== #
  # Build jags.params (monitors)
  # ========================================================================================== #

  jags.params <- c()

  if (has_mm) {
    # RE variance/effect monitors per mmid group
    for (g in seq_along(all_mmid_names)) {
      block_indices <- mmid_to_blocks[[all_mmid_names[g]]]
      re_idx <- block_indices[sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$RE))]
      if (length(re_idx) > 0) {
        re_spec <- mm_blocks[[re_idx[1]]]$RE
        if (re_spec$intercept || any(sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$ar)))) {
          jags.params <- c(jags.params, paste0("sigma.mm.", g))
        }
        for (s in seq_along(re_spec$slopes)) {
          jags.params <- c(jags.params, paste0("sigma.mm.", g, ".s", s))
        }
        if (isTRUE(re_spec$cor) && length(re_spec$slopes) == 1) {
          jags.params <- c(jags.params, paste0("rho.mm.", g))
        }
        if (monitor) {
          jags.params <- c(jags.params, paste0("re.mm.", g))
          for (s in seq_along(re_spec$slopes)) {
            jags.params <- c(jags.params, paste0("re.mm.", g, ".s", s))
          }
        }
      }
    }

    # Feature coefficients, shape params, weight params, weights, FE
    for (i in seq_along(mm_blocks)) {
      block <- mm_blocks[[i]]
      if (length(block$feature_labels) > 0) {
        jags.params <- c(jags.params, paste0("b.fn.", i))
      }
      for (pn in names(block$fn$est_params)) {
        jags.params <- c(jags.params, paste0("fn.", pn, ".", i))
      }
      if (length(block$w$params) > 0) {
        jags.params <- c(jags.params, paste0("b.w.", i))
      }
      # weights can only be monitored when they are model nodes (not pre-computed data)
      if (monitor && w_in_jags(block)) {
        jags.params <- c(jags.params, paste0("w.", i))
      }
      if (!is.null(block$FE) && block$FE$showFE) {
        if (block$FE$intercept) jags.params <- c(jags.params, paste0("alpha.mm.", i))
        for (s in seq_along(block$FE$slopes)) {
          jags.params <- c(jags.params, paste0("alphas.mm.", i, ".s", s))
        }
      }
    }
  }

  # Interaction coefficients
  if (has_int) {
    jags.params <- c(jags.params, "b.int")
  }

  # Main-level parameters
  if (family %in% c("Gaussian", "Weibull")) {
    jags.params <- c(jags.params, "sigma")
  }
  if (n.Xmain > 0) {
    jags.params <- c(jags.params, "b")
  }
  if (monitor) {
    jags.params <- c(jags.params, "pred", "mu")
    if (family %in% c("Gaussian", "Binomial")) {
      jags.params <- c(jags.params, "log_lik")
    }
  }

  # HM-level parameters
  if (has_hm) {
    for (i in seq_along(hm_blocks)) {
      block <- hm_blocks[[i]]
      if (!is.null(block$RE)) {
        if (block$RE$intercept || !is.null(block$ar)) {
          jags.params <- c(jags.params, paste0("sigma.hm.", i))
        }
        for (s in seq_along(block$RE$slopes)) {
          jags.params <- c(jags.params, paste0("sigma.hm.", i, ".s", s))
        }
        if (isTRUE(block$RE$cor) && length(block$RE$slopes) == 1) {
          jags.params <- c(jags.params, paste0("rho.hm.", i))
        }
        if (monitor) {
          jags.params <- c(jags.params, paste0("re.hm.", i))
          for (s in seq_along(block$RE$slopes)) {
            jags.params <- c(jags.params, paste0("re.hm.", i, ".s", s))
          }
        }
      } else if (!is.null(block$FE) && block$FE$showFE) {
        if (block$FE$intercept) jags.params <- c(jags.params, paste0("alpha.hm.", i))
        for (s in seq_along(block$FE$slopes)) {
          jags.params <- c(jags.params, paste0("alphas.hm.", i, ".s", s))
        }
      }
    }
  }

  # ========================================================================================== #
  # Family-specific data and inits
  # ========================================================================================== #

  if (family %in% c("Gaussian", "Binomial")) {

    Y <- maindat %>% dplyr::rename(Y = all_of(lhs)) %>% dplyr::pull(Y)
    jags.inits <- list()
    Ys <- list(Y = Y)

  } else if (family == "Weibull") {

    t <- maindat %>%
      dplyr::rename(t = all_of(lhs[1]), ev = all_of(lhs[2])) %>%
      dplyr::mutate(t = dplyr::case_when(ev == 0 ~ NA_real_, TRUE ~ t)) %>%
      dplyr::pull(t)

    ct.lb <- maindat %>%
      dplyr::rename(t = all_of(lhs[1]), ev = all_of(lhs[2])) %>%
      dplyr::mutate(ct.lb = t + ev) %>%
      dplyr::pull(ct.lb)

    event <- maindat %>% dplyr::rename(ev = all_of(lhs[2])) %>% dplyr::pull(ev)
    censored <- 1 - event

    jags.params <- c(jags.params, "shape")

    t.init <- t
    t.init[] <- NA
    t.init[censored == 1] <- ct.lb[censored == 1] + 1

    jags.inits <- list(t = t.init, shape = 1)
    Ys <- list(t = t, ct.lb = ct.lb, event = event, censored = censored)

  } else if (family == "Cox") {

    t <- maindat %>% dplyr::rename(t = all_of(lhs[1]), ev = all_of(lhs[2])) %>% dplyr::pull(t)
    event <- maindat %>% dplyr::rename(ev = all_of(lhs[2])) %>% dplyr::pull(ev)

    if (!is.null(cox_intervals) && is.numeric(cox_intervals) && cox_intervals > 0) {

      n.intervals <- as.integer(cox_intervals)

      event_times <- t[event == 1]
      if (length(event_times) == 0) {
        stop("No events observed in the data. Cox model cannot be estimated.")
      }

      time_breaks <- quantile(event_times, probs = seq(0, 1, length.out = n.intervals + 1))
      time_breaks[1] <- 0
      time_breaks[length(time_breaks)] <- max(t) + 1

      Y_interval <- matrix(data = 0, nrow = n.main, ncol = n.intervals)
      dN_interval <- matrix(data = 0, nrow = n.main, ncol = n.intervals)

      for (j in 1:n.main) {
        t_j <- t[j]
        event_j <- event[j]

        for (k in 1:n.intervals) {
          interval_start <- time_breaks[k]
          interval_end <- time_breaks[k + 1]

          if (t_j <= interval_start) {
            Y_interval[j, k] <- 0
          } else if (t_j > interval_end) {
            Y_interval[j, k] <- interval_end - interval_start
          } else {
            Y_interval[j, k] <- t_j - interval_start
          }

          if (event_j == 1) {
            if (t_j > interval_start && t_j <= interval_end) {
              dN_interval[j, k] <- 1
            }
          }
        }
      }

      jags.inits <- list(lambda0 = rep(0.01, n.intervals))
      Ys <- list(
        Y_interval = Y_interval,
        dN_interval = dN_interval,
        t = t,
        time_breaks = time_breaks,
        event = event,
        c = 0.001,
        d = 0.1,
        n.intervals = n.intervals
      )

    } else {

      t.unique <- c(sort(unique(t)), max(t) + 1)
      n.tu <- length(t.unique) - 1

      Y <- matrix(data = NA, nrow = n.main, ncol = n.tu)
      dN <- matrix(data = NA, nrow = n.main, ncol = n.tu)
      for (j in 1:n.main) {
        for (k in 1:n.tu) {
          Y[j, k] <- as.numeric(t[j] - t.unique[k] + 1e-05 >= 0)
          dN[j, k] <- Y[j, k] * event[j] * as.numeric(t.unique[k + 1] - t[j] >= 1e-05)
        }
      }

      jags.inits <- list(dL0 = rep(1.0, n.tu))
      Ys <- list(Y = Y, dN = dN, t = t, t.unique = t.unique, event = event,
                 c = 0.001, d = 0.1, n.tu = n.tu)
    }
  }

  # ========================================================================================== #
  # Finalize inits
  # ========================================================================================== #

  if (has_mm) {
    for (i in seq_along(mm_blocks)) {
      block <- mm_blocks[[i]]
      if (length(block$w$params) > 0) {
        jags.inits[[paste0("b.w.", i)]] <- rep(0, length(block$w$params))
      }
      for (pn in names(block$fn$est_params)) {
        jags.inits[[paste0("fn.", pn, ".", i)]] <- block$fn$est_params[[pn]]$init
      }
    }
  }

  jags.inits <- c(jags.inits, inits)
  jags.inits <- lapply(1:chains, function(x) jags.inits)

  # ========================================================================================== #
  # Build the data list for JAGS
  # ========================================================================================== #

  jags.data.list <- list(
    n.main = n.main
  )

  if (family %in% c("Gaussian", "Binomial")) {
    jags.data.list$Y <- Ys$Y
  } else if (family == "Weibull") {
    jags.data.list$t <- Ys$t
    jags.data.list$ct.lb <- Ys$ct.lb
    jags.data.list$censored <- Ys$censored
  } else if (family == "Cox") {
    if (!is.null(cox_intervals) && is.numeric(cox_intervals) && cox_intervals > 0) {
      jags.data.list$Y_interval <- Ys$Y_interval
      jags.data.list$dN_interval <- Ys$dN_interval
      jags.data.list$n.intervals <- Ys$n.intervals
      jags.data.list$c <- Ys$c
      jags.data.list$d <- Ys$d
    } else {
      jags.data.list$Y <- Ys$Y
      jags.data.list$dN <- Ys$dN
      jags.data.list$t.unique <- Ys$t.unique
      jags.data.list$n.tu <- Ys$n.tu
      jags.data.list$c <- Ys$c
      jags.data.list$d <- Ys$d
    }
  }

  if (n.Xmain > 0) {
    jags.data.list$X.main <- X.main
  }
  if (n.Xmain.fix > 0) {
    jags.data.list$offset.main <- offset.main
  }

  if (has_mm) {

    # Member-loop types that need group-scalar lookups at the member level
    needs_grp_types <- c("var", "cov", "expr")
    member_loop_types <- c("var", "entropy", "threshold", "smax", "gmean", "cov", "expr")

    for (g in seq_along(all_mmid_names)) {
      block_indices <- mmid_to_blocks[[all_mmid_names[g]]]

      # index ranges are only referenced when something in this group is computed in JAGS
      needs_ranges <- any(sapply(block_indices, function(i) {
        b <- mm_blocks[[i]]
        f_in_jags(b) || !is.null(b$RE) || !is.null(b$FE) ||
          !is.null(offset.mm[[i]]) || (w_in_jags(b) && b$w$constraint)
      }))
      if (needs_ranges) {
        jags.data.list[[paste0("mmi1.", g)]] <- mmi1_list[[g]]
        jags.data.list[[paste0("mmi2.", g)]] <- mmi2_list[[g]]
      }

      needs_n_mm <- any(sapply(block_indices, function(i) {
        b <- mm_blocks[[i]]
        w_in_jags(b) ||
          (f_in_jags(b) && b$fn$type %in% member_loop_types) ||
          !is.null(b$ar)
      }))
      if (needs_n_mm) {
        jags.data.list[[paste0("n.mm.", g)]] <- n.mm_list[[g]]
      }

      # RE / FE need mmid indexing
      has_re_in_group <- any(sapply(block_indices, function(i) {
        !is.null(mm_blocks[[i]]$RE) || !is.null(mm_blocks[[i]]$FE)
      }))
      if (has_re_in_group) {
        jags.data.list[[paste0("mmid.", g)]] <- mmid_list[[g]]
        jags.data.list[[paste0("n.umm.", g)]] <- n.umm_list[[g]]
      }

      needs_grp <- any(sapply(block_indices, function(i) {
        b <- mm_blocks[[i]]
        (w_in_jags(b) && b$w$constraint) ||
          (f_in_jags(b) && b$fn$type %in% needs_grp_types)
      }))
      if (needs_grp) {
        jags.data.list[[paste0("grp.mm.", g)]] <- grp.mm_list[[g]]
      }

      has_ar_in_group <- any(sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$ar)))
      if (has_ar_in_group) {
        jags.data.list[[paste0("n.GPn.", g)]] <- n.GPn_list[[g]]
        jags.data.list[[paste0("n.GPNi.", g)]] <- n.GPNi_list[[g]]
        if (!is.null(gap.mm_list[[g]])) {
          jags.data.list[[paste0("gap.mm.", g)]] <- gap.mm_list[[g]]
        }
      }

      # RE slope covariates
      if (!is.null(X.re[[g]])) {
        jags.data.list[[paste0("X.re.", g)]] <- X.re[[g]]
      }
    }

    for (i in seq_along(mm_blocks)) {
      block <- mm_blocks[[i]]

      # Attribute design matrix: needed only when the feature is computed in JAGS
      if (n.Xmm[[i]] > 0 && f_in_jags(block)) {
        jags.data.list[[paste0("X.mm.", i)]] <- X.mm[[i]]
        if (block$fn$type == "sum" && length(block$feature_labels) > 1) {
          jags.data.list[[paste0("n.Xmm.", i)]] <- n.Xmm[[i]]
        }
      }

      # Precomputed features
      if (F.is.precomp[[i]]) {
        jags.data.list[[paste0("F.", i)]] <- F.precomp[[i]]
      }

      # Fixed-coefficient offsets
      if (!is.null(offset.mm[[i]])) {
        jags.data.list[[paste0("offset.mm.", i)]] <- offset.mm[[i]]
      }

      # Weights: precomputed as data, or covariates for in-JAGS computation
      if (w.is.precomp[[i]]) {
        # w.k needed as data whenever the model string references it
        needs_w <- f_in_jags(block) || !is.null(block$RE) || !is.null(block$FE) ||
          !is.null(offset.mm[[i]])
        if (needs_w) {
          jags.data.list[[paste0("w.", i)]] <- w.precomp[[i]]
        }
      } else if (n.Xw[[i]] > 0) {
        jags.data.list[[paste0("X.w.", i)]] <- X.w[[i]]
      }

      # FE slope covariates
      if (!is.null(X.fe[[i]])) {
        jags.data.list[[paste0("X.fe.", i)]] <- X.fe[[i]]
      }
    }
  }

  if (has_hm) {
    jags.data.list$hmid <- hmid
    jags.data.list$n.hm <- n.hm

    for (i in seq_along(hm_blocks)) {
      if (!is.null(X.re.hm[[i]])) {
        jags.data.list[[paste0("X.re.hm.", i)]] <- X.re.hm[[i]]
      }
    }

    if (any(sapply(hm_blocks, function(b) !is.null(b$ar)))) {
      jags.data.list$n.HMn <- n.HMn
      jags.data.list$n.HMNi <- n.HMNi
      for (i in seq_along(hm_blocks)) {
        if (!is.null(gap.hm_list[[i]])) {
          jags.data.list[[paste0("gap.hm.", i)]] <- gap.hm_list[[i]]
        }
      }
    }
  }

  # ========================================================================================== #
  # Return
  # ========================================================================================== #

  return(
    list(
      ids = list(
        mmid = mmid, mainid = mainid, hmid = hmid,
        mmi1 = mmi1, mmi2 = mmi2,
        mmid_list = mmid_list, mmi1_list = mmi1_list, mmi2_list = mmi2_list,
        mmn_list = mmn_list, grp.mm_list = grp.mm_list
      ),
      Ns = list(
        n.mm = n.mm, mmn = mmn, n.umm = n.umm,
        n.main = n.main, n.hm = n.hm,
        n.Xmm = n.Xmm, n.Xmain = n.Xmain, n.Xw = n.Xw,
        n.mmblocks = n.mmblocks,
        n.GPN = n.GPN, n.GPNi = n.GPNi, n.GPn = n.GPn,
        n.HMN = n.HMN, n.HMNi = n.HMNi, n.HMn = n.HMn,
        n.mm_list = n.mm_list, n.umm_list = n.umm_list, mmn_list = mmn_list,
        n.GPN_list = n.GPN_list, n.GPNi_list = n.GPNi_list, n.GPn_list = n.GPn_list,
        n_mmid_groups = n_mmid_groups,
        all_mmid_names = all_mmid_names, mmid_to_blocks = mmid_to_blocks
      ),
      Xs = list(
        X.mm = X.mm, X.main = X.main, X.w = X.w, X.re = X.re, X.re.hm = X.re.hm,
        w.is.precomp = w.is.precomp, w.precomp = w.precomp,
        F.is.precomp = F.is.precomp, F.precomp = F.precomp
      ),
      Ys = Ys,
      jags.params = unique(jags.params),
      jags.inits = jags.inits,
      jags.data = jags.data.list
    )
  )
}
