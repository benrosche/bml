# ================================================================================================ #
# Function createJagsVars
# ================================================================================================ #

createJagsVars <- function(data, family, mm_blocks, main, hm_blocks, mm, hm, monitor, modelfile, chains, inits) {

  # Unpack main ------------------------------------------------------------------------------ #

  lhs      <- main$lhs
  mainvars <- main$vars
  maindat  <- main$dat

  # Flags -------------------------------------------------------------------------------------- #

  has_mm <- !is.null(mm_blocks) && length(mm_blocks) > 0
  has_hm <- !is.null(hm_blocks) && length(hm_blocks) > 0
  has_mm_RE <- has_mm && attr(mm_blocks, "has_RE")
  has_mm_vars <- has_mm && attr(mm_blocks, "has_vars")

  # ========================================================================================== #
  # Create IDs
  # ========================================================================================== #

  mmid   <- if (has_mm) data %>% dplyr::arrange(mainid, mmid) %>% dplyr::pull(mmid) else c()
  mainid <- maindat %>% dplyr::pull(mainid)
  hmid   <- if (has_hm) maindat %>% dplyr::pull(hmid) else c()

  mmn <- if (has_mm) maindat %>% dplyr::pull(mmn) else c()

  mmi1 <- if (has_mm) maindat %>% dplyr::pull(mmi1) else c()
  mmi1.mm <- if (has_mm) rep(mmi1, mmn) else c()

  mmi2 <- if (has_mm) maindat %>% dplyr::pull(mmi2) else c()
  mmi2.mm <- if (has_mm) rep(mmi2, mmn) else c()

  # ========================================================================================== #
  # Create Ns
  # ========================================================================================== #

  n.mm   <- length(mmid)
  n.umm  <- length(unique(mmid))
  n.main <- length(mainid)
  n.hm   <- if (has_hm) length(unique(hmid)) else 0

  # For autoregressive RE - mm level
  n.GPN  <- if (has_mm) data %>% dplyr::count(mmid) %>% dplyr::pull(n) %>% max() else c()
  n.GPNi <- if (has_mm) data %>% dplyr::count(mmid) %>% dplyr::pull(n) else c()
  n.GPn  <- if (has_mm) data %>% dplyr::group_by(mmid) %>% dplyr::mutate(n = row_number()) %>% dplyr::pull(n) else c()

  # For autoregressive RE - hm level
  n.HMN  <- if (has_hm) maindat %>% dplyr::count(hmid) %>% dplyr::pull(n) %>% max() else c()
  n.HMNi <- if (has_hm) maindat %>% dplyr::count(hmid) %>% dplyr::pull(n) else c()
  n.HMn  <- if (has_hm) maindat %>% dplyr::group_by(hmid) %>% dplyr::mutate(n = row_number()) %>% dplyr::pull(n) else c()

  # Number of mm blocks
  n.mmblocks <- if (has_mm) length(mm_blocks) else 0

  # ========================================================================================== #
  # Create X matrices for mm() blocks
  # ========================================================================================== #

  X.mm      <- list()
  X.mm.fix  <- list()
  fix.mm    <- list()
  X.w       <- list()
  n.Xmm     <- list()
  n.Xmm.fix <- list()
  n.Xw      <- list()

  if (has_mm) {
    for (i in seq_along(mm_blocks)) {
      block <- mm_blocks[[i]]

      # Free variables matrix for this block
      if (!is.null(block$vars) && length(block$vars) > 0) {
        X.mm[[i]] <- block$dat %>% dplyr::select(all_of(block$vars)) %>% as.matrix()
        n.Xmm[[i]] <- ncol(X.mm[[i]])
      } else {
        X.mm[[i]] <- NULL
        n.Xmm[[i]] <- 0
      }

      # Fixed variables matrix for this block
      if (!is.null(block$dat_fixed)) {
        var_names <- sapply(block$vars_fixed, function(x) x$var)
        X.mm.fix[[i]] <- block$dat_fixed %>% dplyr::select(all_of(var_names)) %>% as.matrix()
        fix.mm[[i]] <- block$fix_values
        n.Xmm.fix[[i]] <- length(fix.mm[[i]])
      } else {
        X.mm.fix[[i]] <- NULL
        fix.mm[[i]] <- NULL
        n.Xmm.fix[[i]] <- 0
      }

      # Weight matrix for this block
      wvars <- block$fn$vars
      if (length(wvars) > 0) {
        X.w[[i]] <- block$wdat %>% dplyr::select(all_of(wvars)) %>% as.matrix()
        n.Xw[[i]] <- ncol(X.w[[i]])
      } else {
        X.w[[i]] <- NULL
        n.Xw[[i]] <- 0
      }
    }
  }

  # ========================================================================================== #
  # Create X matrices for main and hm levels
  # ========================================================================================== #

  # Main level - free variables
  mainvars_clean <- mainvars[mainvars != "X0"]
  if (length(mainvars_clean) > 0) {
    X.main <- maindat %>% dplyr::select(all_of(mainvars_clean)) %>% as.matrix()
    n.Xmain <- ncol(X.main)
  } else {
    X.main <- NULL
    n.Xmain <- 0
  }

  # Main level - fixed variables
  if (!is.null(main$dat_fixed)) {
    var_names <- sapply(main$vars_fixed, function(x) x$var)
    # Handle X0 specially
    var_names_data <- if ("X0" %in% var_names) c("X0", var_names[var_names != "X0"]) else var_names
    X.main.fix <- main$dat_fixed %>% dplyr::select(all_of(var_names_data)) %>% as.matrix()
    fix.main <- main$fix_values
    n.Xmain.fix <- length(fix.main)
  } else {
    X.main.fix <- NULL
    fix.main <- NULL
    n.Xmain.fix <- 0
  }

  # HM level
  X.hm      <- list()
  X.hm.fix  <- list()
  fix.hm    <- list()
  n.Xhm     <- list()
  n.Xhm.fix <- list()

  if (has_hm) {
    for (i in seq_along(hm_blocks)) {
      block <- hm_blocks[[i]]

      # Free variables
      if (!is.null(block$vars) && length(block$vars) > 0) {
        X.hm[[i]] <- block$dat %>% dplyr::select(all_of(block$vars)) %>% as.matrix()
        n.Xhm[[i]] <- ncol(X.hm[[i]])
      } else {
        X.hm[[i]] <- NULL
        n.Xhm[[i]] <- 0
      }

      # Fixed variables
      if (!is.null(block$dat_fixed)) {
        var_names <- sapply(block$vars_fixed, function(x) x$var)
        X.hm.fix[[i]] <- block$dat_fixed %>% dplyr::select(all_of(var_names)) %>% as.matrix()
        fix.hm[[i]] <- block$fix_values
        n.Xhm.fix[[i]] <- length(fix.hm[[i]])
      } else {
        X.hm.fix[[i]] <- NULL
        fix.hm[[i]] <- NULL
        n.Xhm.fix[[i]] <- 0
      }
    }
  }

  # ========================================================================================== #
  # Build jags.params
  # ========================================================================================== #

  jags.params <- c()

  # MM-level parameters
  if (has_mm_RE) {
    jags.params <- c(jags.params, "sigma.mm")
    if (monitor) jags.params <- c(jags.params, "re.mm")
  }

  # b.mm for each mm block with variables
  if (has_mm) {
    for (i in seq_along(mm_blocks)) {
      if (!is.null(mm_blocks[[i]]$vars) && length(mm_blocks[[i]]$vars) > 0) {
        jags.params <- c(jags.params, paste0("b.mm.", i))
      }
    }
  }

  # Main-level parameters
  if (family != "Cox") {
    jags.params <- c(jags.params, "sigma")
  }
  if (n.Xmain > 0) {
    jags.params <- c(jags.params, "b")
  }
  if (monitor) {
    jags.params <- c(jags.params, "pred")
  }

  # HM-level parameters
  if (has_hm) {
    for (i in seq_along(hm_blocks)) {
      block <- hm_blocks[[i]]
      if (block$type == "RE") {
        jags.params <- c(jags.params, paste0("sigma.hm.", i))
        if (monitor) jags.params <- c(jags.params, paste0("re.hm.", i))
      }
      if (n.Xhm[[i]] > 0 && (block$type == "RE" || block$showFE)) {
        jags.params <- c(jags.params, paste0("b.hm.", i))
      }
    }
  }

  # Weight parameters for each mm block
  if (has_mm) {
    for (i in seq_along(mm_blocks)) {
      block <- mm_blocks[[i]]
      if (length(block$fn$params) > 0) {
        jags.params <- c(jags.params, paste0("b.w.", i))
      }
      if (monitor) {
        jags.params <- c(jags.params, paste0("w.", i))
      }
    }
  }

  # ========================================================================================== #
  # Build jags.data
  # ========================================================================================== #

  jags.data <- c()

  # MM-level data
  if (has_mm) {
    jags.data <- c(jags.data, "mmid", "mmi1", "mmi2", "n.mm", "n.umm", "n.mmblocks")

    for (i in seq_along(mm_blocks)) {
      block <- mm_blocks[[i]]

      # X.mm for this block
      if (n.Xmm[[i]] > 0) {
        jags.data <- c(jags.data, paste0("X.mm.", i), paste0("n.Xmm.", i))
      }

      # X.w for this block
      if (n.Xw[[i]] > 0) {
        jags.data <- c(jags.data, paste0("X.w.", i))
      }

      # Constraint indices if needed
      if (block$fn$constraint) {
        jags.data <- c(jags.data, "mmi1.mm", "mmi2.mm")
      }
    }

    # AR data if any block uses it
    if (any(sapply(mm_blocks, function(b) b$ar))) {
      jags.data <- c(jags.data, "n.GPn", "n.GPNi")
    }
  }

  # Main-level data
  jags.data <- c(jags.data, "n.main")
  if (n.Xmain > 0) {
    jags.data <- c(jags.data, "X.main", "n.Xmain")
  }

  # HM-level data
  if (has_hm) {
    jags.data <- c(jags.data, "hmid", "n.hm")
    for (i in seq_along(hm_blocks)) {
      if (n.Xhm[[i]] > 0) {
        jags.data <- c(jags.data, paste0("X.hm.", i), paste0("n.Xhm.", i))
      }
    }

    # AR data for hm blocks
    if (any(sapply(hm_blocks, function(b) b$ar))) {
      jags.data <- c(jags.data, "n.HMn", "n.HMNi")
    }
  }

  # Remove duplicates
  jags.data <- unique(jags.data)

  # ========================================================================================== #
  # Family-specific data and inits
  # ========================================================================================== #

  if (family %in% c("Gaussian", "Binomial")) {

    Y <- maindat %>% dplyr::rename(Y = all_of(lhs)) %>% dplyr::pull(Y)
    jags.data <- c(jags.data, "Y")
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

    jags.data <- c(jags.data, "t", "ct.lb", "censored")
    jags.params <- c(jags.params, "shape")

    t.init <- t
    t.init[] <- NA
    t.init[censored == 1] <- ct.lb[censored == 1] + 1

    jags.inits <- list(t = t.init, shape = 1)
    Ys <- list(t = t, ct.lb = ct.lb, event = event, censored = censored)

  } else if (family == "Cox") {

    t <- maindat %>% dplyr::rename(t = all_of(lhs[1]), ev = all_of(lhs[2])) %>% dplyr::pull(t)
    t.unique <- c(sort(unique(t)), max(t) + 1)
    n.tu <- length(t.unique) - 1
    event <- maindat %>% dplyr::rename(ev = all_of(lhs[2])) %>% dplyr::pull(ev)

    Y <- matrix(data = NA, nrow = n.main, ncol = n.tu)
    dN <- matrix(data = NA, nrow = n.main, ncol = n.tu)
    for (j in 1:n.main) {
      for (k in 1:n.tu) {
        Y[j, k] <- as.numeric(t[j] - t.unique[k] + 1e-05 >= 0)
        dN[j, k] <- Y[j, k] * event[j] * as.numeric(t.unique[k + 1] - t[j] >= 1e-05)
      }
    }

    jags.data <- c(jags.data, "Y", "dN", "t.unique", "n.tu", "c", "d")
    jags.inits <- list(dL0 = rep(1.0, n.tu))
    Ys <- list(Y = Y, dN = dN, t = t, t.unique = t.unique, event = event, c = 0.001, d = 0.1, n.tu = n.tu)

  }

  # ========================================================================================== #
  # Finalize inits
  # ========================================================================================== #

  jags.inits <- c(jags.inits, inits)
  jags.inits <- lapply(1:chains, function(x) jags.inits)

  # ========================================================================================== #
  # Build the actual data list for JAGS
  # ========================================================================================== #

  # Start with scalars and basic vectors
  jags.data.list <- list(
    n.main = n.main
  )

  # Add Y or survival data
  if (family %in% c("Gaussian", "Binomial")) {
    jags.data.list$Y <- Ys$Y
  } else if (family == "Weibull") {
    jags.data.list$t <- Ys$t
    jags.data.list$ct.lb <- Ys$ct.lb
    jags.data.list$censored <- Ys$censored
  } else if (family == "Cox") {
    jags.data.list$Y <- Ys$Y
    jags.data.list$dN <- Ys$dN
    jags.data.list$t.unique <- Ys$t.unique
    jags.data.list$n.tu <- Ys$n.tu
    jags.data.list$c <- Ys$c
    jags.data.list$d <- Ys$d
  }

  # Main-level data
  if (n.Xmain > 0) {
    jags.data.list$X.main <- X.main
    # n.Xmain not passed to JAGS (unused in model; b range is hardcoded)
  }

  # Main-level fixed data
  if (n.Xmain.fix > 0) {
    jags.data.list$X.fix.main <- X.main.fix
    jags.data.list$fix.main <- fix.main
  }

  # MM-level data
  if (has_mm) {
    jags.data.list$mmi1 <- mmi1
    jags.data.list$mmi2 <- mmi2
    jags.data.list$n.mm <- n.mm
    # n.mmblocks not passed to JAGS (unused in model; available in R via Ns$n.mmblocks)

    # Only pass mmid and n.umm if we have random effects
    if (has_mm_RE) {
      jags.data.list$mmid <- mmid
      jags.data.list$n.umm <- n.umm
    }

    # Per-block data
    for (i in seq_along(mm_blocks)) {
      block <- mm_blocks[[i]]

      # Free variables
      if (n.Xmm[[i]] > 0) {
        jags.data.list[[paste0("X.mm.", i)]] <- X.mm[[i]]
        jags.data.list[[paste0("n.Xmm.", i)]] <- n.Xmm[[i]]
      }

      # Fixed variables
      if (n.Xmm.fix[[i]] > 0) {
        jags.data.list[[paste0("X.fix.mm.", i)]] <- X.mm.fix[[i]]
        jags.data.list[[paste0("fix.mm.", i)]] <- fix.mm[[i]]
      }

      # Weight variables
      if (n.Xw[[i]] > 0) {
        jags.data.list[[paste0("X.w.", i)]] <- X.w[[i]]
      }
    }

    # Constraint indices
    if (any(sapply(mm_blocks, function(b) b$fn$constraint))) {
      jags.data.list$mmi1.mm <- mmi1.mm
      jags.data.list$mmi2.mm <- mmi2.mm
    }

    # AR data
    if (any(sapply(mm_blocks, function(b) b$ar))) {
      jags.data.list$n.GPn <- n.GPn
      jags.data.list$n.GPNi <- n.GPNi
    }
  }

  # HM-level data
  if (has_hm) {
    jags.data.list$hmid <- hmid
    jags.data.list$n.hm <- n.hm

    for (i in seq_along(hm_blocks)) {
      # Free variables
      if (n.Xhm[[i]] > 0) {
        jags.data.list[[paste0("X.hm.", i)]] <- X.hm[[i]]
        jags.data.list[[paste0("n.Xhm.", i)]] <- n.Xhm[[i]]
      }

      # Fixed variables
      if (n.Xhm.fix[[i]] > 0) {
        jags.data.list[[paste0("X.fix.hm.", i)]] <- X.hm.fix[[i]]
        jags.data.list[[paste0("fix.hm.", i)]] <- fix.hm[[i]]
      }
    }

    # AR data for hm blocks
    if (any(sapply(hm_blocks, function(b) b$ar))) {
      jags.data.list$n.HMn <- n.HMn
      jags.data.list$n.HMNi <- n.HMNi
    }
  }

  # ========================================================================================== #
  # Return
  # ========================================================================================== #

  return(
    list(
      ids = list(
        mmid = mmid, mainid = mainid, hmid = hmid,
        mmi1 = mmi1, mmi1.mm = mmi1.mm,
        mmi2 = mmi2, mmi2.mm = mmi2.mm
      ),
      Ns = list(
        n.mm = n.mm, mmn = mmn, n.umm = n.umm,
        n.main = n.main, n.hm = n.hm,
        n.Xmm = n.Xmm, n.Xmain = n.Xmain, n.Xhm = n.Xhm, n.Xw = n.Xw,
        n.Xmm.fix = n.Xmm.fix, n.Xmain.fix = n.Xmain.fix, n.Xhm.fix = n.Xhm.fix,
        n.mmblocks = n.mmblocks,
        n.GPN = n.GPN, n.GPNi = n.GPNi, n.GPn = n.GPn,
        n.HMN = n.HMN, n.HMNi = n.HMNi, n.HMn = n.HMn
      ),
      Xs = list(
        X.mm = X.mm, X.main = X.main, X.hm = X.hm, X.w = X.w,
        X.mm.fix = X.mm.fix, X.main.fix = X.main.fix, X.hm.fix = X.hm.fix,
        fix.mm = fix.mm, fix.main = fix.main, fix.hm = fix.hm
      ),
      Ys = Ys,
      jags.params = jags.params,
      jags.inits = jags.inits,
      jags.data = jags.data.list
    )
  )

}
