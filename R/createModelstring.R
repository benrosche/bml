# ================================================================================================ #
# Function createModelstring
# ================================================================================================ #
#
# Programmatically constructs the complete JAGS model string.
#   - MM level: weight functions, group-level feature nodes F.k (one per block), member REs/FEs
#   - HM level: unit random/fixed effects (re()/fe() grammar)
#   - Main level: linear predictor (features enter with coefficients b.fn.k; named-feature
#     interactions with b.int), likelihood, log_lik for loo/waic
#   - Priors: emitted from the resolved prior table (prior() system); raw JAGS strings are
#     applied as regex replacements at the end (escape hatch)
#
# ================================================================================================ #

createModelstring <- function(family, prior_spec, mm_blocks, main, hm_blocks, interactions,
                              monitor_spec, cox_intervals = NULL) {

  overrides <- prior_spec$overrides %||% list()
  raw_priors <- prior_spec$raw %||% character(0)

  # ========================================================================================== #
  # Prior helpers
  # ========================================================================================== #

  # Distribution for a node: user override or default
  pdist <- function(node, default) {
    overrides[[node]] %||% default
  }

  # Emit an sd-type prior: default is precision-scale dscaled.gamma; an override is a
  # distribution on the SD itself, with the precision derived.
  sd_prior_lines <- function(sigma_node, tau_node) {
    ov <- overrides[[sigma_node]]
    if (is.null(ov)) {
      c(paste0("  ", tau_node, " ~ dscaled.gamma(25, 1)"),
        paste0("  ", sigma_node, " <- 1/sqrt(", tau_node, ")"))
    } else {
      c(paste0("  ", sigma_node, " ~ ", ov),
        paste0("  ", tau_node, " <- pow(", sigma_node, ", -2)"))
    }
  }

  # Direct-sigma prior (for the correlated-RE separation strategy): SD sampled directly.
  sd_direct_line <- function(sigma_node) {
    ov <- overrides[[sigma_node]]
    if (is.null(ov)) {
      # half-Cauchy(25): same marginal as dscaled.gamma(25, 1)
      paste0("  ", sigma_node, " ~ dt(0, ", format(1 / 25^2), ", 1)T(0,)")
    } else {
      paste0("  ", sigma_node, " ~ ", ov)
    }
  }

  # Correlation prior: rho = 2*z - 1, z ~ dbeta(eta, eta); lkj(eta) maps onto eta.
  cor_prior_lines <- function(rho_node) {
    ov <- overrides[[rho_node]] %||% "__lkj__(1)"
    eta <- sub("^__lkj__\\((.*)\\)$", "\\1", ov)
    if (identical(eta, ov)) {
      # not an lkj spec: user gave a direct distribution on rho
      paste0("  ", rho_node, " ~ ", ov)
    } else {
      c(paste0("  z.", rho_node, " ~ dbeta(", eta, ", ", eta, ")"),
        paste0("  ", rho_node, " <- 2 * z.", rho_node, " - 1"))
    }
  }

  # ========================================================================================== #
  # Flags
  # ========================================================================================== #

  has_mm <- !is.null(mm_blocks) && length(mm_blocks) > 0
  has_hm <- !is.null(hm_blocks) && length(hm_blocks) > 0
  has_int <- length(interactions %||% list()) > 0

  mainvars <- main$vars
  lhs <- main$lhs

  # Per-block computation mode
  w_in_jags <- function(block) length(block$w$params) > 0
  # Features must be computed in JAGS when the weights are parametric or fn has estimated params
  f_in_jags <- function(block) {
    (length(block$feature_labels) > 0) &&
      (w_in_jags(block) || length(block$fn$est_params) > 0)
  }
  has_features <- function(block) length(block$feature_labels) > 0

  # Shape-parameter reference: literal for fixed values, node for estimated
  fn_param_ref <- function(block, k, pname) {
    if (pname %in% names(block$fn$est_params)) {
      paste0("fn.", pname, ".", k)
    } else {
      format(block$fn$fixed_params[[pname]])
    }
  }

  # ========================================================================================== #
  # Build model string programmatically
  # ========================================================================================== #

  lines <- c()
  add <- function(...) lines <<- c(lines, paste0(...))

  add("model {")
  add("")

  # ------------------------------------------------------------------------------------------ #
  # MM level
  # ------------------------------------------------------------------------------------------ #

  if (has_mm) {

    add("  # ==================== MM Level: Multiple Membership ==================== #")
    add("")

    all_mmid_names <- attr(mm_blocks, "all_mmid_names")
    mmid_to_blocks <- attr(mm_blocks, "mmid_to_blocks")

    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]
      g <- block$mmid_group
      rng <- function(v) paste0(v, ".", g, "[j]")
      idx_range <- paste0("mmi1.", g, "[j]:mmi2.", g, "[j]")

      # ------------------------- Weights ------------------------- #
      if (!w_in_jags(block)) {
        add("  # Weights for mm block ", k, " (pre-computed in R, passed as data)")
      } else {
        add("  # Weight function for mm block ", k, " (mmid group ", g, ")")

        fn_string <- block$w$string
        for (p in seq_along(block$w$params)) {
          fn_string <- gsub(paste0("\\b", block$w$params[p], "\\b"),
                            paste0("b.w.", k, "[", p, "]"), fn_string)
        }
        for (v in seq_along(block$w$vars)) {
          fn_string <- gsub(paste0("\\b", block$w$vars[v], "\\b"),
                            paste0("X.w.", k, "[i,", v, "]"), fn_string)
        }

        if (block$w$constraint) {
          add("  for (i in 1:n.mm.", g, ") {")
          add("    uw.", k, "[i] <- ", fn_string)
          add("  }")
          add("")
          add("  # Accumulator: cumulative sums for efficient group sums")
          add("  cum.uw.", k, "[1] <- 0")
          add("  for (i in 1:n.mm.", g, ") {")
          add("    cum.uw.", k, "[i+1] <- cum.uw.", k, "[i] + uw.", k, "[i]")
          add("  }")
          add("")
          add("  for (j in 1:n.main) {")
          add("    sum.uw.", k, "[j] <- cum.uw.", k, "[mmi2.", g, "[j]+1] - cum.uw.", k, "[mmi1.", g, "[j]]")
          add("  }")
          add("")
          add("  for (i in 1:n.mm.", g, ") {")
          add("    w.", k, "[i] <- uw.", k, "[i] / sum.uw.", k, "[grp.mm.", g, "[i]]")
          add("  }")
        } else {
          add("  for (i in 1:n.mm.", g, ") {")
          add("    w.", k, "[i] <- ", fn_string)
          add("  }")
        }
      }
      add("")

      # ------------------------- Features ------------------------- #
      if (has_features(block)) {
        n_feat <- length(block$feature_labels)

        if (!f_in_jags(block)) {
          add("  # Feature(s) for mm block ", k, " [",
              paste(block$feature_labels, collapse = ", "), "] (pre-computed in R)")
        } else {
          add("  # Feature(s) for mm block ", k, " [",
              paste(block$feature_labels, collapse = ", "), "]")
          type <- block$fn$type
          xref <- function(col = 1) paste0("X.mm.", k, "[i,", col, "]")
          xref_j <- function(col = 1) paste0("X.mm.", k, "[", idx_range, ",", col, "]")
          w_i <- paste0("w.", k, "[i]")
          w_j <- paste0("w.", k, "[", idx_range, "]")
          grp_i <- paste0("grp.mm.", g, "[i]")

          if (type == "sum") {
            add("  for (j in 1:n.main) {")
            if (n_feat == 1) {
              add("    F.", k, "[j] <- inprod(", w_j, ", ", xref_j(1), ")")
            } else {
              add("    for (x in 1:n.Xmm.", k, ") {")
              add("      F.", k, "[j,x] <- inprod(", w_j, ", X.mm.", k, "[", idx_range, ", x])")
              add("    }")
            }
            add("  }")

          } else if (type == "var") {
            m <- block$fn$moment %||% 2
            add("  for (j in 1:n.main) {")
            add("    A.", k, "[j] <- inprod(", w_j, ", ", xref_j(1), ")")
            add("  }")
            add("  for (i in 1:n.mm.", g, ") {")
            add("    v.", k, "[i] <- ", w_i, " * pow(", xref(1), " - A.", k, "[", grp_i, "], ", m, ")")
            add("  }")
            add("  for (j in 1:n.main) {")
            add("    F.", k, "[j] <- sum(v.", k, "[", idx_range, "])")
            add("  }")

          } else if (type %in% c("hhi", "effn")) {
            add("  for (j in 1:n.main) {")
            if (type == "hhi") {
              add("    F.", k, "[j] <- inprod(", w_j, ", ", w_j, ")")
            } else {
              add("    F.", k, "[j] <- 1 / inprod(", w_j, ", ", w_j, ")")
            }
            add("  }")

          } else if (type == "entropy") {
            add("  for (i in 1:n.mm.", g, ") {")
            add("    e.", k, "[i] <- ", w_i, " * log(", w_i, ")")
            add("  }")
            add("  for (j in 1:n.main) {")
            add("    F.", k, "[j] <- -sum(e.", k, "[", idx_range, "])")
            add("  }")

          } else if (type == "threshold") {
            cref <- fn_param_ref(block, k, "c")
            kref <- fn_param_ref(block, k, "kappa")
            add("  for (i in 1:n.mm.", g, ") {")
            add("    s.", k, "[i] <- ", w_i, " * ilogit(", kref, " * (", xref(1), " - ", cref, "))")
            add("  }")
            add("  for (j in 1:n.main) {")
            add("    F.", k, "[j] <- sum(s.", k, "[", idx_range, "])")
            add("  }")

          } else if (type == "smax") {
            kref <- fn_param_ref(block, k, "kappa")
            add("  for (i in 1:n.mm.", g, ") {")
            add("    e.", k, "[i] <- ", w_i, " * exp(", kref, " * ", xref(1), ")")
            add("  }")
            add("  for (j in 1:n.main) {")
            add("    F.", k, "[j] <- (1 / ", kref, ") * log(sum(e.", k, "[", idx_range, "]))")
            add("  }")

          } else if (type == "gmean") {
            pref <- fn_param_ref(block, k, "p")
            add("  for (i in 1:n.mm.", g, ") {")
            add("    e.", k, "[i] <- ", w_i, " * pow(", xref(1), ", ", pref, ")")
            add("  }")
            add("  for (j in 1:n.main) {")
            add("    F.", k, "[j] <- pow(sum(e.", k, "[", idx_range, "]), 1 / ", pref, ")")
            add("  }")

          } else if (type == "cov") {
            add("  for (j in 1:n.main) {")
            add("    A1.", k, "[j] <- inprod(", w_j, ", ", xref_j(1), ")")
            add("    A2.", k, "[j] <- inprod(", w_j, ", ", xref_j(2), ")")
            add("  }")
            add("  for (i in 1:n.mm.", g, ") {")
            add("    v.", k, "[i] <- ", w_i, " * (", xref(1), " - A1.", k, "[", grp_i, "]) * (",
                xref(2), " - A2.", k, "[", grp_i, "])")
            add("  }")
            add("  for (j in 1:n.main) {")
            add("    F.", k, "[j] <- sum(v.", k, "[", idx_range, "])")
            add("  }")

          } else if (type == "expr") {
            # DSL: one accumulator pass per E() node (inner nodes first)
            graph <- block$fn$graph
            subs <- c()
            for (a in seq_along(block$attr_cols)) {
              col <- match(block$attr_cols[a], block$vars)
              subs[block$attr_cols[a]] <- paste0("X.mm.", k, "[i,", col, "]")
            }
            subs["w"] <- w_i
            for (pn in names(block$fn$est_params)) {
              subs[pn] <- paste0("fn.", pn, ".", k)
            }
            for (pn in names(block$fn$fixed_params %||% list())) {
              subs[pn] <- format(block$fn$fixed_params[[pn]])
            }
            for (e in seq_along(graph$enodes)) {
              subs[paste0(".E", e)] <- paste0("E", e, ".", k, "[", grp_i, "]")
            }

            for (e in seq_along(graph$enodes)) {
              body <- dsl_deparse_jags(graph$enodes[[e]]$body, subs)
              add("  for (i in 1:n.mm.", g, ") {")
              add("    t", e, ".", k, "[i] <- ", w_i, " * ", body)
              add("  }")
              add("  for (j in 1:n.main) {")
              add("    E", e, ".", k, "[j] <- sum(t", e, ".", k, "[", idx_range, "])")
              add("  }")
            }
            # Top level: group scalar
            subs_top <- subs
            for (e in seq_along(graph$enodes)) {
              subs_top[paste0(".E", e)] <- paste0("E", e, ".", k, "[j]")
            }
            top <- dsl_deparse_jags(graph$top, subs_top)
            add("  for (j in 1:n.main) {")
            add("    F.", k, "[j] <- ", top)
            add("  }")
          }
        }
        add("")

        # Priors for feature coefficients
        add("  # Prior(s) for the feature coefficient(s) of mm block ", k)
        for (x in seq_len(n_feat)) {
          node <- paste0("b.fn.", k, "[", x, "]")
          add("  ", node, " ~ ", pdist(node, "dnorm(0, 0.0001)"))
        }
        add("")
      }

      # ------------------------- fn shape-parameter priors ------------------------- #
      if (length(block$fn$est_params) > 0) {
        add("  # Priors for fn shape parameters of mm block ", k)
        for (pn in names(block$fn$est_params)) {
          node <- paste0("fn.", pn, ".", k)
          add("  ", node, " ~ ", pdist(node, block$fn$est_params[[pn]]$default))
        }
        add("")
      }

      # ------------------------- Weight-parameter priors ------------------------- #
      if (length(block$w$params) > 0) {
        add("  # Priors for weight parameters of mm block ", k)
        for (p in seq_along(block$w$params)) {
          node <- paste0("b.w.", k, "[", p, "]")
          add("  ", node, " ~ ", pdist(node, "dnorm(0, 0.0001)"))
        }
        add("")
      }

      # ------------------------- Member fixed effects (fe) ------------------------- #
      if (!is.null(block$FE)) {
        fe <- block$FE
        add("  # Member-level fixed effects for mm block ", k, " (unit 1 = reference)")
        if (fe$intercept) {
          add("  alpha.mm.", k, "[1] <- 0")
          add("  for (u in 2:n.umm.", g, ") {")
          add("    alpha.mm.", k, "[u] ~ dnorm(0, 0.0001)")
          add("  }")
        }
        for (s in seq_along(fe$slopes)) {
          add("  alphas.mm.", k, ".s", s, "[1] <- 0")
          add("  for (u in 2:n.umm.", g, ") {")
          add("    alphas.mm.", k, ".s", s, "[u] ~ dnorm(0, 0.0001)")
          add("  }")
        }
        add("")
      }
    }

    # ------------------------- Random effects (per mmid group) ------------------------- #
    for (g in seq_along(all_mmid_names)) {
      block_indices <- mmid_to_blocks[[all_mmid_names[g]]]
      re_block_idx <- block_indices[sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$RE))]
      if (length(re_block_idx) == 0) next
      re_block <- mm_blocks[[re_block_idx[1]]]
      re_spec <- re_block$RE
      any_ar_in_group <- any(sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$ar)))
      n_slopes <- length(re_spec$slopes)

      add("  # MM-level random effects (mmid group ", g, ")")

      if (any_ar_in_group) {
        # AR random walk (intercept-only, validated upstream). For a time-indexed
        # walk the step precision is scaled by the normalized time gap:
        # u_t ~ N(u_{t-1}, sigma^2 * gap[i,t])  <=>  precision tau / gap[i,t]
        step_prec <- if (!is.null(re_spec$ar) && !is.null(re_spec$ar$time)) {
          paste0("tau.mm.", g, " / gap.mm.", g, "[i,t]")
        } else {
          paste0("tau.mm.", g)
        }
        add("  for (i in 1:n.umm.", g, ") {")
        add("    re.mm.", g, "[i,1] ~ dnorm(0, tau.mm.", g, ")")
        add("    for (t in 2:n.GPNi.", g, "[i]) {")
        add("      re.mm.", g, "[i,t] ~ dnorm(re.mm.", g, "[i,t-1], ", step_prec, ")")
        add("    }")
        add("  }")
        add("")
        add("  # Extract AR random effects for each mm observation (mmid group ", g, ")")
        add("  for (i in 1:n.mm.", g, ") {")
        add("    re.mm.", g, ".i[i] <- re.mm.", g, "[mmid.", g, "[i], n.GPn.", g, "[i]]")
        add("  }")
        add("")
        add(sd_prior_lines(paste0("sigma.mm.", g), paste0("tau.mm.", g)))

      } else if (isTRUE(re_spec$cor) && n_slopes == 1) {
        # Correlated intercept + slope (2x2, separation strategy)
        add("  zero.mm.", g, "[1] <- 0")
        add("  zero.mm.", g, "[2] <- 0")
        add("  Sigma.u.", g, "[1,1] <- pow(sigma.mm.", g, ", 2)")
        add("  Sigma.u.", g, "[1,2] <- rho.mm.", g, " * sigma.mm.", g, " * sigma.mm.", g, ".s1")
        add("  Sigma.u.", g, "[2,1] <- Sigma.u.", g, "[1,2]")
        add("  Sigma.u.", g, "[2,2] <- pow(sigma.mm.", g, ".s1, 2)")
        add("  for (i in 1:n.umm.", g, ") {")
        add("    u.mm.", g, "[i,1:2] ~ dmnorm.vcov(zero.mm.", g, ", Sigma.u.", g, ")")
        add("    re.mm.", g, "[i] <- u.mm.", g, "[i,1]")
        add("    re.mm.", g, ".s1[i] <- u.mm.", g, "[i,2]")
        add("  }")
        add(sd_direct_line(paste0("sigma.mm.", g)))
        add(sd_direct_line(paste0("sigma.mm.", g, ".s1")))
        add(cor_prior_lines(paste0("rho.mm.", g)))

      } else {
        # Independent intercept and/or slopes
        if (re_spec$intercept) {
          add("  for (i in 1:n.umm.", g, ") {")
          add("    re.mm.", g, "[i] ~ dnorm(0, tau.mm.", g, ")")
          add("  }")
          add(sd_prior_lines(paste0("sigma.mm.", g), paste0("tau.mm.", g)))
        }
        for (s in seq_len(n_slopes)) {
          add("  for (i in 1:n.umm.", g, ") {")
          add("    re.mm.", g, ".s", s, "[i] ~ dnorm(0, tau.mm.", g, ".s", s, ")")
          add("  }")
          add(sd_prior_lines(paste0("sigma.mm.", g, ".s", s), paste0("tau.mm.", g, ".s", s)))
        }
      }
      add("")
    }
  }

  # ------------------------------------------------------------------------------------------ #
  # HM level
  # ------------------------------------------------------------------------------------------ #

  if (has_hm) {

    add("  # ==================== HM Level: Hierarchical Membership ==================== #")
    add("")

    for (k in seq_along(hm_blocks)) {
      block <- hm_blocks[[k]]

      if (!is.null(block$RE)) {
        re_spec <- block$RE
        n_slopes <- length(re_spec$slopes)

        add("  # HM-level random effects (hm block ", k, ")")

        if (!is.null(block$ar)) {
          # Time-indexed walk: step precision scaled by the normalized time gap
          step_prec <- if (!is.null(block$ar$time)) {
            paste0("tau.hm.", k, " / gap.hm.", k, "[c,t]")
          } else {
            paste0("tau.hm.", k)
          }
          add("  for (c in 1:n.hm) {")
          add("    re.hm.", k, "[c,1] ~ dnorm(0, tau.hm.", k, ")")
          add("    for (t in 2:n.HMNi[c]) {")
          add("      re.hm.", k, "[c,t] ~ dnorm(re.hm.", k, "[c,t-1], ", step_prec, ")")
          add("    }")
          add("  }")
          add(sd_prior_lines(paste0("sigma.hm.", k), paste0("tau.hm.", k)))

        } else if (isTRUE(re_spec$cor) && n_slopes == 1) {
          add("  zero.hm.", k, "[1] <- 0")
          add("  zero.hm.", k, "[2] <- 0")
          add("  Sigma.u.hm.", k, "[1,1] <- pow(sigma.hm.", k, ", 2)")
          add("  Sigma.u.hm.", k, "[1,2] <- rho.hm.", k, " * sigma.hm.", k, " * sigma.hm.", k, ".s1")
          add("  Sigma.u.hm.", k, "[2,1] <- Sigma.u.hm.", k, "[1,2]")
          add("  Sigma.u.hm.", k, "[2,2] <- pow(sigma.hm.", k, ".s1, 2)")
          add("  for (c in 1:n.hm) {")
          add("    u.hm.", k, "[c,1:2] ~ dmnorm.vcov(zero.hm.", k, ", Sigma.u.hm.", k, ")")
          add("    re.hm.", k, "[c] <- u.hm.", k, "[c,1]")
          add("    re.hm.", k, ".s1[c] <- u.hm.", k, "[c,2]")
          add("  }")
          add(sd_direct_line(paste0("sigma.hm.", k)))
          add(sd_direct_line(paste0("sigma.hm.", k, ".s1")))
          add(cor_prior_lines(paste0("rho.hm.", k)))

        } else {
          if (re_spec$intercept) {
            add("  for (c in 1:n.hm) {")
            add("    re.hm.", k, "[c] ~ dnorm(0, tau.hm.", k, ")")
            add("  }")
            add(sd_prior_lines(paste0("sigma.hm.", k), paste0("tau.hm.", k)))
          }
          for (s in seq_len(n_slopes)) {
            add("  for (c in 1:n.hm) {")
            add("    re.hm.", k, ".s", s, "[c] ~ dnorm(0, tau.hm.", k, ".s", s, ")")
            add("  }")
            add(sd_prior_lines(paste0("sigma.hm.", k, ".s", s), paste0("tau.hm.", k, ".s", s)))
          }
        }
        add("")

      } else if (!is.null(block$FE)) {
        fe <- block$FE
        add("  # HM-level fixed effects (hm block ", k, "; unit 1 = reference)")
        if (fe$intercept) {
          add("  alpha.hm.", k, "[1] <- 0")
          add("  for (c in 2:n.hm) {")
          add("    alpha.hm.", k, "[c] ~ dnorm(0, 0.0001)")
          add("  }")
        }
        for (s in seq_along(fe$slopes)) {
          add("  alphas.hm.", k, ".s", s, "[1] <- 0")
          add("  for (c in 2:n.hm) {")
          add("    alphas.hm.", k, ".s", s, "[c] ~ dnorm(0, 0.0001)")
          add("  }")
        }
        add("")
      }
    }
  }

  # ------------------------------------------------------------------------------------------ #
  # Main level
  # ------------------------------------------------------------------------------------------ #

  add("  # ==================== Main Level: Main Model ==================== #")
  add("")

  add("  for (j in 1:n.main) {")
  add("")

  # Aggregate mm RE/FE/offset contributions (member-level parts that stay inside the sum)
  mm_agg_needed <- FALSE
  if (has_mm) {
    all_mmid_names <- attr(mm_blocks, "all_mmid_names")
    mmid_to_blocks <- attr(mm_blocks, "mmid_to_blocks")

    mm_terms_lines <- c()
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]
      g <- block$mmid_group
      idx_range <- paste0("mmi1.", g, "[j]:mmi2.", g, "[j]")
      block_indices <- mmid_to_blocks[[all_mmid_names[g]]]
      any_ar_in_group <- any(sapply(block_indices, function(i) !is.null(mm_blocks[[i]]$ar)))

      terms_k <- c()

      # fixed-coefficient offsets (fix() in vars)
      if (!is.null(block$vars_fixed) && length(block$vars_fixed) > 0) {
        terms_k <- c(terms_k, paste0("w.", k, "[", idx_range, "] * offset.mm.", k, "[", idx_range, "]"))
      }

      # random effects (weighted aggregation)
      if (!is.null(block$RE)) {
        re_spec <- block$RE
        if (any_ar_in_group) {
          terms_k <- c(terms_k, paste0("w.", k, "[", idx_range, "] * re.mm.", g, ".i[", idx_range, "]"))
        } else {
          if (re_spec$intercept) {
            terms_k <- c(terms_k, paste0("w.", k, "[", idx_range, "] * re.mm.", g,
                                         "[mmid.", g, "[", idx_range, "]]"))
          }
          for (s in seq_along(re_spec$slopes)) {
            terms_k <- c(terms_k, paste0("w.", k, "[", idx_range, "] * re.mm.", g, ".s", s,
                                         "[mmid.", g, "[", idx_range, "]] * X.re.", g,
                                         "[", idx_range, ",", s, "]"))
          }
        }
      }

      # member fixed effects
      if (!is.null(block$FE)) {
        fe <- block$FE
        if (fe$intercept) {
          terms_k <- c(terms_k, paste0("w.", k, "[", idx_range, "] * alpha.mm.", k,
                                       "[mmid.", g, "[", idx_range, "]]"))
        }
        for (s in seq_along(fe$slopes)) {
          terms_k <- c(terms_k, paste0("w.", k, "[", idx_range, "] * alphas.mm.", k, ".s", s,
                                       "[mmid.", g, "[", idx_range, "]] * X.fe.", k,
                                       "[", idx_range, ",", s, "]"))
        }
      }

      if (length(terms_k) > 0) {
        mm_terms_lines <- c(mm_terms_lines, paste0("      ", paste(terms_k, collapse = " + ")))
      }
    }

    if (length(mm_terms_lines) > 0) {
      mm_agg_needed <- TRUE
      add("    # Aggregate mm-level random/fixed member contributions")
      add("    mm.agg[j] <- sum(")
      add(paste(mm_terms_lines, collapse = " +\n"), ")")
      add("")
    }
  }

  # Build mu[j]
  mu_terms <- c()

  if (length(mainvars) > 0) {
    mu_terms <- c(mu_terms, "inprod(X.main[j,], b)")
  }
  if (!is.null(main$vars_fixed) && length(main$vars_fixed) > 0) {
    mu_terms <- c(mu_terms, "offset.main[j]")
  }

  # Block features with their coefficients
  if (has_mm) {
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]
      if (!has_features(block)) next
      n_feat <- length(block$feature_labels)
      if (n_feat == 1) {
        mu_terms <- c(mu_terms, paste0("b.fn.", k, "[1] * F.", k, "[j]"))
      } else {
        mu_terms <- c(mu_terms, paste0("inprod(F.", k, "[j,], b.fn.", k, ")"))
      }
    }
  }

  # Named-feature interactions
  if (has_int) {
    for (t in seq_along(interactions)) {
      int <- interactions[[t]]
      if (int$type == "macro") {
        col_idx <- match(int$var, mainvars)
        mu_terms <- c(mu_terms, paste0("b.int[", t, "] * F.", int$block1,
                                       "[j] * X.main[j,", col_idx, "]"))
      } else {
        mu_terms <- c(mu_terms, paste0("b.int[", t, "] * F.", int$block1,
                                       "[j] * F.", int$block2, "[j]"))
      }
    }
  }

  if (mm_agg_needed) {
    mu_terms <- c(mu_terms, "mm.agg[j]")
  }

  # HM contributions
  if (has_hm) {
    for (k in seq_along(hm_blocks)) {
      block <- hm_blocks[[k]]
      if (!is.null(block$RE)) {
        re_spec <- block$RE
        if (!is.null(block$ar)) {
          mu_terms <- c(mu_terms, paste0("re.hm.", k, "[hmid[j], n.HMn[j]]"))
        } else {
          if (re_spec$intercept) {
            mu_terms <- c(mu_terms, paste0("re.hm.", k, "[hmid[j]]"))
          }
          for (s in seq_along(re_spec$slopes)) {
            mu_terms <- c(mu_terms, paste0("re.hm.", k, ".s", s, "[hmid[j]] * X.re.hm.", k,
                                           "[j,", s, "]"))
          }
        }
      } else if (!is.null(block$FE)) {
        fe <- block$FE
        if (fe$intercept) {
          mu_terms <- c(mu_terms, paste0("alpha.hm.", k, "[hmid[j]]"))
        }
        for (s in seq_along(fe$slopes)) {
          mu_terms <- c(mu_terms, paste0("alphas.hm.", k, ".s", s, "[hmid[j]] * X.re.hm.", k,
                                         "[j,", s, "]"))
        }
      }
    }
  }

  add("    mu[j] <- ", paste(mu_terms, collapse = " + "))
  add("")

  # Likelihood based on family
  if (family == "Gaussian") {
    add("    Y[j] ~ dnorm(mu[j], tau)")
    if (monitor_has(monitor_spec, "predictive")) {
      add("    pred[j] ~ dnorm(mu[j], tau)")
    }
    if (monitor_has(monitor_spec, "log_lik")) {
      add("    log_lik[j] <- logdensity.norm(Y[j], mu[j], tau)")
    }
  } else if (family == "Binomial") {
    add("    logit(p[j]) <- mu[j]")
    add("    Y[j] ~ dbern(p[j])")
    if (monitor_has(monitor_spec, "predictive")) {
      add("    pred[j] ~ dbern(p[j])")
    }
    if (monitor_has(monitor_spec, "log_lik")) {
      add("    log_lik[j] <- logdensity.bern(Y[j], p[j])")
    }
  } else if (family == "Weibull") {
    add("    lambda[j] <- exp(-mu[j] * shape)")
    add("    t[j] ~ dweib(shape, lambda[j])")
    add("    censored[j] ~ dinterval(t[j], ct.lb[j])")
    if (monitor_has(monitor_spec, "predictive")) {
      add("    pred[j] ~ dweib(shape, lambda[j])")
    }
  } else if (family == "Cox") {
    if (!is.null(cox_intervals) && is.numeric(cox_intervals) && cox_intervals > 0) {
      add("    for (k in 1:n.intervals) {")
      add("      dN_interval[j,k] ~ dpois(Idt[j,k])")
      add("      Idt[j,k] <- Y_interval[j,k] * exp(mu[j]) * lambda0[k]")
      add("    }")
    } else {
      add("    for (k in 1:n.tu) {")
      add("      dN[j,k] ~ dpois(Idt[j,k])")
      add("      Idt[j,k] <- Y[j,k] * exp(mu[j]) * dL0[k]")
      add("    }")
    }
  }

  add("  }")
  add("")

  # ------------------------------------------------------------------------------------------ #
  # Main-level priors
  # ------------------------------------------------------------------------------------------ #

  add("  # Main-level priors")

  n_main_params <- length(mainvars)
  if (n_main_params > 0) {
    for (x in seq_len(n_main_params)) {
      node <- paste0("b[", x, "]")
      add("  ", node, " ~ ", pdist(node, "dnorm(0, 0.0001)"))
    }
  }

  # Interaction coefficients
  if (has_int) {
    for (t in seq_along(interactions)) {
      node <- paste0("b.int[", t, "]")
      add("  ", node, " ~ ", pdist(node, "dnorm(0, 0.0001)"))
    }
  }

  # Family parameters
  if (family == "Gaussian") {
    add(sd_prior_lines("sigma", "tau"))
  } else if (family == "Weibull") {
    add("  shape ~ ", pdist("shape", "dexp(0.001)"))
    add("  tau <- pow(shape, -2)  # approximate")
    add("  sigma <- 1/shape")
  } else if (family == "Cox") {
    if (!is.null(cox_intervals) && is.numeric(cox_intervals) && cox_intervals > 0) {
      add("  for (k in 1:n.intervals) {")
      add("    lambda0[k] ~ dgamma(c, d)")
      add("  }")
    } else {
      add("  for (k in 1:n.tu) {")
      add("    dL0[k] ~ dgamma(c, d)")
      add("  }")
    }
  }
  add("")

  add("}")

  modelstring <- paste(lines, collapse = "\n")

  # ========================================================================================== #
  # Raw-string prior escape hatch (regex replacement, as before)
  # ========================================================================================== #

  if (length(raw_priors) > 0) {
    for (i in seq_along(raw_priors)) {
      full_param <- stringr::str_extract(raw_priors[i], "^[^~<]+") %>% trimws()

      if (stringr::str_detect(full_param, "\\[")) {
        full_param_escaped <- stringr::str_replace_all(full_param, "(\\[|\\]|\\.)", "\\\\\\1")
        pattern <- paste0("(?m)^([ \\t]*)", full_param_escaped, "\\s*(~|<-)\\s*[^\\n]*")
        modelstring <- stringr::str_replace(modelstring, pattern, paste0("\\1", raw_priors[i]))
      } else {
        base_param <- stringr::str_extract(raw_priors[i], "^[^~<]+") %>% trimws()
        operator <- stringr::str_extract(raw_priors[i], "(~|<-)")
        rhs <- stringr::str_extract(raw_priors[i], "(?<=~|<-).*") %>% trimws()
        escaped <- stringr::str_replace_all(base_param, "(\\.|\\[|\\])", "\\\\\\1")
        pattern <- paste0("(?m)^([ \\t]*)(", escaped, "(?:\\[[^\\]]+\\])?)\\s*(~|<-)\\s*[^\\n]*")
        modelstring <- stringr::str_replace_all(modelstring, pattern,
                                                paste0("\\1\\2 ", operator, " ", rhs))
      }
    }
  }

  return(modelstring)
}
