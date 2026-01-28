# ================================================================================================ #
# Function editModelstring
# ================================================================================================ #

editModelstring <- function(family, priors, mm_blocks, main, hm_blocks, mm, hm, DIR, monitor, modelfile) {

  # ========================================================================================== #
  # Flags
  # ========================================================================================== #

  has_mm <- !is.null(mm_blocks) && length(mm_blocks) > 0
  has_hm <- !is.null(hm_blocks) && length(hm_blocks) > 0
  has_mm_RE <- has_mm && attr(mm_blocks, "has_RE")
  has_mm_vars <- has_mm && attr(mm_blocks, "has_vars")

  n.mmblocks <- if (has_mm) length(mm_blocks) else 0
  n.hmblocks <- if (has_hm) length(hm_blocks) else 0

  mainvars <- main$vars
  lhs <- main$lhs

  # ========================================================================================== #
  # Build model string programmatically
  # ========================================================================================== #

  lines <- c()
  add <- function(...) lines <<- c(lines, paste0(...))

  add("model {")
  add("")

  # ------------------------------------------------------------------------------------------ #
  # MM level: mm() blocks - Weight functions and contributions
  # ------------------------------------------------------------------------------------------ #

  if (has_mm) {

    add("  # ==================== MM Level: Multiple Membership ==================== #")
    add("")

    # Weight functions for each mm block
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]
      fn <- block$fn

      add("  # Weight function for mm block ", k)

      # Build the weight function string for JAGS
      fn_string <- fn$string

      # Replace parameters with indexed b.w.k
      if (length(fn$params) > 0) {
        for (p in seq_along(fn$params)) {
          fn_string <- gsub(paste0("\\b", fn$params[p], "\\b"),
                           paste0("b.w.", k, "[", p, "]"), fn_string)
        }
      }

      # Replace variables with indexed X.w.k
      for (v in seq_along(fn$vars)) {
        fn_string <- gsub(paste0("\\b", fn$vars[v], "\\b"),
                         paste0("X.w.", k, "[i,", v, "]"), fn_string)
      }

      add("  for (i in 1:n.mm) {")

      if (fn$constraint) {
        add("    uw.", k, "[i] <- ", fn_string)
        add("    w.", k, "[i] <- uw.", k, "[i] / sum(uw.", k, "[mmi1.mm[i]:mmi2.mm[i]])")
      } else {
        add("    w.", k, "[i] <- ", fn_string)
      }

      # mm contribution for this block (free variables)
      if (!is.null(block$vars) && length(block$vars) > 0) {
        add("    mm.vars.", k, "[i] <- inprod(X.mm.", k, "[i,], b.mm.", k, ")")
      }

      # mm contribution for fixed variables
      if (!is.null(block$vars_fixed) && length(block$vars_fixed) > 0) {
        add("    mm.vars.fix.", k, "[i] <- inprod(X.fix.mm.", k, "[i,], fix.mm.", k, ")")
      }

      add("  }")
      add("")
    }

    # Random effects (shared across blocks that use RE)
    if (has_mm_RE) {
      # Check if any block uses AR
      any_ar <- any(sapply(mm_blocks, function(b) b$ar))

      add("  # MM-level random effects")
      add("  for (i in 1:n.umm) {")
      if (any_ar) {
        add("    re.mm[i,1] ~ dnorm(0, tau.mm)")
        add("    for (t in 2:n.GPNi[i]) {")
        add("      re.mm[i,t] ~ dnorm(re.mm[i,t-1], tau.mm)")
        add("    }")
      } else {
        add("    re.mm[i] ~ dnorm(0, tau.mm)")
      }
      add("  }")
      add("")

      # If AR, extract the appropriate RE value for each mm observation
      if (any_ar) {
        add("  # Extract AR random effects for each mm observation")
        add("  for (i in 1:n.mm) {")
        add("    re.mm.i[i] <- re.mm[mmid[i], n.GPn[i]]")
        add("  }")
        add("")
      }

      add("  # MM-level variance")
      add("  tau.mm ~ dscaled.gamma(25, 1)")
      add("  sigma.mm <- 1/sqrt(tau.mm)")
      add("")
    }

    # Priors for b.mm and b.w for each block
    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]

      if (!is.null(block$vars) && length(block$vars) > 0) {
        add("  # Priors for b.mm.", k)
        add("  for (x in 1:n.Xmm.", k, ") {")
        add("    b.mm.", k, "[x] ~ dnorm(0, 0.0001)")
        add("  }")
        add("")
      }

      if (length(block$fn$params) > 0) {
        add("  # Priors for b.w.", k)
        for (p in seq_along(block$fn$params)) {
          add("  b.w.", k, "[", p, "] ~ dnorm(0, 0.0001)")
        }
        add("")
      }
    }
  }

  # ------------------------------------------------------------------------------------------ #
  # HM level: hm() blocks
  # ------------------------------------------------------------------------------------------ #

  if (has_hm) {

    add("  # ==================== HM Level: Hierarchical Membership ==================== #")
    add("")

    for (k in seq_along(hm_blocks)) {
      block <- hm_blocks[[k]]

      if (block$type == "RE") {
        add("  # HM-level random effects (hm block ", k, ")")

        # Check if this block uses AR
        if (block$ar) {
          # AR: 2D indexing for time-varying effects
          add("  for (k in 1:n.hm) {")
          add("    re.hm.", k, "[k,1] ~ dnorm(0, tau.hm.", k, ")")
          add("    for (t in 2:n.HMNi[k]) {")
          add("      re.hm.", k, "[k,t] ~ dnorm(re.hm.", k, "[k,t-1], tau.hm.", k, ")")
          add("    }")
          add("  }")
          add("")

          # For AR, compute hm.k at main level (to access time-varying RE)
          add("  # HM-level effects computed at main level (AR)")
          add("  for (j in 1:n.main) {")

          # Build hm.k expression at main level
          hm_terms <- c()
          if (!is.null(block$vars) && length(block$vars) > 0) {
            hm_terms <- c(hm_terms, paste0("inprod(X.hm.", k, "[hmid[j],], b.hm.", k, ")"))
          }
          if (!is.null(block$vars_fixed) && length(block$vars_fixed) > 0) {
            hm_terms <- c(hm_terms, paste0("inprod(X.fix.hm.", k, "[hmid[j],], fix.hm.", k, ")"))
          }
          hm_terms <- c(hm_terms, paste0("re.hm.", k, "[hmid[j], n.HMn[j]]"))

          add("    hm.", k, "[j] <- ", paste(hm_terms, collapse = " + "))
          add("  }")
        } else {
          # Non-AR: 1D indexing for constant effects at hm level
          add("  for (k in 1:n.hm) {")

          # Build hm.k expression
          hm_terms <- c()
          if (!is.null(block$vars) && length(block$vars) > 0) {
            hm_terms <- c(hm_terms, paste0("inprod(X.hm.", k, "[k,], b.hm.", k, ")"))
          }
          if (!is.null(block$vars_fixed) && length(block$vars_fixed) > 0) {
            hm_terms <- c(hm_terms, paste0("inprod(X.fix.hm.", k, "[k,], fix.hm.", k, ")"))
          }
          hm_terms <- c(hm_terms, paste0("re.hm.", k, "[k]"))

          add("    hm.", k, "[k] <- ", paste(hm_terms, collapse = " + "))
          add("    re.hm.", k, "[k] ~ dnorm(0, tau.hm.", k, ")")

          add("  }")
        }
        add("")

        add("  # HM-level variance (hm block ", k, ")")
        add("  tau.hm.", k, " ~ dscaled.gamma(25, 1)")
        add("  sigma.hm.", k, " <- 1/sqrt(tau.hm.", k, ")")
        add("")

        if (!is.null(block$vars) && length(block$vars) > 0) {
          add("  # Priors for b.hm.", k)
          add("  for (x in 1:n.Xhm.", k, ") {")
          add("    b.hm.", k, "[x] ~ dnorm(0, 0.0001)")
          add("  }")
          add("")
        }

      } else {
        # Fixed effects
        add("  # HM-level fixed effects (hm block ", k, ")")
        add("  for (k in 1:n.hm) {")

        # Build hm.k expression
        hm_terms <- c()
        if (!is.null(block$vars) && length(block$vars) > 0) {
          hm_terms <- c(hm_terms, paste0("inprod(X.hm.", k, "[k,], b.hm.", k, ")"))
        }
        if (!is.null(block$vars_fixed) && length(block$vars_fixed) > 0) {
          hm_terms <- c(hm_terms, paste0("inprod(X.fix.hm.", k, "[k,], fix.hm.", k, ")"))
        }

        add("    hm.", k, "[k] <- ", paste(hm_terms, collapse = " + "))
        add("  }")
        add("")

        if (!is.null(block$vars) && length(block$vars) > 0) {
          add("  # Priors for b.hm.", k, " (fixed effects)")
          add("  for (x in 1:n.Xhm.", k, ") {")
          add("    b.hm.", k, "[x] ~ dnorm(0, 0.0001)")
          add("  }")
          add("")
        }
      }
    }
  }

  # ------------------------------------------------------------------------------------------ #
  # Main level: Main model
  # ------------------------------------------------------------------------------------------ #

  add("  # ==================== Main Level: Main Model ==================== #")
  add("")

  # Build the linear predictor
  add("  for (j in 1:n.main) {")
  add("")

  # Aggregate mm contributions
  if (has_mm) {
    add("    # Aggregate mm-level contributions")
    add("    mm.agg[j] <- sum(")

    mm_terms <- c()
    # Check once if any mm block uses AR (determines re.mm dimensionality)
    any_ar <- any(sapply(mm_blocks, function(b) b$ar))

    for (k in seq_along(mm_blocks)) {
      block <- mm_blocks[[k]]

      terms_k <- c()

      # Variables contribution (free)
      if (!is.null(block$vars) && length(block$vars) > 0) {
        terms_k <- c(terms_k, paste0("w.", k, "[mmi1[j]:mmi2[j]] * mm.vars.", k, "[mmi1[j]:mmi2[j]]"))
      }

      # Variables contribution (fixed)
      if (!is.null(block$vars_fixed) && length(block$vars_fixed) > 0) {
        terms_k <- c(terms_k, paste0("w.", k, "[mmi1[j]:mmi2[j]] * mm.vars.fix.", k, "[mmi1[j]:mmi2[j]]"))
      }

      # RE contribution
      # Note: All mm blocks share re.mm, so indexing must be consistent across blocks
      if (block$RE) {
        if (any_ar) {
          # Use pre-extracted re.mm.i which has AR values at mm-level
          terms_k <- c(terms_k, paste0("w.", k, "[mmi1[j]:mmi2[j]] * re.mm.i[mmi1[j]:mmi2[j]]"))
        } else {
          terms_k <- c(terms_k, paste0("w.", k, "[mmi1[j]:mmi2[j]] * re.mm[mmid[mmi1[j]:mmi2[j]]]"))
        }
      }

      if (length(terms_k) > 0) {
        mm_terms <- c(mm_terms, paste0("      ", paste(terms_k, collapse = " + ")))
      }
    }

    add(paste(mm_terms, collapse = " +\n"), ")")
    add("")
  }

  # Build mu[j]
  mu_terms <- c()

  # Intercept (estimated)
  if ("X0" %in% mainvars) {
    mu_terms <- c(mu_terms, "b[1]")
  }

  # Main-level covariates (excluding intercept)
  mainvars_no_intercept <- mainvars[mainvars != "X0"]
  if (length(mainvars_no_intercept) > 0) {
    if ("X0" %in% mainvars) {
      # Intercept is b[1], other vars start at 2
      mu_terms <- c(mu_terms, paste0("inprod(X.main[j,], b[2:", length(mainvars), "])"))
    } else {
      mu_terms <- c(mu_terms, "inprod(X.main[j,], b)")
    }
  }

  # Main-level fixed covariates (including fixed intercept if present)
  if (!is.null(main$vars_fixed) && length(main$vars_fixed) > 0) {
    mu_terms <- c(mu_terms, "inprod(X.fix.main[j,], fix.main)")
  }

  # MM-level aggregated
  if (has_mm) {
    mu_terms <- c(mu_terms, "mm.agg[j]")
  }

  # HM-level contributions
  if (has_hm) {
    for (k in seq_along(hm_blocks)) {
      # For AR blocks, hm.k is computed at main level (indexed by j)
      # For non-AR blocks, hm.k is computed at hm level (indexed by hmid[j])
      if (hm_blocks[[k]]$ar && hm_blocks[[k]]$type == "RE") {
        mu_terms <- c(mu_terms, paste0("hm.", k, "[j]"))
      } else {
        mu_terms <- c(mu_terms, paste0("hm.", k, "[hmid[j]]"))
      }
    }
  }

  add("    mu[j] <- ", paste(mu_terms, collapse = " + "))
  add("")

  # Likelihood based on family
  if (family == "Gaussian") {
    add("    Y[j] ~ dnorm(mu[j], tau)")
    if (monitor) {
      add("    pred[j] ~ dnorm(mu[j], tau)")
    }
  } else if (family == "Binomial") {
    add("    logit(p[j]) <- mu[j]")
    add("    Y[j] ~ dbern(p[j])")
    if (monitor) {
      add("    pred[j] ~ dbern(p[j])")
    }
  } else if (family == "Weibull") {
    add("    lambda[j] <- exp(-mu[j] * shape)")
    add("    t[j] ~ dweib(shape, lambda[j])")
    add("    censored[j] ~ dinterval(t[j], ct.lb[j])")
    if (monitor) {
      add("    pred[j] ~ dweib(shape, lambda[j])")
    }
  } else if (family == "Cox") {
    add("    for (k in 1:n.tu) {")
    add("      dN[j,k] ~ dpois(Idt[j,k])")
    add("      Idt[j,k] <- Y[j,k] * exp(mu[j]) * dL0[k]")
    add("    }")
  }

  add("  }")
  add("")

  # ------------------------------------------------------------------------------------------ #
  # Main-level priors
  # ------------------------------------------------------------------------------------------ #

  add("  # Main-level priors")

  # Intercept and coefficients
  n_main_params <- length(mainvars)
  if (n_main_params > 0) {
    add("  for (x in 1:", n_main_params, ") {")
    add("    b[x] ~ dnorm(0, 0.0001)")
    add("  }")
  }

  # Variance
  if (family == "Gaussian") {
    add("  tau ~ dscaled.gamma(25, 1)")
    add("  sigma <- 1/sqrt(tau)")
  } else if (family == "Weibull") {
    add("  shape ~ dexp(0.001)")
    add("  tau <- pow(shape, -2)  # approximate")
    add("  sigma <- 1/shape")
  } else if (family == "Cox") {
    add("  for (k in 1:n.tu) {")
    add("    dL0[k] ~ dgamma(c, d)")
    add("  }")
  }
  add("")

  add("}")

  # Join lines
  modelstring <- paste(lines, collapse = "\n")

  # ========================================================================================== #
  # Apply user-specified priors
  # ========================================================================================== #

  if (!is.null(priors)) {
    for (i in seq_along(priors)) {
      # Extract full param (e.g., b.w.1 or b.w.1[1])
      full_param <- stringr::str_extract(priors[i], "^[^~<]+") %>% trimws()

      # Check if full_param is indexed (e.g., b.w.1[1])
      if (stringr::str_detect(full_param, "\\[")) {
        # Replace prior for one specific parameter
        full_param_escaped <- stringr::str_replace_all(full_param, "(\\[|\\]|\\.)", "\\\\\\1")
        pattern <- paste0("(?m)^([ \\t]*)", full_param_escaped, "\\s*(~|<-)\\s*[^\\n]*")
        modelstring <- stringr::str_replace(modelstring, pattern, paste0("\\1", priors[i]))
      } else {
        # Replace prior for all indices of a parameter
        base_param <- stringr::str_extract(priors[i], "^[^~<]+") %>% trimws()
        operator <- stringr::str_extract(priors[i], "(~|<-)")
        rhs <- stringr::str_extract(priors[i], "(?<=~|<-).*") %>% trimws()
        escaped <- stringr::str_replace_all(base_param, "(\\.|\\[|\\])", "\\\\\\1")
        pattern <- paste0("(?m)^([ \\t]*)(", escaped, "(?:\\[[^\\]]+\\])?)\\s*(~|<-)\\s*[^\\n]*")
        modelstring <- stringr::str_replace_all(modelstring, pattern, paste0("\\1\\2 ", operator, " ", rhs))
      }
    }
  }

  # ========================================================================================== #
  # Save or return
  # ========================================================================================== #

  return(modelstring)

}
