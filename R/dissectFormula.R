# ================================================================================================ #
# Function dissectFormula
# ================================================================================================ #

dissectFormula <- function(formula, family, data) {

  # This function takes the formula object and turns it into structures for further processing.
  # Arguments:
  # - formula (formula) formula object
  # - family (str) internal family string ("Gaussian", "Binomial", "Weibull", "Cox")
  # - data (df) data.frame to resolve variables, parameters, and block names against
  # Return:
  # - list(lhs, mainvars, mainvars_fixed, main_formula, mm, mm_groups, hm, interactions)

  # --------------------------------------------------------------------------------------------- #
  # Extract left- and right-hand side from formula
  # --------------------------------------------------------------------------------------------- #

  lhs <- formula[[2]] %>% all.vars()
  rhs_terms <- attr(terms(formula), "term.labels")

  # Validate LHS based on family
  if (family %in% c("Gaussian", "Binomial") && length(lhs) > 1) {
    stop("family = ", tolower(family), " takes only one variable on the left-hand side.",
         call. = FALSE)
  }
  if (family %in% c("Weibull", "Cox") && length(lhs) != 2) {
    stop("family = ", tolower(family),
         " takes two variables on the left-hand side: 'Surv(survtime, event)'", call. = FALSE)
  }

  # --------------------------------------------------------------------------------------------- #
  # Separate mm(), hm(), fix(), and regular terms
  # --------------------------------------------------------------------------------------------- #

  mm_terms <- rhs_terms[startsWith(rhs_terms, "mm(")]
  hm_terms <- rhs_terms[startsWith(rhs_terms, "hm(")]
  fix_terms <- rhs_terms[startsWith(rhs_terms, "fix(")]
  regular_terms <- rhs_terms[!startsWith(rhs_terms, "mm(") &
                               !startsWith(rhs_terms, "hm(") &
                               !startsWith(rhs_terms, "fix(")]

  # Parse fix() terms at main level
  eval_env_fix <- list2env(list(fix = fix), parent = baseenv())
  mainvars_fixed <- list()

  if (length(fix_terms) > 0) {
    for (term in fix_terms) {
      parsed <- tryCatch(
        eval(parse(text = term), envir = eval_env_fix),
        error = function(e) {
          stop("Error parsing fix() term: ", term, "\n", e$message, call. = FALSE)
        }
      )
      if (inherits(parsed, "bml_fix")) {
        mainvars_fixed[[length(mainvars_fixed) + 1]] <- parsed
      }
    }
  }

  # --------------------------------------------------------------------------------------------- #
  # Parse mm() blocks
  # --------------------------------------------------------------------------------------------- #

  eval_env <- list2env(list(
    id = id, vars = vars, w = w, fn = fn, mm = mm, re = re, fe = fe, est = est
  ), parent = baseenv())

  mm_list <- list()
  mmid_groups <- NULL
  param_msgs <- character(0)

  if (length(mm_terms) > 0) {

    mm_list <- lapply(mm_terms, function(term) {
      tryCatch(
        eval(parse(text = term), envir = eval_env),
        error = function(e) {
          stop("Error parsing mm() block: ", term, "\n", conditionMessage(e), call. = FALSE)
        }
      )
    })

    # Validate: all mm() blocks must share the same mainid
    first_mainid <- mm_list[[1]]$id[2]
    for (i in seq_along(mm_list)) {
      if (mm_list[[i]]$id[2] != first_mainid) {
        stop("All mm() blocks must share the same mainid (second element of id()). ",
             "Found: mainid='", first_mainid, "' and mainid='", mm_list[[i]]$id[2], "'",
             call. = FALSE)
      }
    }

    # Group blocks by mmid
    mmid_groups <- list()
    for (i in seq_along(mm_list)) {
      mmid_name <- mm_list[[i]]$id[1]
      mmid_groups[[mmid_name]] <- c(mmid_groups[[mmid_name]], i)
    }
    all_mmid_names <- names(mmid_groups)

    # Validate: id variables exist in data
    all_id_vars <- unique(c(sapply(mm_list, function(m) m$id[1]), first_mainid))
    if (!all(all_id_vars %in% names(data))) {
      missing_ids <- all_id_vars[!all_id_vars %in% names(data)]
      stop("ID variable(s) not found in data: ", paste(missing_ids, collapse = ", "),
           call. = FALSE)
    }

    # ------------------------------------------------------------------------------------------- #
    # Per-block resolution: vars, weights (one parameter rule), fn, labels, names
    # ------------------------------------------------------------------------------------------- #

    for (k in seq_along(mm_list)) {
      block <- mm_list[[k]]

      # ----- vars: validate against data, compute design column names --------------------------- #
      attr_cols <- character(0)
      if (!is.null(block$vars)) {
        base_vars <- c(
          if (!is.null(block$vars$formula)) all.vars(block$vars$formula) else block$vars$free,
          if (!is.null(block$vars$fixed)) sapply(block$vars$fixed, function(x) x$var)
        )
        missing_vars <- setdiff(base_vars, names(data))
        if (length(missing_vars) > 0) {
          stop("mm() variable(s) not found in data: ", paste(missing_vars, collapse = ", "),
               call. = FALSE)
        }
        if (!is.null(block$vars$formula)) {
          mmdat_k <- data[!is.na(data[[block$id[1]]]), , drop = FALSE]
          attr_cols <- colnames(model.matrix(block$vars$formula, data = mmdat_k))
        }
      }
      block$attr_cols <- attr_cols

      # ----- fn: arity, attribute resolution, DSL graph, shape parameters ----------------------- #
      fnobj <- block$fn
      arity <- if (fnobj$type %in% names(FN_ARITY)) FN_ARITY[[fnobj$type]] else NA
      if (!is.na(arity) && fnobj$type != "expr") {
        if (arity == 0 && length(attr_cols) > 0) {
          stop("fn(\"", fnobj$type, "\") takes no attributes; remove vars().", call. = FALSE)
        }
        if (arity > 0 && length(attr_cols) != arity) {
          stop("fn(\"", fnobj$type, "\") requires exactly ", arity, " member attribute",
               if (arity > 1) "s", " (found ", length(attr_cols), " design column",
               if (length(attr_cols) != 1) "s", " in vars()).", call. = FALSE)
        }
        if (!is.null(block$vars$fixed)) {
          stop("fix() variables are not supported with fn(\"", fnobj$type, "\").", call. = FALSE)
        }
      }

      dsl_params <- character(0)
      if (fnobj$type == "expr") {
        syms <- dsl_symbols(fnobj$expr)
        block_attr_names <- if (!is.null(block$vars)) block$vars$free else character(0)
        fn_attrs <- intersect(syms, block_attr_names)
        dsl_params <- setdiff(syms, c(block_attr_names, "w"))
        if (length(fn_attrs) == 0 && !"w" %in% syms) {
          stop("fn() expression references no member-level quantity: use the vars() ",
               "attribute names or the reserved symbol w inside E().", call. = FALSE)
        }
        graph <- fn_expr_graph(fnobj$expr)
        validate_dsl_resolved(graph, attrs = fn_attrs, params = dsl_params)
        fnobj$graph <- graph
        fnobj$attrs <- fn_attrs
        # design columns actually used by the expression, in vars() order
        attr_cols_expr <- intersect(block_attr_names, fn_attrs)
        block$attr_cols <- attr_cols_expr
        attr_cols <- attr_cols_expr
      }

      # x range for the threshold cutpoint default prior
      x_range <- NULL
      if (fnobj$type %in% c("threshold", "smax", "gmean", "var") && length(attr_cols) >= 1 &&
          attr_cols[1] %in% names(data)) {
        xv <- data[[attr_cols[1]]]
        if (is.numeric(xv)) x_range <- range(xv, na.rm = TRUE)
      }

      fnobj$est_params <- fn_est_params(fnobj, dsl_params = dsl_params, x_range = x_range)
      fnobj$fixed_params <- Filter(function(p) !is_est(p), fnobj$params %||% list())
      block$fn <- fnobj

      # ----- weights: one parameter rule ---------------------------------------------------------- #
      wobj <- block$w
      agg_cols <- if (!is.null(wobj$agg_funcs)) sapply(wobj$agg_funcs, function(a) a$col_name) else character(0)
      reserved <- c("n", "X0", agg_cols)
      wsyms <- wobj$symbols

      # base variables inside aggregation shortcuts must exist in the data
      if (!is.null(wobj$agg_vars)) {
        missing_agg <- setdiff(wobj$agg_vars, names(data))
        if (length(missing_agg) > 0) {
          stop("Variable(s) used in weight aggregation functions not found in data: ",
               paste(missing_agg, collapse = ", "), call. = FALSE)
        }
      }

      wobj$params <- setdiff(wsyms, c(names(data), reserved))
      wobj$vars <- setdiff(wsyms, wobj$params)

      # variables that have parameters (param * var pattern) - for labeling
      if (length(wobj$params) > 0) {
        vp <- stringr::str_match_all(
          wobj$string,
          paste0("\\b(", paste(wobj$params, collapse = "|"), ")\\s*\\*\\s*([a-zA-Z_][a-zA-Z0-9_\\.]*)")
        )[[1]]
        wobj$vars_p <- if (nrow(vp) > 0) vp[, 3] else character(0)
      } else {
        wobj$vars_p <- character(0)
      }
      block$w <- wobj

      # tilted-mean near-equivalence: same variable in weights and vars()
      overlap <- intersect(wobj$vars, if (!is.null(block$vars)) block$vars$free else character(0))
      overlap <- c(overlap, intersect(wobj$agg_vars %||% character(0),
                                      if (!is.null(block$vars)) block$vars$free else character(0)))
      if (length(overlap) > 0 && length(wobj$params) > 0) {
        warning("mm() block ", k, ": variable(s) ", paste(unique(overlap), collapse = ", "),
                " appear in both the weight formula and vars(). Attribute-dependent weights ",
                "are near-equivalent to a smooth-max in fn() (tilted mean); prefer ",
                "fn(\"smax\", kappa = est()) so weights keep their measurement role.",
                call. = FALSE)
      }

      # ----- RE slopes: member-level columns ------------------------------------------------------ #
      if (!is.null(block$RE) && length(block$RE$slopes) > 0) {
        missing_sl <- setdiff(block$RE$slopes, names(data))
        if (length(missing_sl) > 0) {
          stop("re() slope variable(s) not found in data: ", paste(missing_sl, collapse = ", "),
               call. = FALSE)
        }
        not_in_vars <- setdiff(block$RE$slopes,
                               if (!is.null(block$vars)) block$vars$free else character(0))
        if (length(not_in_vars) > 0) {
          warning("mm() block ", k, ": re() slope variable(s) ",
                  paste(not_in_vars, collapse = ", "),
                  " have no mean effect in this block's vars(). A random slope is the residual ",
                  "half of the effect heterogeneity; usually the mean effect belongs in vars().",
                  call. = FALSE)
        }
      }
      if (!is.null(block$FE) && length(block$FE$slopes) > 0) {
        missing_sl <- setdiff(block$FE$slopes, names(data))
        if (length(missing_sl) > 0) {
          stop("fe() slope variable(s) not found in data: ", paste(missing_sl, collapse = ", "),
               call. = FALSE)
        }
      }

      # ----- time-indexed AR walk: validate the time column ---------------------------------- #
      if (!is.null(block$RE) && !is.null(block$RE$ar) && !is.null(block$RE$ar$time)) {
        tv <- block$RE$ar$time
        if (!tv %in% names(data)) {
          stop("re(ar = ", tv, "): time variable '", tv, "' not found in data.", call. = FALSE)
        }
        tvals <- data[[tv]][!is.na(data[[block$id[1]]])]
        if (!is.numeric(tvals)) {
          stop("re(ar = ", tv, "): the time variable must be numeric ",
               "(convert Dates, e.g. as.numeric(format(date, \"%Y\"))).", call. = FALSE)
        }
        if (any(is.na(tvals))) {
          stop("re(ar = ", tv, "): the time variable contains missing values.", call. = FALSE)
        }
        dup <- data[!is.na(data[[block$id[1]]]), c(block$id[1], tv)]
        if (anyDuplicated(dup)) {
          stop("re(ar = ", tv, "): duplicated time values within a member ",
               "(a zero time gap makes the random-walk step variance zero). ",
               "Each member's participations must have unique times.", call. = FALSE)
        }
      }

      # ----- feature labels and block name -------------------------------------------------------- #
      block$feature_labels <- fn_feature_labels(fnobj, attr_cols)
      if (fnobj$type == "expr" && !is.null(block$name)) {
        block$feature_labels <- block$name
      } else if (fnobj$type == "expr") {
        block$feature_labels <- paste0("f_", k)
      }
      if (is.null(block$name)) {
        block$name <- if (length(block$feature_labels) == 1) block$feature_labels else paste0("mm", k)
      } else if (length(block$feature_labels) == 1) {
        block$feature_labels <- block$name
      }
      if (block$name %in% names(data)) {
        stop("mm() block name '", block$name, "' collides with a data column. ",
             "Choose a different name.", call. = FALSE)
      }

      # ----- data-driven identification warnings --------------------------------------------------- #
      if (fnobj$type %in% c("var", "smax", "threshold", "gmean", "cov") && length(attr_cols) >= 1) {
        xcol <- attr_cols[1]
        if (xcol %in% names(data)) {
          mainid_var <- block$id[2]
          wg_var <- tapply(data[[xcol]], data[[mainid_var]], stats::var, na.rm = TRUE)
          wg_var <- wg_var[!is.na(wg_var)]
          if (length(wg_var) > 0 && all(wg_var < 1e-12)) {
            warning("mm() block ", k, ": '", xcol, "' is constant within groups, so the ",
                    "fn(\"", fnobj$type, "\") feature has no variation and its coefficient ",
                    "is not identified.", call. = FALSE)
          }
        }
      }
      if (fnobj$type == "gmean" && length(attr_cols) >= 1 && attr_cols[1] %in% names(data)) {
        if (any(data[[attr_cols[1]]] <= 0, na.rm = TRUE)) {
          stop("fn(\"gmean\") requires strictly positive attribute values ",
               "(the power mean is undefined otherwise); '", attr_cols[1],
               "' has values <= 0.", call. = FALSE)
        }
      }

      # collect the "Estimating parameters" message parts
      if (length(wobj$params) > 0) {
        param_msgs <- c(param_msgs,
                        paste0(paste(wobj$params, collapse = ", "), " (w.", k, ")"))
      }
      if (length(fnobj$est_params) > 0) {
        param_msgs <- c(param_msgs,
                        paste0(paste(names(fnobj$est_params), collapse = ", "), " (fn.", k, ")"))
      }

      block$mmid_group <- which(all_mmid_names == block$id[1])
      mm_list[[k]] <- block
    }

    # ------------------------------------------------------------------------------------------- #
    # Cross-block validations
    # ------------------------------------------------------------------------------------------- #

    # Duplicate block names
    block_names <- sapply(mm_list, function(m) m$name)
    if (anyDuplicated(block_names)) {
      dups <- unique(block_names[duplicated(block_names)])
      stop("Duplicate mm() block name(s): ", paste(dups, collapse = ", "),
           ". Use name= to disambiguate.", call. = FALSE)
    }

    # At least one block must contribute something
    has_effect <- any(sapply(mm_list, function(m) {
      !is.null(m$RE) || !is.null(m$FE) || length(m$feature_labels) > 0
    }))
    if (!has_effect) {
      stop("At least one mm() block must have RE/FE or emit a feature (vars/fn).", call. = FALSE)
    }

    # RE = one per mmid group
    for (mmid_name in names(mmid_groups)) {
      block_indices <- mmid_groups[[mmid_name]]
      n_re_in_group <- sum(sapply(block_indices, function(i) !is.null(mm_list[[i]]$RE)))
      if (n_re_in_group > 1) {
        stop("Random effects can only be specified for one mm() block per mmid. ",
             "Found ", n_re_in_group, " blocks with RE sharing mmid='", mmid_name, "'",
             call. = FALSE)
      }
    }
  }

  # --------------------------------------------------------------------------------------------- #
  # Parse hm() blocks
  # --------------------------------------------------------------------------------------------- #

  eval_env_hm <- list2env(list(id = id, hm = hm, re = re, fe = fe), parent = baseenv())

  hm_list <- list()

  if (length(hm_terms) > 0) {

    hm_list <- lapply(hm_terms, function(term) {
      tryCatch(
        eval(parse(text = term), envir = eval_env_hm),
        error = function(e) {
          stop("Error parsing hm() block: ", term, "\n", conditionMessage(e), call. = FALSE)
        }
      )
    })

    for (i in seq_along(hm_list)) {
      h <- hm_list[[i]]
      if (!h$id %in% names(data)) {
        stop("hm() id variable not found in data: ", h$id, call. = FALSE)
      }
      if (!is.null(h$labels) && !h$labels %in% names(data)) {
        stop("hm() labels variable not found in data: ", h$labels, call. = FALSE)
      }
      # slope variables must exist; they are main-level columns
      slopes <- c(if (!is.null(h$RE)) h$RE$slopes, if (!is.null(h$FE)) h$FE$slopes)
      missing_sl <- setdiff(slopes, names(data))
      if (length(missing_sl) > 0) {
        stop("hm() slope variable(s) not found in data: ", paste(missing_sl, collapse = ", "),
             call. = FALSE)
      }

      # time-indexed AR walk: validate the time column (uniqueness within unit is
      # checked at the main level, i.e. one time value per main unit observation)
      if (!is.null(h$RE) && !is.null(h$RE$ar) && !is.null(h$RE$ar$time)) {
        tv <- h$RE$ar$time
        if (!tv %in% names(data)) {
          stop("re(ar = ", tv, "): time variable '", tv, "' not found in data.", call. = FALSE)
        }
        if (!is.numeric(data[[tv]])) {
          stop("re(ar = ", tv, "): the time variable must be numeric ",
               "(convert Dates, e.g. as.numeric(format(date, \"%Y\"))).", call. = FALSE)
        }
        if (any(is.na(data[[tv]]))) {
          stop("re(ar = ", tv, "): the time variable contains missing values.", call. = FALSE)
        }
      }
    }
  }

  # --------------------------------------------------------------------------------------------- #
  # Named-block references in the main formula: interactions
  # --------------------------------------------------------------------------------------------- #

  block_names <- if (length(mm_list) > 0) sapply(mm_list, function(m) m$name) else character(0)
  interactions <- list()
  keep <- rep(TRUE, length(regular_terms))

  for (t in seq_along(regular_terms)) {
    term <- regular_terms[t]
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]

    if (length(parts) == 1) {
      # bare block-name reference is an error (the block's own term already enters mu)
      if (term %in% block_names) {
        stop("Block '", term, "' is referenced bare in the main formula. Its main effect ",
             "already enters the model through the mm() block itself - reference it only in ",
             "interactions (e.g. ", term, ":x).", call. = FALSE)
      }
      next
    }

    if (length(parts) == 2 && any(parts %in% block_names)) {
      keep[t] <- FALSE
      b1 <- which(parts %in% block_names)[1]
      other <- parts[-b1]
      k1 <- which(block_names == parts[b1])

      if (length(mm_list[[k1]]$feature_labels) != 1) {
        stop("Block '", parts[b1], "' emits ", length(mm_list[[k1]]$feature_labels),
             " features; only single-feature blocks can be referenced in interactions.",
             call. = FALSE)
      }

      if (other %in% block_names) {
        # block x block
        k2 <- which(block_names == other)
        if (length(mm_list[[k2]]$feature_labels) != 1) {
          stop("Block '", other, "' emits multiple features; only single-feature blocks can ",
               "be referenced in interactions.", call. = FALSE)
        }
        interactions[[length(interactions) + 1]] <- list(
          type = "block",
          block1 = k1, block2 = k2,
          label = paste0(mm_list[[k1]]$feature_labels, ":", mm_list[[k2]]$feature_labels)
        )
      } else {
        # cross-level: block feature x macro variable
        if (!other %in% names(data)) {
          stop("Interaction term '", term, "': '", other, "' is neither a data column nor a ",
               "block name.", call. = FALSE)
        }
        # the macro variable must also appear as a main effect
        if (!other %in% regular_terms) {
          stop("Cross-level interaction '", term, "' requires '", other,
               "' as a main effect too (add + ", other, " to the formula).", call. = FALSE)
        }
        # validate the interacting covariate is macro-level (constant within mainid)
        mainid_var <- mm_list[[k1]]$id[2]
        wg <- tapply(data[[other]], data[[mainid_var]], function(v) length(unique(v)))
        if (any(wg > 1, na.rm = TRUE)) {
          warning("Cross-level interaction '", term, "': '", other, "' varies within ",
                  mainid_var, " groups. The cross-level collapse sum_k w_k x_k z = z A_x only ",
                  "holds for macro-level (group-constant) variables; a within-group-varying ",
                  "moderator is a member-paired interaction and belongs inside vars() ",
                  "(e.g. vars(x:", other, ")).", call. = FALSE)
        }
        interactions[[length(interactions) + 1]] <- list(
          type = "macro",
          block1 = k1, var = other,
          label = paste0(mm_list[[k1]]$feature_labels, ":", other)
        )
      }
    }
  }
  regular_terms <- regular_terms[keep]

  # --------------------------------------------------------------------------------------------- #
  # Main-level formula
  # --------------------------------------------------------------------------------------------- #

  mainvars <- c(
    if (attr(terms(formula), "intercept") == 1) "X0",
    regular_terms
  )

  has_intercept <- attr(terms(formula), "intercept") == 1
  if (length(regular_terms) > 0) {
    main_formula_str <- paste(regular_terms, collapse = " + ")
    if (has_intercept) {
      main_formula <- as.formula(paste("~", main_formula_str))
    } else {
      main_formula <- as.formula(paste("~ 0 +", main_formula_str))
    }
  } else {
    if (has_intercept) {
      main_formula <- as.formula("~ 1")
    } else {
      main_formula <- NULL
    }
  }

  # Validate mainvars exist in data
  if (!is.null(main_formula)) {
    mainvars_base <- all.vars(main_formula)
  } else {
    mainvars_base <- c()
  }
  if (length(mainvars_base) > 0 && !all(mainvars_base %in% names(data))) {
    missing_vars <- mainvars_base[!mainvars_base %in% names(data)]
    stop("Main-level variable(s) not found in data: ", paste(missing_vars, collapse = ", "),
         call. = FALSE)
  }

  # Validate fixed main vars exist in data
  if (length(mainvars_fixed) > 0) {
    fixed_var_names <- sapply(mainvars_fixed, function(x) x$var)
    fixed_var_names_check <- fixed_var_names[fixed_var_names != "X0"]
    if (length(fixed_var_names_check) > 0 && !all(fixed_var_names_check %in% names(data))) {
      missing_vars <- fixed_var_names_check[!fixed_var_names_check %in% names(data)]
      stop("Fixed main-level variable(s) not found in data: ",
           paste(missing_vars, collapse = ", "), call. = FALSE)
    }
  }

  # hm slope variables should normally have a main-formula fixed effect (lme4 convention)
  for (h in hm_list) {
    slopes <- c(if (!is.null(h$RE)) h$RE$slopes, if (!is.null(h$FE)) h$FE$slopes)
    no_main <- setdiff(slopes, regular_terms)
    if (length(no_main) > 0) {
      warning("hm(", h$id, "): slope variable(s) ", paste(no_main, collapse = ", "),
              " have no main-formula fixed effect. This is legal (as in lme4) but usually ",
              "a mistake - add the variable to the main formula.", call. = FALSE)
    }
  }

  # --------------------------------------------------------------------------------------------- #
  # One message listing every estimated parameter (the typo fence)
  # --------------------------------------------------------------------------------------------- #

  if (length(param_msgs) > 0) {
    message("Estimating parameters: ", paste(param_msgs, collapse = "; "))
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
      hm             = hm_list,
      interactions   = interactions
    )
  )
}
