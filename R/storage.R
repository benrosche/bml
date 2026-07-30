# ================================================================================================ #
# Posterior storage specifications and compact draw helpers
# ================================================================================================ #

.bml_monitor_capabilities <- c(
  "parameters", "random_effects", "weights", "fitted", "predictive", "log_lik"
)

normalize_monitor <- function(monitor, family) {
  if (is.logical(monitor)) {
    stop(
      "`monitor` is now capability-based. Use monitor = \"summary\" for a compact fit, ",
      "monitor = \"full\" for all posterior draws, or monitor = \"jags\" for the raw ",
      "R2jags object.",
      call. = FALSE
    )
  }
  if (!is.character(monitor) || length(monitor) == 0L || anyNA(monitor)) {
    stop("`monitor` must be a non-empty character vector of storage capabilities.",
         call. = FALSE)
  }

  monitor <- unique(monitor)
  allowed <- c("summary", .bml_monitor_capabilities, "full", "jags")
  unknown <- setdiff(monitor, allowed)
  if (length(unknown) > 0L) {
    stop(
      "Unknown `monitor` capability/capabilities: ", paste(unknown, collapse = ", "),
      ". Available values are: ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if ("jags" %in% monitor && length(setdiff(monitor, c("summary", "jags"))) > 0L) {
    stop("`monitor = \"jags\"` is exclusive; do not combine it with other capabilities.",
         call. = FALSE)
  }
  if ("full" %in% monitor && length(setdiff(monitor, c("summary", "full"))) > 0L) {
    stop("`monitor = \"full\"` is a preset; do not combine it with other capabilities.",
         call. = FALSE)
  }

  family_caps <- .bml_monitor_capabilities
  if (family == "Cox") {
    family_caps <- setdiff(family_caps, c("predictive", "log_lik"))
  } else if (family == "Weibull") {
    family_caps <- setdiff(family_caps, "log_lik")
  }

  if ("log_lik" %in% monitor && !family %in% c("Gaussian", "Binomial")) {
    stop("`monitor = \"log_lik\"` is currently available only for gaussian and ",
         "bernoulli families.", call. = FALSE)
  }
  if ("predictive" %in% monitor && family == "Cox") {
    stop("`monitor = \"predictive\"` is not currently available for Cox models.",
         call. = FALSE)
  }

  raw_jags <- "jags" %in% monitor
  preset_full <- "full" %in% monitor
  capabilities <- if (raw_jags || preset_full) {
    family_caps
  } else {
    intersect(.bml_monitor_capabilities, monitor)
  }

  list(
    requested = monitor,
    capabilities = capabilities,
    raw_jags = raw_jags,
    preset = if (raw_jags) "jags" else if (preset_full) "full" else NULL
  )
}

monitor_has <- function(spec, capability) {
  capability %in% spec$capabilities
}

bml_storage_label <- function(spec) {
  if (isTRUE(spec$raw_jags)) return("jags")
  if (identical(spec$preset, "full")) return("full")
  caps <- spec$capabilities
  if (length(caps) == 0L) "summary" else paste(caps, collapse = " + ")
}

bml_raw_sims_array <- function(jags.out) {
  jags.out$BUGSoutput$sims.array
}

bml_select_draws <- function(jags.out, variables) {
  arr <- bml_raw_sims_array(jags.out)
  variables <- intersect(unique(variables), dimnames(arr)[[3]])
  if (length(variables) == 0L) return(NULL)
  posterior::as_draws_array(arr[, , variables, drop = FALSE])
}

bml_indexed_nodes <- function(variables, prefix) {
  variables[grepl(paste0("^", prefix, "(\\[|$)"), variables)]
}

bml_reported_draw_variables <- function(reg.table, variables) {
  intersect(rownames(reg.table), variables)
}

bml_retained_draw_variables <- function(jags.out, reg.table, monitor_spec) {
  variables <- dimnames(bml_raw_sims_array(jags.out))[[3]]
  keep <- character()

  if (monitor_has(monitor_spec, "parameters")) {
    keep <- c(keep, bml_reported_draw_variables(reg.table, variables))
  }
  if (monitor_has(monitor_spec, "random_effects")) {
    keep <- c(
      keep,
      variables[grepl("^re\\.(mm|hm)\\.\\d+(\\.s\\d+)?\\[", variables)]
    )
  }
  if (monitor_has(monitor_spec, "weights")) {
    keep <- c(keep, variables[grepl("^w\\.\\d+\\[", variables)])
  }
  if (monitor_has(monitor_spec, "fitted")) {
    keep <- c(keep, bml_indexed_nodes(variables, "mu"))
  }
  if (monitor_has(monitor_spec, "predictive")) {
    keep <- c(keep, bml_indexed_nodes(variables, "pred"))
  }
  if (monitor_has(monitor_spec, "log_lik")) {
    keep <- c(keep, bml_indexed_nodes(variables, "log_lik"))
  }

  unique(keep)
}

bml_parameter_diagnostics <- function(jags.out, reg.table, chains) {
  variables <- dimnames(bml_raw_sims_array(jags.out))[[3]]
  reported <- bml_reported_draw_variables(reg.table, variables)
  if (length(reported) == 0L) {
    return(tibble::tibble(
      node = character(), Parameter = character(), rhat = numeric(),
      ess_bulk = numeric(), ess_tail = numeric(), mcse_mean = numeric(),
      mcse_sd = numeric(), convergence = logical()
    ))
  }

  draws <- bml_select_draws(jags.out, reported)
  diagnostics <- suppressWarnings(
    posterior::summarise_draws(
      draws,
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail,
      mcse_mean = posterior::mcse_mean,
      mcse_sd = posterior::mcse_sd
    )
  )
  names(diagnostics)[names(diagnostics) == "variable"] <- "node"
  labels <- reg.table$Parameter[match(diagnostics$node, rownames(reg.table))]
  threshold <- 100 * chains
  finite <- is.finite(diagnostics$rhat) &
    is.finite(diagnostics$ess_bulk) &
    is.finite(diagnostics$ess_tail)

  diagnostics |>
    dplyr::mutate(
      Parameter = labels,
      convergence = finite &
        .data$rhat <= 1.01 &
        .data$ess_bulk >= threshold &
        .data$ess_tail >= threshold
    ) |>
    dplyr::select(
      node, Parameter, rhat, ess_bulk, ess_tail, mcse_mean, mcse_sd, convergence
    )
}

bml_available_capabilities <- function(object) {
  object$storage$capabilities %||% character()
}

bml_require_capabilities <- function(object, capabilities, method) {
  if (!inherits(object, "bml")) {
    stop("`object` must be a fitted bml model.", call. = FALSE)
  }
  available <- bml_available_capabilities(object)
  if (!is.null(object$jags.out)) {
    available <- union(available, .bml_monitor_capabilities)
  }
  missing <- setdiff(capabilities, available)
  if (length(missing) > 0L) {
    capabilities <- unique(capabilities)
    spec <- paste(capabilities, collapse = "\", \"")
    setting <- if (length(capabilities) == 1L) {
      paste0("monitor = \"", capabilities, "\"")
    } else {
      paste0("monitor = c(\"", spec, "\")")
    }
    stop(
      method, "() requires posterior capability/capabilities \"", spec,
      "\". Refit with ", setting, ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

bml_diagnostic_overview <- function(object) {
  diagnostics <- object$diagnostics
  if (is.null(diagnostics) || nrow(diagnostics) == 0L) {
    return(list(
      max_rhat = NA_real_, min_ess_bulk = NA_real_, min_ess_tail = NA_real_,
      n_flagged = NA_integer_, convergence = NA_character_
    ))
  }
  finite_max <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
  finite_min <- function(x) if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_
  n_flagged <- sum(!diagnostics$convergence)
  list(
    max_rhat = finite_max(diagnostics$rhat),
    min_ess_bulk = finite_min(diagnostics$ess_bulk),
    min_ess_tail = finite_min(diagnostics$ess_tail),
    n_flagged = n_flagged,
    convergence = if (n_flagged == 0L) "OK" else paste0("WARNING (", n_flagged, " flagged)")
  )
}

bml_draws_array <- function(object, capabilities = NULL, method = "as_draws") {
  if (!is.null(capabilities)) {
    bml_require_capabilities(object, capabilities, method)
  }
  if (!is.null(object$draws)) {
    return(posterior::as_draws_array(object$draws))
  }
  if (!is.null(object$jags.out) &&
      !is.null(object$jags.out$BUGSoutput$sims.array)) {
    return(posterior::as_draws_array(object$jags.out$BUGSoutput$sims.array))
  }
  stop(
    "No posterior draws were retained. Refit with an appropriate `monitor` capability, ",
    "for example monitor = \"parameters\".",
    call. = FALSE
  )
}

bml_draws_matrix <- function(object, capabilities = NULL, method = "as_draws") {
  posterior::as_draws_matrix(bml_draws_array(object, capabilities, method))
}
