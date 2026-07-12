# ================================================================================================ #
# prior(): structured, class-targeted priors translated to JAGS
# ================================================================================================ #

#' Structured prior specifications
#'
#' @description
#' Specify priors brms-style: a distribution call plus a parameter \code{class},
#' optionally narrowed to a single coefficient (\code{coef}) or block
#' (\code{group}). Distributions are written on the natural (standard-deviation)
#' scale and translated to JAGS's precision parameterization automatically.
#'
#' \preformatted{
#' prior(normal(0, 10),    class = "b")                      # all fixed coefficients
#' prior(normal(0, 1),     class = "b", coef = "education")  # one coefficient
#' prior(exponential(1),   class = "sd")                     # RE standard deviations
#' prior(lkj(2),           class = "cor")                    # RE correlation (re(cor=TRUE))
#' prior(student_t(3,0,5), class = "Intercept")
#' prior(normal(0, 1),     class = "w")                      # weight-model parameters
#' prior(uniform(0, 1),    class = "fn", coef = "c")         # fn shape parameters
#' }
#'
#' Combine several with \code{c()} and pass as \code{bml(..., prior = c(...))}.
#' Raw JAGS strings remain a power-user escape hatch: character elements such as
#' \code{"b.fn.1[1] ~ dnorm(0, 0.1)"} are passed through untranslated.
#'
#' @param dist A distribution call. Supported: \code{normal(mu, sigma)},
#'   \code{student_t(nu, mu, sigma)}, \code{cauchy(mu, sigma)},
#'   \code{exponential(rate)}, \code{gamma(shape, rate)}, \code{uniform(lower,
#'   upper)}, \code{lkj(eta)} (class \code{"cor"} only). Scale parameters are
#'   standard deviations, not precisions.
#' @param class Parameter class. One of:
#'   \itemize{
#'     \item \code{"Intercept"}: the main intercept
#'     \item \code{"b"}: fixed coefficients (main-formula terms and block
#'       features); narrow with \code{coef} (term label, e.g. \code{"education"}
#'       or a feature name like \code{"A_x"})
#'     \item \code{"sd"}: random-effect standard deviations; narrow with
#'       \code{group} (\code{"mm"}/\code{"hm"} or a block id like \code{"mm.1"})
#'       and \code{coef} (\code{"Intercept"} or a slope variable)
#'     \item \code{"cor"}: the RE intercept-slope correlation (only exists when
#'       the user opts into \code{re(..., cor = TRUE)}) — use \code{lkj(eta)}
#'     \item \code{"sigma"}: the Gaussian residual SD
#'     \item \code{"shape"}: the Weibull shape parameter
#'     \item \code{"w"}: weight-model parameters (free symbols in \code{w()})
#'     \item \code{"fn"}: \code{fn()} shape parameters (\code{c}, \code{kappa},
#'       \code{p}, DSL free symbols); narrow with \code{coef}
#'   }
#' @param coef Optional coefficient/parameter name to narrow the prior to.
#' @param group Optional block identifier (e.g. \code{"mm.1"}, \code{"hm.2"},
#'   or just \code{"mm"}/\code{"hm"}).
#'
#' @return A \code{bml_prior} data frame (one row per prior); combine with \code{c()}.
#'
#' @seealso \code{\link{get_prior}}, \code{\link{bml}}
#' @export
prior <- function(dist, class = "b", coef = "", group = "") {
  dist_call <- substitute(dist)
  dist_str <- paste(deparse(dist_call, width.cutoff = 500), collapse = "")
  new_bml_prior(dist_str, class = class, coef = coef, group = group)
}

new_bml_prior <- function(prior, class, coef = "", group = "") {
  out <- data.frame(
    prior = prior,
    class = class,
    coef = coef %||% "",
    group = group %||% "",
    stringsAsFactors = FALSE
  )
  class(out) <- c("bml_prior", "data.frame")
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
c.bml_prior <- function(...) {
  parts <- list(...)
  parts <- lapply(parts, function(p) {
    if (inherits(p, "bml_prior")) return(p)
    if (is.character(p)) return(new_bml_prior(p, class = "raw"))
    stop("prior specifications must be prior() calls or raw JAGS strings.", call. = FALSE)
  })
  out <- do.call(rbind, parts)
  class(out) <- c("bml_prior", "data.frame")
  out
}

#' @export
print.bml_prior <- function(x, ...) {
  cat("bml prior specification:\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

# ------------------------------------------------------------------------------------------------ #
# Distribution translation: user scale (SD) -> JAGS scale (precision)
# ------------------------------------------------------------------------------------------------ #

# Translate a user distribution string (e.g. "normal(0, 10)") into a JAGS
# distribution string (e.g. "dnorm(0, 0.01)"). `positive = TRUE` appends a
# T(0,) truncation for SD-type parameters where the distribution has support
# below zero.
translate_dist <- function(dist_str, positive = FALSE) {

  expr <- tryCatch(parse(text = dist_str)[[1]], error = function(e) NULL)
  if (is.null(expr) || !is.call(expr)) {
    stop("Could not parse prior distribution: '", dist_str, "'", call. = FALSE)
  }

  fname <- as.character(expr[[1]])
  args <- as.list(expr)[-1]
  num <- function(i) {
    v <- tryCatch(eval(args[[i]], envir = baseenv()), error = function(e) NULL)
    if (!is.numeric(v) || length(v) != 1) {
      stop("Prior '", dist_str, "': argument ", i, " must be a single number.", call. = FALSE)
    }
    v
  }

  trunc_pos <- function(s) if (positive) paste0(s, "T(0,)") else s

  out <- switch(
    fname,
    "normal" = {
      mu <- num(1); sigma <- num(2)
      if (sigma <= 0) stop("normal(): sigma must be > 0.", call. = FALSE)
      trunc_pos(sprintf("dnorm(%s, %s)", format(mu), format(1 / sigma^2)))
    },
    "student_t" = {
      nu <- num(1); mu <- num(2); sigma <- num(3)
      if (sigma <= 0) stop("student_t(): sigma must be > 0.", call. = FALSE)
      trunc_pos(sprintf("dt(%s, %s, %s)", format(mu), format(1 / sigma^2), format(nu)))
    },
    "cauchy" = {
      mu <- num(1); sigma <- num(2)
      if (sigma <= 0) stop("cauchy(): sigma must be > 0.", call. = FALSE)
      trunc_pos(sprintf("dt(%s, %s, 1)", format(mu), format(1 / sigma^2)))
    },
    "exponential" = {
      rate <- num(1)
      sprintf("dexp(%s)", format(rate))
    },
    "gamma" = {
      shape <- num(1); rate <- num(2)
      sprintf("dgamma(%s, %s)", format(shape), format(rate))
    },
    "uniform" = {
      lo <- num(1); hi <- num(2)
      sprintf("dunif(%s, %s)", format(lo), format(hi))
    },
    "lkj" = {
      # Handled specially for class "cor" (2x2): rho = 2*z - 1, z ~ dbeta(eta, eta).
      eta <- num(1)
      sprintf("__lkj__(%s)", format(eta))
    },
    # Already-JAGS distributions pass through (dnorm, dgamma, ...):
    if (grepl("^d[a-z.]+$", fname)) dist_str else
      stop("Unsupported prior distribution: '", fname, "'. Supported: normal, student_t, ",
           "cauchy, exponential, gamma, uniform, lkj (or a raw JAGS d* call).", call. = FALSE)
  )

  out
}

# ------------------------------------------------------------------------------------------------ #
# Prior table: enumerate every settable prior of a model
# ------------------------------------------------------------------------------------------------ #

# Build the table of settable priors from the dissected model. Used by both
# get_prior() (returned to the user) and resolve_priors() (matching).
# Each row: class, coef, group, node (JAGS node name), default (JAGS dist),
# positive (needs T(0,) when overridden by an unbounded dist).
build_prior_table <- function(formula_parts, family) {

  rows <- list()
  addrow <- function(class, coef, group, node, default, positive = FALSE) {
    rows[[length(rows) + 1]] <<- data.frame(
      class = class, coef = coef, group = group, node = node,
      default = default, positive = positive, stringsAsFactors = FALSE
    )
  }

  mainvars <- formula_parts$mainvars
  mm_list <- formula_parts$mm
  hm_list <- formula_parts$hm

  # Intercept + main coefficients ------------------------------------------------------------- #
  if (length(mainvars) > 0) {
    for (x in seq_along(mainvars)) {
      v <- mainvars[x]
      if (v == "X0") {
        addrow("Intercept", "", "", paste0("b[", x, "]"), "dnorm(0, 0.0001)")
      } else {
        addrow("b", v, "", paste0("b[", x, "]"), "dnorm(0, 0.0001)")
      }
    }
  }

  # mm blocks ---------------------------------------------------------------------------------- #
  for (k in seq_along(mm_list)) {
    m <- mm_list[[k]]
    group_id <- paste0("mm.", k)

    # feature coefficients
    feats <- m$feature_labels %||% character(0)
    for (x in seq_along(feats)) {
      addrow("b", feats[x], group_id, paste0("b.fn.", k, "[", x, "]"), "dnorm(0, 0.0001)")
    }

    # weight parameters
    for (p in seq_along(m$w$params %||% character(0))) {
      addrow("w", m$w$params[p], group_id, paste0("b.w.", k, "[", p, "]"), "dnorm(0, 0.0001)")
    }

    # fn shape parameters
    for (pn in names(m$fn$est_params %||% list())) {
      addrow("fn", pn, group_id, paste0("fn.", pn, ".", k),
             m$fn$est_params[[pn]]$default %||% "dunif(-1e6, 1e6)")
    }

    # RE sd / cor
    if (!is.null(m$RE)) {
      g <- m$mmid_group %||% 1
      addrow("sd", "Intercept", group_id, paste0("sigma.mm.", g), "half_t(25)", positive = TRUE)
      for (s in seq_along(m$RE$slopes %||% character(0))) {
        addrow("sd", m$RE$slopes[s], group_id, paste0("sigma.mm.", g, ".s", s),
               "half_t(25)", positive = TRUE)
      }
      if (isTRUE(m$RE$cor) && length(m$RE$slopes %||% character(0)) > 0) {
        addrow("cor", paste0("Intercept, ", m$RE$slopes[1]), group_id,
               paste0("rho.mm.", g), "lkj(1)")
      }
    }
  }

  # hm blocks ---------------------------------------------------------------------------------- #
  for (k in seq_along(hm_list)) {
    h <- hm_list[[k]]
    group_id <- paste0("hm.", k)
    if (!is.null(h$RE)) {
      addrow("sd", "Intercept", group_id, paste0("sigma.hm.", k), "half_t(25)", positive = TRUE)
      for (s in seq_along(h$RE$slopes %||% character(0))) {
        addrow("sd", h$RE$slopes[s], group_id, paste0("sigma.hm.", k, ".s", s),
               "half_t(25)", positive = TRUE)
      }
      if (isTRUE(h$RE$cor) && length(h$RE$slopes %||% character(0)) > 0) {
        addrow("cor", paste0("Intercept, ", h$RE$slopes[1]), group_id,
               paste0("rho.hm.", k), "lkj(1)")
      }
    }
  }

  # interaction coefficients ------------------------------------------------------------------- #
  for (t in seq_along(formula_parts$interactions %||% list())) {
    addrow("b", formula_parts$interactions[[t]]$label, "",
           paste0("b.int[", t, "]"), "dnorm(0, 0.0001)")
  }

  # family parameters -------------------------------------------------------------------------- #
  if (family == "Gaussian") {
    addrow("sigma", "", "", "sigma", "half_t(25)", positive = TRUE)
  } else if (family == "Weibull") {
    addrow("shape", "", "", "shape", "dexp(0.001)")
  }

  if (length(rows) == 0) {
    return(data.frame(class = character(), coef = character(), group = character(),
                      node = character(), default = character(), positive = logical(),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' List every settable prior of a bml model
#'
#' @description
#' Returns the table of all parameters whose priors can be set via
#' \code{\link{prior}}, together with their default priors — the bml analogue of
#' \code{brms::get_prior()}. Use it to discover valid \code{class}/\code{coef}/
#' \code{group} combinations before overriding.
#'
#' @param formula The model formula (same as for \code{\link{bml}}).
#' @param data The data (member-level long format).
#' @param family Model family (function, object, or string). Default \code{gaussian()}.
#'
#' @return A data frame with columns \code{class}, \code{coef}, \code{group},
#'   \code{node} (the internal JAGS node), and \code{default}.
#'
#' @examples
#' \dontrun{
#' get_prior(
#'   y ~ x + mm(id = id(pid, gid), vars = vars(z), w = w(~ 1/n)),
#'   data = dat, family = gaussian()
#' )
#' }
#'
#' @seealso \code{\link{prior}}, \code{\link{bml}}
#' @export
get_prior <- function(formula, data, family = stats::gaussian()) {
  fam <- validate_family(family)
  formula_parts <- dissectFormula(formula, fam$family, data)
  tab <- build_prior_table(formula_parts, fam$family)
  tab$positive <- NULL
  # Show the user-facing default spelling for the half-t default
  tab$default <- sub("^half_t\\((\\d+)\\)$", "student_t(1, 0, \\1)", tab$default)
  tab
}

# ------------------------------------------------------------------------------------------------ #
# Resolution: match user priors against the table -> node -> JAGS dist string
# ------------------------------------------------------------------------------------------------ #

# Returns list(overrides = named list node -> JAGS dist string,
#              raw = character vector of raw JAGS prior lines)
resolve_priors <- function(prior, prior_table) {

  overrides <- list()
  raw <- character(0)

  if (is.null(prior)) {
    return(list(overrides = overrides, raw = raw))
  }

  # Legacy/raw input: plain character vector or list of strings
  if (!inherits(prior, "bml_prior")) {
    if (is.character(prior)) {
      return(list(overrides = overrides, raw = prior))
    }
    if (is.list(prior) && all(vapply(prior, is.character, logical(1)))) {
      return(list(overrides = overrides, raw = unlist(prior)))
    }
    stop("'prior' must be built with prior() (combine with c()), ",
         "or raw JAGS strings.", call. = FALSE)
  }

  for (i in seq_len(nrow(prior))) {
    p <- prior[i, ]

    if (p$class == "raw") {
      raw <- c(raw, p$prior)
      next
    }

    matches <- prior_table$class == p$class
    if (nzchar(p$coef)) matches <- matches & prior_table$coef == p$coef
    if (nzchar(p$group)) {
      # allow "mm" to match "mm.1", "mm.2", ...
      matches <- matches & (prior_table$group == p$group |
                              startsWith(prior_table$group, paste0(p$group, ".")))
    }

    if (!any(matches)) {
      stop("prior(", p$prior, ", class = \"", p$class, "\"",
           if (nzchar(p$coef)) paste0(", coef = \"", p$coef, "\""),
           if (nzchar(p$group)) paste0(", group = \"", p$group, "\""),
           ") matches no parameter. Run get_prior() to see what is settable.", call. = FALSE)
    }

    idx <- which(matches)
    for (j in idx) {
      node <- prior_table$node[j]
      jags_dist <- translate_dist(p$prior, positive = isTRUE(prior_table$positive[j]))
      overrides[[node]] <- jags_dist
    }
  }

  list(overrides = overrides, raw = raw)
}
