# ================================================================================================ #
# Temporal helpers: aggregation over membership spells (observed-panel case)
# ================================================================================================ #
#
# With observed weights and attributes, temporal summaries of a composition are pre-computable
# covariates - there is no in-model uncertainty to preserve. These helpers build the matched
# spell table (union roster across two snapshots, per-snapshot normalized weights) and return
# group-level covariates:
#
#   bml_delta():    dA_x = A_{x,t} - A_{x,t-1}, plus its matched decomposition
#                   (reallocation + attribute change + entry/exit; adding-up holds exactly)
#   bml_realloc():  the reallocation component sum_k dw_k xbar_k, and its intensity E((dw)^2)
#   bml_turnover(): entry/exit weight mass and membership Jaccard distance
#
# The in-model spell path (estimated weights dw(omega); reallocation toward members with high
# random effects sum_k dw_k u_k) is future work; see the package vignettes.
#
# ================================================================================================ #

# Internal: build the spell table for one pair of snapshots.
# Returns one row per (group, member) on the union roster with w1, w2 (normalized within
# group and snapshot; 0 when absent), x1, x2 (NA when absent), dw, wbar, and status.
build_spell_table <- function(data, group, member, time, weight = NULL, x = NULL,
                              t1, t2) {

  d <- dplyr::ungroup(data)

  for (v in c(group, member, time, weight, x)) {
    if (!is.null(v) && !v %in% names(d)) {
      stop("Variable '", v, "' not found in data.", call. = FALSE)
    }
  }

  snap <- function(tt, suffix) {
    di <- d[d[[time]] == tt, , drop = FALSE]
    if (nrow(di) == 0) {
      stop("No rows found for time value '", tt, "'.", call. = FALSE)
    }
    out <- data.frame(
      group = di[[group]],
      member = di[[member]],
      stringsAsFactors = FALSE
    )
    wv <- if (!is.null(weight)) di[[weight]] else rep(1, nrow(di))
    # normalize within group and snapshot: each snapshot's weights sum to 1,
    # so dw sums to 0 over the union roster (reallocation is a pure shift term)
    out$w <- as.numeric(wv) / stats::ave(as.numeric(wv), out$group, FUN = sum)
    if (!is.null(x)) out$x <- di[[x]]
    names(out)[names(out) == "w"] <- paste0("w", suffix)
    if (!is.null(x)) names(out)[names(out) == "x"] <- paste0("x", suffix)
    out
  }

  s1 <- snap(t1, "1")
  s2 <- snap(t2, "2")

  spells <- merge(s1, s2, by = c("group", "member"), all = TRUE)
  spells$w1[is.na(spells$w1)] <- 0
  spells$w2[is.na(spells$w2)] <- 0
  spells$dw <- spells$w2 - spells$w1
  spells$wbar <- (spells$w1 + spells$w2) / 2
  spells$status <- ifelse(spells$w1 == 0, "entered",
                          ifelse(spells$w2 == 0, "exited", "stayed"))
  spells
}

#' Temporal composition change: matched decomposition of a feature difference
#'
#' @description
#' Computes the change in a group's weighted mean attribute between two
#' snapshots, \eqn{\Delta A_x = A_{x,t_2} - A_{x,t_1}}, together with its
#' matched (shift-share / Foster-Haltiwanger-Krizan style) decomposition over
#' membership spells:
#' \deqn{\Delta A_x = \underbrace{\sum_k \Delta w_k \bar x_k}_{reallocation}
#'   + \underbrace{\sum_k \bar w_k \Delta x_k}_{attribute\ change}
#'   + entry/exit\ terms}
#' The components sum to \eqn{\Delta A_x} exactly (adding-up identity). The
#' returned group-level covariates can be entered into a \code{\link{bml}}
#' main formula.
#'
#' @param data Member-level long data covering both snapshots.
#' @param group,member,time Column names (strings): group id, member id, and
#'   the time variable.
#' @param weight Column name (string) of the raw weight, or \code{NULL} for
#'   equal weights.
#' @param x Column name (string) of the member attribute.
#' @param t1,t2 The two time values to compare.
#'
#' @return A data frame with one row per group: \code{dA} (total change),
#'   \code{realloc} (reallocation among stayers), \code{attr_change}
#'   (attribute change among stayers), \code{entry_exit} (net entry/exit
#'   contribution).
#'
#' @details
#' Weights are normalized within each group and snapshot, so \eqn{\Delta w}
#' sums to zero over the union roster. For stayers, \eqn{\bar x} and
#' \eqn{\bar w} are the across-snapshot means; entering/exiting members
#' contribute their full \eqn{w x} term to \code{entry_exit}.
#'
#' @seealso \code{\link{bml_realloc}}, \code{\link{bml_turnover}}
#' @export
bml_delta <- function(data, group, member, time, weight = NULL, x, t1, t2) {

  spells <- build_spell_table(data, group, member, time, weight, x, t1, t2)

  by_group <- split(spells, spells$group)
  out <- lapply(names(by_group), function(gname) {
    s <- by_group[[gname]]
    stay <- s$status == "stayed"

    A2 <- sum(s$w2 * s$x2, na.rm = TRUE)
    A1 <- sum(s$w1 * s$x1, na.rm = TRUE)
    dA <- A2 - A1

    realloc <- sum(s$dw[stay] * (s$x1[stay] + s$x2[stay]) / 2)
    attr_change <- sum(s$wbar[stay] * (s$x2[stay] - s$x1[stay]))
    entry <- sum(s$w2[s$status == "entered"] * s$x2[s$status == "entered"], na.rm = TRUE)
    exitc <- -sum(s$w1[s$status == "exited"] * s$x1[s$status == "exited"], na.rm = TRUE)

    data.frame(
      group = gname,
      dA = dA,
      realloc = realloc,
      attr_change = attr_change,
      entry_exit = entry + exitc,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  names(out)[1] <- group
  rownames(out) <- NULL
  out
}

#' Reallocation of membership weights toward high-attribute members
#'
#' @description
#' The headline temporal micro-macro term: did the group shift weight mass
#' toward members with high \eqn{x}? Computes \eqn{\sum_k \Delta w_k x_k} over
#' the union roster (using each member's observed attribute) plus the
#' reallocation \emph{intensity} \eqn{\sum_k \bar w_k (\Delta w_k)^2}-style
#' spread measures.
#'
#' @inheritParams bml_delta
#'
#' @return A data frame with one row per group: \code{realloc_x}
#'   (\eqn{\sum_k \Delta w_k x_k}; reallocation toward high-\code{x} members)
#'   and \code{realloc_intensity} (\eqn{\sum_k (\Delta w_k)^2}; how much
#'   reshuffling happened, regardless of direction).
#'
#' @seealso \code{\link{bml_delta}}, \code{\link{bml_turnover}}
#' @export
bml_realloc <- function(data, group, member, time, weight = NULL, x, t1, t2) {

  spells <- build_spell_table(data, group, member, time, weight, x, t1, t2)

  xbar <- ifelse(is.na(spells$x1), spells$x2,
                 ifelse(is.na(spells$x2), spells$x1, (spells$x1 + spells$x2) / 2))

  by_group <- split(data.frame(dw = spells$dw, xbar = xbar), spells$group)
  out <- lapply(names(by_group), function(gname) {
    s <- by_group[[gname]]
    data.frame(
      group = gname,
      realloc_x = sum(s$dw * s$xbar, na.rm = TRUE),
      realloc_intensity = sum(s$dw^2),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  names(out)[1] <- group
  rownames(out) <- NULL
  out
}

#' Membership turnover between two snapshots
#'
#' @description
#' Set-level composition dynamics: how much of the group's membership mass
#' entered or exited between the two snapshots.
#'
#' @inheritParams bml_delta
#' @param x Not used (turnover is a property of the roster and weights alone).
#'
#' @return A data frame with one row per group: \code{entry_mass} (weight mass
#'   of entrants at \code{t2}), \code{exit_mass} (weight mass of exiters at
#'   \code{t1}), \code{jaccard_dist} (1 - |intersection| / |union| of the
#'   member sets), and \code{realloc_l1} (\eqn{\frac{1}{2}\sum_k |\Delta w_k|},
#'   the total variation distance between the two weight measures).
#'
#' @seealso \code{\link{bml_delta}}, \code{\link{bml_realloc}}
#' @export
bml_turnover <- function(data, group, member, time, weight = NULL, t1, t2, x = NULL) {

  spells <- build_spell_table(data, group, member, time, weight, x = NULL, t1, t2)

  by_group <- split(spells, spells$group)
  out <- lapply(names(by_group), function(gname) {
    s <- by_group[[gname]]
    n_union <- nrow(s)
    n_inter <- sum(s$status == "stayed")
    data.frame(
      group = gname,
      entry_mass = sum(s$w2[s$status == "entered"]),
      exit_mass = sum(s$w1[s$status == "exited"]),
      jaccard_dist = 1 - n_inter / n_union,
      realloc_l1 = 0.5 * sum(abs(s$dw)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  names(out)[1] <- group
  rownames(out) <- NULL
  out
}
