#' @title summary() method for an bml object
#'
#' @description summary() method for an bml object
#' 
#' @return Returns a table with regression results.
#'
#' @examples data(coalgov)
#' m1 <- bml(Surv(govdur, earlyterm, govmaxdur) ~ 1 + mm(id(pid, gid), mmc(fdep), mmw(w ~ 1/n, constraint=1)) + majority + hm(id=cid, name=cname, type=RE, showFE=F),
#'           family="Weibull", monitor=T, data=coalgov)
#' summary(m1)
#'
#' @exportS3Method summary bml
#' @author Benjamin Rosche <benrosche@@nyu.edu>

summary.bml <- function(bml, r=3) {

  rounded_table <- bml$reg.table %>% dplyr::mutate(across(where(is.numeric), ~round(.,r)))

  # Preserve metadata
  attr(rounded_table, "estimate_type") <- attr(bml$reg.table, "estimate_type")
  attr(rounded_table, "credible_interval") <- attr(bml$reg.table, "credible_interval")
  attr(rounded_table, "DIC") <- attr(bml$reg.table, "DIC")
  attr(rounded_table, "level_spec") <- attr(bml$reg.table, "level_spec")

  # Add class for custom printing
  class(rounded_table) <- c("bml_summary", class(rounded_table))

  return(rounded_table)

}

#' @exportS3Method print bml_summary
print.bml_summary <- function(x, ...) {
  # Print header with metadata
  estimate_type <- attr(x, "estimate_type")
  ci_info <- attr(x, "credible_interval")
  dic_value <- attr(x, "DIC")
  level_spec <- attr(x, "level_spec")

  if (!is.null(estimate_type)) {
    cat("Estimates:", estimate_type, "\n")
  }
  if (!is.null(ci_info)) {
    cat("Uncertainty:", ci_info, "\n")
  }
  if (!is.null(level_spec)) {
    cat("\nLevel specification:\n")
    cat(level_spec, "\n")
  }
  if (!is.null(estimate_type) || !is.null(ci_info) || !is.null(level_spec)) {
    cat("\n")
  }

  # Remove custom class and print as data frame
  class(x) <- setdiff(class(x), "bml_summary")
  print(x, ...)

  # Print DIC at the bottom
  if (!is.null(dic_value)) {
    cat("\nModel fit: DIC =", dic_value, "\n")
  }
}




