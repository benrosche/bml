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
  
  return( bml$reg.table %>% dplyr::mutate(across(where(is.numeric), ~round(.,r))) )
  
}




