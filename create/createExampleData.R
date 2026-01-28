library(dplyr)
library(haven)
library(labelled)
library(MASS)
library(lqmm)
library(bml)
library(Hmisc)

## ATTENTION! THIS IS JUST A COPY. THE MAIN VERSION OF THIS FILE IS AND SHOULD BE KEPT IN GOVSURVIVAL

# ================================================================================================ #
# Function to create Weibull data
# ================================================================================================ #

crWeib <- function(N, lambda, rho, LP, rateC) {
  
  # Bender et al 2005: Generating survival times to simulate Cox proportional hazards models
  # https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.2059
  # https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring
  
  # Function arguments:
  # N = sample size    
  # lambda = scale parameter in h0()
  # rho = shape parameter in h0()
  # LP = linear predictor
  # rateC = rate parameter of the exponential distribution of C
  
  # Weibull latent event times
  v <- runif(N, 0, 1)
  Tlat <- (- log(v) / (lambda * exp(LP)))^(1 / rho)
  
  # censoring times
  C <- rexp(n=N, rate=rateC)
  
  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)
  
  return(cbind(survtime=time, event=status))
}

# ================================================================================================ #
# Function to create data
# ================================================================================================ #

crDat <- function(party=3, gov=c(1,3), country=3, weight=c(0,0,0), Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3), mm=F, level=1, missing=FALSE, transform=FALSE, seed=NULL) {
  
  # party=3; gov=c(1,3); country=3; weight=c(0,0,0); Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3); mm=F; level=1; transform=FALSE; missing=T; seed=NULL
  
  if(!level %in% c(1,2)) stop("Level can only be 1 or 2")
  
  # Function to center or standardize data
  cen_std <- function(x, transform) { 
    if(transform=="center") {
      if(is.numeric(x) & dim(table(x))>2) x-mean(x) else x 
    } else if(transform=="std") {
      if(is.numeric(x) & dim(table(x))>2) (x-mean(x))/sqrt(var(x)) else x 
    } else if(transform=="std2") {
      if(is.numeric(x) & dim(table(x))>2) (x-mean(x))/(2*sqrt(var(x))) else x 
    } else {
      x
    }
  }
  
  # Create X ------------------------------------------------------------------------------------- #
  
  # Sigma = covariance matrix of party, gov, and country random effects 
  
  if(!is.null(seed)) set.seed(seed)
  
  # Load data 
  dat.party   <- read_dta("./data/Rosche2017-party.dta")
  dat.gov     <- read_dta("./data/Rosche2017-government.dta")
  dat.country <- read_dta("./data/Rosche2017-country.dta")
  
  # Variable selection
  dat <- 
    dat.party %>% 
    dplyr::select(pid, partyname, gid, cid, prime, partyfam, rile, pipd, finance1, pseatrel) %>% dplyr::rename(pname=partyname, pfam=partyfam, ipd=pipd, fdep=finance1) %>%
    dplyr::inner_join(dat.gov %>% 
                        dplyr::select(gid, gname, cid, country, gstart, gend, nPG, majority, mwc, rilegov2) %>% 
                        dplyr::rename(n=nPG, cname=country, hetero=rilegov2), by=c("gid", "cid")) %>% # add gov vars
    dplyr::inner_join(dat.country %>% dplyr::select(cid, investiture, pmpower), by=c("cid")) %>% # add country vars
    dplyr::relocate(c(gname), .after=gid) %>%
    dplyr::relocate(c(cname, gstart, gend, n), .after=cid)
  var_label(dat) <- c("Unique party ID", "Party name",
                      "Unique government ID", "Government name", 
                      "Unique country ID", "Country name", 
                      "Government start date", "Government end date", 
                      "# government parties", "Prime minister party", "Party family", "Right-left position", "Intra-party democracy", "Financial dependency", "Party's relative seat share within coalition", 
                      "Majority government", "Minimal winning coalition", "SD(rile) of goverment / SD(rile) of parliament", 
                      "Investiture vote", "Prime ministerial powers")
  
  # Control the degree of multiple memberships?
  if(!isFALSE(mm)) {
    n <- dim(dat)[1]
    pid <- rep(1:n, each=mm)[1:n]
    dat <- dat %>% mutate(pid=!!pid)
  }
  
  # Create RE
  
  if(det(Sigma)<0) {
    Sigma <- make.positive.definite(Sigma)
    message("Sigma changed")
  }
  
  re.party   <- mvrnorm(n = length(unique(dat %>% .$pid)), mu=c(0,0,0), Sigma^2)[,1]
  re.gov     <- mvrnorm(n = length(unique(dat %>% .$gid)), mu=c(0,0,0), Sigma^2)[,2]
  re.country <- mvrnorm(n = length(unique(dat %>% .$cid)), mu=c(0,0,0), Sigma^2)[,3]
  
  # Add RE
  dat <- 
    dat %>%
    dplyr::group_by(pid) %>% dplyr::mutate(mmid = cur_group_id()) %>% dplyr::ungroup() %>%  dplyr::mutate(re.party=re.party[mmid]) %>%
    dplyr::group_by(gid) %>% dplyr::mutate(mainid = cur_group_id()) %>% dplyr::ungroup() %>%  dplyr::mutate(re.gov=re.gov[mainid]) %>%
    dplyr::group_by(cid) %>% dplyr::mutate(hmid = cur_group_id()) %>% dplyr::ungroup() %>%  dplyr::mutate(re.country=re.country[hmid]) %>%
    dplyr::select(-mmid, -mainid, -hmid)
  
  # Create Y ------------------------------------------------------------------------------------- #
  
  # Create linear predictor
  crLP <- 
    dat %>%
    dplyr::inner_join(dat.gov %>% 
                        dplyr::select(gid, cid, event_wkb, dur_wkb, maxdur_wkb) %>% 
                        dplyr::rename(earlyterm=event_wkb, govdur=dur_wkb, govmaxdur=maxdur_wkb), by=c("gid", "cid")) %>%
    dplyr::mutate(govmaxdur = case_when(govdur>govmaxdur ~ govdur,
                                        TRUE ~ govmaxdur)) %>%
    dplyr::mutate(partyeffect=party*fdep+re.party) %>%
    dplyr::mutate(w=1/n^exp(-(weight[1]*pseatrel+weight[2]*hetero+weight[3]*pmpower)), ng=max(gid), w=w*ng/sum(w)) %>%
    dplyr::group_by(gid) %>%
    dplyr::mutate(aggpartyeffect=sum(w*partyeffect)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(LP=(gov[1]+gov[2]*majority+re.gov) + (aggpartyeffect) + (country*investiture+re.country)) %>%
    dplyr::select(-partyeffect, -aggpartyeffect, -re.party, -re.gov, -re.country) %>%
    dplyr::rename(sim.w=w)
  
  # Create outcomes
  finalDat <-
    crLP %>%
    dplyr::mutate(sim.y=LP) %>%
    cbind(crWeib(N=dim(crLP)[1], lambda=0.01, rho=1, LP=crLP %>% .$LP, rateC=0.001)) %>% 
    dplyr::rename(sim.st=survtime, sim.e=event) %>%
    dplyr::select(-LP, -ng)
  
  var_label(finalDat)$earlyterm <- "Discretionary early termination"
  var_label(finalDat)$govdur <- "Government duration"
  var_label(finalDat)$govmaxdur <- "Maximum possible government duration"
  var_label(finalDat)$sim.w <- "Simulated weights"
  var_label(finalDat)$sim.y <- "Simulated linear outcome"
  var_label(finalDat)$sim.st <- "Simulated survival time"
  var_label(finalDat)$sim.e <- "Simulated event status"
  
  # Induce missingness 
  if(missing==TRUE) {
    
    # Induce missingness based on ipd: the higher ipd, the more likely fdep is observed.
    ipd_dist <- finalDat %>% .$ipd
    is_observed <- as.logical(rbinom(dim(finalDat)[1],1, sapply(finalDat %>% .$ipd, function(x) { sum(ipd_dist<=x)/length(ipd_dist) })))
    
    crMissing <-
      finalDat %>% 
      dplyr::mutate(fdep.mi=case_when(is_observed==TRUE ~ fdep,
                                      TRUE ~ NA_real_)) %>%
      dplyr::mutate(fdep.imp=fdep.mi) %>%
      dplyr::relocate(c(fdep.mi, fdep.imp), .after = fdep)
    
    # Impute values again with Hmisc: aregImpute
    imputed <- aregImpute(fdep.mi~ipd, n.impute = 1, type = "pmm", data = crMissing)$imputed$fdep.mi
    crMissing$fdep.imp[which(is.na(crMissing$fdep.im))] <- as.vector(imputed)
    
    # MICE: Impute fdep.mi and save in fdep.imp
    # imputed_mice <- complete(mice(data=crY %>% dplyr::select(gid, pid, ipd, fdep.mi), m = 1, defaultMethod = "pmm"))
    # crY$fdep.imp <- imputed_mice$fdep.mi
    
    finalDat <- crMissing
    
    # Add labels and add vars to aggregation list
    var_label(finalDat)$fdep.mi <- "fdep with missing values"
    var_label(finalDat)$fdep.imp <- "fdep, missing values imputed"
    
    partyvars <- append(partyvars, c("fdep.mi", "fdep.imp"))
  }
  
  # Aggregate to second level 
  if(level==2) {
    finalDat <- 
      finalDat %>%
      dplyr::group_by(gid) %>%
      dplyr::mutate(across(!!partyvars, ~mean(.x, na.rm = TRUE))) %>%
      dplyr::filter(row_number()==1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(pseatrel=finalDat %>% dplyr::group_by(gid) %>% summarise(pseatrel=max(pseatrel)) %>% .$pseatrel) %>%
      dplyr::select(-pid, -pname, -prime, -pfam, -sim.w) 
    
    var_label(finalDat)$ipd <- "Mean: Intra-party democracy"
    var_label(finalDat)$fdep <- "Mean: Financial dependency"
    var_label(finalDat)$fdep.mi <- "Mean: fdep with missing values"
    var_label(finalDat)$fdep.imp <- "Mean: fdep, missing values imputed"
    var_label(finalDat)$pseatrel <- "Max: Party's relative seat share within coalition"
    
  }
  
  return(finalDat)
  
}

# ================================================================================================ #
# Save data for bml() package
# ================================================================================================ #

# coalgov <- crDat(party=3, gov=c(0,3), country=3, weight=c(1,1,1), Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3), level=1, seed=1)
# save(coalgov, file="data/coalgov.RData")
