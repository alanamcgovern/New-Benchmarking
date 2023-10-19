
library(LaplacesDemon)
library(Rfast)
library(MASS)
library(INLA)
library(VGAM)

# Simulate data for MCMC (fixed spatial effect)  -------

simdat_fsp <- function(n_areas,n_clusters,n_births,b){
  #observed data
  dat <- data.frame(
    # (true) number of births in each cluster
    N = rpois(sum(n_clusters),n_births),
    # admin area of each cluster
    A = c(sapply(1:n_areas,function(i){rep(i,n_clusters[i])})))
  dat$Y <- sapply(1:nrow(dat),function(i){rpois(1,dat$N[i]*exp(b[dat$A[i]]+rnorm(1,0,0.05)))})
  dat$cluster <- 1:nrow(dat)
  return(dat)
}

# Simulate data for MCMC (IID spatial effect) -------  

simdat_iid <- function(n_areas,n_clusters,n_births,a,tau,b=NULL){
  
  if(is.null(b)==T){
    #random effects
    b <- Rfast::rmvnorm(1,rep(0,n_areas),diag(1/tau,n_areas))
  }
  
  #observed data
  dat <- data.frame(
    # (true) number of births in each cluster
    N = rpois(sum(n_clusters),n_births),
    # admin area of each cluster
    A = c(sapply(1:n_areas,function(i){rep(i,n_clusters[i])})))
  dat$Y <- sapply(1:nrow(dat),function(i){rpois(1,dat$N[i]*exp(a + b[dat$A[i]]+rnorm(1,0,0.05)))})
  dat$cluster <- 1:nrow(dat)
  return(list(dat=dat,b=b,tau=tau))
}

# Simulate data for MCMC (BYM2 spatial effect) -------  

#simulate data
simdat_bym2sp <- function(n_areas, # number of areas (integer)
                          n_clusters=NULL, # number of clusters in each area (area-length vector)
                          n_births=NULL, # average number of births in each cluster (integer)
                          n_clusters_urban=NULL, # number of urban clusters in each area (area-length vector)
                          n_clusters_rural=NULL, # number of rural clusters in each area (area-length vector)
                          n_births_urban=NULL, # average number of births in each urban cluster (integer)
                          n_births_rural=NULL, # average number of births in each rural cluster (integer)
                          a=NULL, Var_b,b=NULL,
                          alphaU=NULL, alphaR=NULL, e=NULL){
  
  if(is.null(b)==T){
    #random effects
    b <- Rfast::rmvnorm(1,rep(0,n_areas),Var_b)
  }
  
  #generate without stratification
  if(!is.null(n_births)){
    #observed data
    dat <- data.frame(
      # (true) number of births in each cluster
      N = rpois(sum(n_clusters),n_births),
      # admin area of each cluster
      A = c(sapply(1:n_areas,function(i){rep(i,n_clusters[i])})))
    dat$Y <- sapply(1:nrow(dat),function(i){rpois(1,dat$N[i]*exp(a + b[dat$A[i]]+rnorm(1,0,0.05)))})
    dat$cluster <- 1:nrow(dat)
    return(list(dat=dat,b=b, Var_b=Var_b))
    
  #generate with stratification
  }else{
    dat <- data.frame(
      # (true) number of births in each cluster
      N = c(rpois(sum(n_clusters_urban),n_births_urban),rpois(sum(n_clusters_rural),n_births_rural)),
      # urban or rural strata (U=1 if urban, otherwise 0)
      U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
      # admin area of each cluster
      A = c(unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_rural[i])}))))
    dat$Y <- sapply(1:nrow(dat),function(i){rpois(1,dat$N[i]*exp(alphaU*dat$U[i] + alphaR*(1-dat$U[i]) + b[dat$A[i]] + e[i] +rnorm(1,0,0.05)))})
    dat$cluster <- 1:nrow(dat)
    return(list(dat=dat,b=b, Var_b=Var_b, e=e))
  }
}


# Simulate observed data (Zs) from ground truth data (Ys) -------  
# SRS sample (combining all clusters in an area together and then randomly sampling n_births_samp births from each area)
get_srs <- function(truth_dat, #data simulated from one of above functions (data frame)
                    n_births_samp, #number of births sampled from each area (area-length vector)
                    area_label){ #column indicated area level we are sampling at
  obs_dat <- NULL
  area_col <- which(names(truth_dat)==area_label)
  areas <- unique(truth_dat[,area_col])
  for(area in areas){
    N_area <- sum(truth_dat[truth_dat[,area_col]==area,]$N)
    y_area <- sum(truth_dat[truth_dat[,area_col]==area,]$Y)
    z_area <- rhyper(1,y_area,N_area-y_area,n_births_samp[area])
    
    obs_dat <- rbind(obs_dat,c(N_area,n_births_samp[area],area,y_area,z_area))
  }
  obs_dat <- data.frame(obs_dat)
  colnames(obs_dat) <- c('N','n',area_label,'Y','Z')
  return(obs_dat)
}










