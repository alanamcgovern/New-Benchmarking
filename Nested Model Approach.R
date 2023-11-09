library(INLA)
library(SUMMER)
library(tidyverse)
library(ggpubr)
library(Rfast)

all.ests <- NULL

country_t <- 'Guinea' # country being used
shapeFiles_folder_name <- 'shapeFiles/gadm41_GIN_shp'
beg.years <- c(2000,2004,2008,2012,2016)
end.years <- c(2003,2007,2011,2015,2018)

country_t <- 'Malawi' # country being used
shapeFiles_folder_name <- 'shapeFiles_gadm'
beg.years <- c(2000,2004,2007,2010,2013)
end.years <- c(2003,2006,2009,2012,2015)

# load data ------
data.dir <- (paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,'/'))
res.dir <- paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Results/",country_t,'/')

load(paste0(data.dir,country_t,"_cluster_dat_1frame.rda"))
dat <- mod.dat[mod.dat$age==0,]
dat$died <- dat$Y
dat$years.int <- as.integer(dat$years)
dat$years <- as.numeric(as.character(dat$years))

dat$v005 <- dat$v005/1e6
admin_key <- dat %>% dplyr::select(admin1.char,admin2.char) %>% unique()
admin_key_num <- dat %>% dplyr::select(admin1,admin2) %>% unique()

#adjacency matrix
load(paste0(data.dir,shapeFiles_folder_name,"/",country_t,"_Amat.rda"))
# make separate islands for each admin1 area
admin2.mat.nested <- matrix(0,nrow(admin2.mat),nrow(admin2.mat))
for(area in 1:nrow(admin1.mat)){
  adm2_areas_tmp <- sort(admin_key_num[admin_key_num$admin1==area,]$admin2)
  for(area2 in adm2_areas_tmp){
    if(sum(admin2.mat[area2,adm2_areas_tmp])==0){
      admin2.mat.nested[area2,adm2_areas_tmp] <- admin2.mat[area2,adm2_areas_tmp]
      }else{
       admin2.mat.nested[area2,adm2_areas_tmp] <- admin2.mat[area2,adm2_areas_tmp]/sum(admin2.mat[area2,adm2_areas_tmp])}
  }
}

load(paste0(data.dir,shapeFiles_folder_name,"/",country_t,"_Amat_names.rda"))

# get weights for aggregating from admin2 to admin 1
load(paste0(data.dir,"worldpop/adm1_weights_u1.rda"))
adm1_weights <- weight.adm1.u1
load(paste0(data.dir,"worldpop/adm2_weights_u1.rda"))
#proportion of admin2 relative to national
adm2_weights <- merge(weight.adm2.u1,admin_key,by.x='region',by.y='admin2.char')
#add to get proportion of admin1 relative to national and merge
adm1_props <- adm2_weights %>% group_by(admin1.char,years) %>% summarise(adm1_prop = sum(proportion))
adm2_to_adm1_weights_t <- left_join(adm2_weights,adm1_props,by=c('admin1.char','years'))
# get proportion of admin2 relative to admin 1
adm2_to_adm1_weights <- adm2_to_adm1_weights_t %>% group_by(admin1.char,years) %>% mutate(adm2_to_adm1_prop = proportion/adm1_prop)
adm2_to_adm1_weights <- adm2_to_adm1_weights[,c(1,3,4,6)]
colnames(adm2_to_adm1_weights)[4] <- 'proportion'

# get urban/rural weights for admin1 and admin2
adm1_UR_weights <- readRDS(paste0(res.dir,"UR/U1_fraction/admin1_u1_urban_weights.rds"))
adm1_UR_weights <- gather(adm1_UR_weights,strata,proportion,urban:rural)
adm2_UR_weights <- readRDS(paste0(res.dir,"UR/U1_fraction/admin2_u1_urban_weights.rds"))
adm2_UR_weights <- gather(adm2_UR_weights,strata,proportion,urban:rural)

## ------
for(period in 1:length(beg.years)){
# organize data for modeling -----
data_years <- beg.years[period]:end.years[period] # years of data included
ref_year <- 2014 # year used from population weights

survey_year <- max(dat$survey)
dat_t <- dat[dat$survey==survey_year & dat$years%in%data_years,]
cluster_totals <- dat_t %>% group_by(cluster) %>% summarise(Ysum=sum(Y),totalsum=sum(total))
dat_t <- left_join(dat_t,cluster_totals,by=join_by(cluster))
dat_t <- dat_t %>% select(cluster,survey,age,v005,urban,admin1,admin2,admin1.char,admin2.char,admin1.name,admin2.name,Ysum,totalsum) %>%
  rename(Y=Ysum, total=totalsum) %>% unique()
dat_t$years <- dat_t$survey <- ref_year
dat_t$died <- dat_t$Y
dat_t <- dat_t %>% mutate(strata = ifelse(urban=='urban',1,0),
                          admin1.strata = ifelse(strata==1,admin1,admin1 + max(dat_t$admin1)),
                          admin2.strata = ifelse(strata==1,admin2,admin2 + max(dat_t$admin2)))

adm1_UR_weights_t <- adm1_UR_weights[adm1_UR_weights$years==ref_year,]
adm1_UR_weights_t <- spread(adm1_UR_weights_t,strata,proportion) %>% select(-years)
adm2_to_adm1_weights_t <- adm2_to_adm1_weights[adm2_to_adm1_weights$years==ref_year,]
adm1_weights_t <- adm1_weights[adm1_weights$years==ref_year,]

# functions to draw from posterior and obtain summaries -----

getSmoothedFullPost <- function(inla_mod,data,nsim=1000,
                            admin1.weights=NULL,admin1.UR.weights=NULL){
  
  # data frame for admin 1 overall estimates
  outFullOverall <- data.frame(admin1 = 1:max(data$admin1))
  outFullOverall$admin1.char <- paste0('admin1_', outFullOverall$admin1)
  
  # sample from posterior
  cs <- inla_mod$misc$configs$contents$tag
  cs <- cs[cs != "Predictor"]
  select <- list()
  for (i in 1:length(cs)) {
    select[[i]] <- 0
    names(select)[i] <- cs[i]
  }
  
  sampFull <- INLA::inla.posterior.sample(n = nsim, result = inla_mod, intern = TRUE, selection = select)
  sampFull.draws <- matrix(NA,nsim,length(sampFull[[1]]$latent))
  for(i in 1:nsim){
    sampFull.draws[i,] <- sampFull[[i]]$latent
  }
  fields <- colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)
  
  ## matrix which specifies which columns to include and how to weight them for admin1 overall
  AA.wt.adm1 <- matrix(0,nrow = length(fields),ncol = max(data$admin1))
  # matrix which specifies which columns to include and how to weight them for national
  AA.wt.natl <- rep(0,ncol(sampFull.draws))
  
  for(k in 1:max(data$admin1)){
  # add admin1 effects (if in model)
    col.id <- which(paste0("factor(admin1)",k,":1") == fields)
    if(length(col.id)==1){
      AA.wt.adm1[col.id,] <- I(outFullOverall$admin1==k)
      AA.wt.natl[col.id] <- admin1.weights[admin1.weights$region==paste0('admin1_',k),]$proportion
    }
    
  # add weights for admin1.strata effects (if in model)
    urban.id <- which(paste0("factor(admin1.strata)",k,":1") == fields)
    if(length(urban.id)==1){
      AA.wt.adm1[urban.id,] <- I(outFullOverall$admin1==k)*admin1.UR.weights$urban[k]
      AA.wt.natl[urban.id] <- admin1.weights$proportion[k]*admin1.UR.weights$urban[k]
    }
    rural.id <- which(paste0("factor(admin1.strata)",k+max(data$admin1),":1") == fields)
    if(length(rural.id)==1){
      AA.wt.adm1[rural.id,] <- I(outFullOverall$admin1==k)*admin1.UR.weights$rural[k]
      AA.wt.natl[rural.id] <- admin1.weights$proportion[k]*admin1.UR.weights$rural[k]
    }
  }
  # add weights for urban/rural (if in model)
  rural.id <- which("factor(strata)0:1" == fields)
  if(length(rural.id)==1){
    AA.wt.adm1[rural.id,] <- admin1.UR.weights$rural
    AA.wt.natl[rural.id] <- sum(admin1.UR.weights$rural*admin1.weights$proportion)
  }
  urban.id <- which("factor(strata)1:1" == fields)
  if(length(urban.id)==1){
    AA.wt.adm1[urban.id,] <- admin1.UR.weights$urban
    AA.wt.natl[urban.id] <- sum(admin1.UR.weights$urban*admin1.weights$proportion)
  }
  # add weights for any admin2 effects (if in model)
  
  ## admin1 estimates
  AA.adm1 <- sampFull.draws %*% AA.wt.adm1
  outFullOverall$log.mean <- colMeans(AA.adm1)
  outFullOverall$mean <- colMeans(exp(AA.adm1))
  outFullOverall$median <- colMedians(exp(AA.adm1))
  outFullOverall$variance <- colVars(exp(AA.adm1))
  outFullOverall$lower <- apply(exp(AA.adm1),2,quantile,0.05)
  outFullOverall$upper <- apply(exp(AA.adm1),2,quantile,0.95)
  
  ## national estimates
  AA.natl <- sampFull.draws %*% AA.wt.natl
  outFullNatl <- data.frame(mean = colMeans(exp(AA.natl)))
  outFullNatl$log.mean <- colMeans(AA.natl)
  outFullNatl$median <- colMedians(exp(AA.natl))
  outFullNatl$variance <- colVars(exp(AA.natl))
  outFullNatl$lower <- apply(exp(AA.natl),2,quantile,0.05)
  outFullNatl$upper <- apply(exp(AA.natl),2,quantile,0.95)
  
  return(list(adm1.est = outFullOverall, natl.est = outFullNatl))
}


# fit models without temporal component (as if all data is 1 year) ------
hyper.bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                   phi = list(prior = "pc",  param = c(0.5, 2/3)))
nsim <- 1000

# direct estimates -----
direct.adm1.out <- getDirect(dat_t,ref_year,regionVar = 'admin1.char',timeVar = 'years',clusterVar = '~cluster',Ntrials = 'total')
direct.natl <- direct.adm1.out[direct.adm1.out$region=='All',] %>% select(mean,lower,upper)
direct.adm1 <- direct.adm1.out[direct.adm1.out$region!='All',] %>% rename(admin1.char=region) %>% select(admin1.char,mean,lower,upper)
direct.adm1$median <- direct.adm1$mean
direct.natl$median <- direct.natl$mean

#smoothed direct estimates -----
sd.adm1.fit <- smoothDirect(direct.adm1.out,Amat = admin1.mat,time.model = NULL,year_label=NULL)
sd.adm1 <- getSmoothed(sd.adm1.fit) %>% rename(admin1.char=region) %>% select(admin1.char,median,lower,upper)
sd.adm1$mean <- NA

# mod1: urban (fixed) + admin1 (fixed) -----
mod1 <- INLA::inla(Y ~ factor(strata) + factor(admin1) -1,
                   data=dat_t, family='nbinomial', E=total,
                   control.predictor = list(compute = T, link = 1),
                   control.family = list(link = 'log'),
                   control.compute = list(return.marginals=T, 
                                          return.marginals.predictor=T,config=T))

# sample from full posterior
mod1.est <- getSmoothedFullPost(mod1,dat_t,admin1.weights = adm1_weights_t,admin1.UR.weights = adm1_UR_weights_t)
mod1.adm1 <- mod1.est$adm1.est
mod1.natl <- mod1.est$natl.est

# mod1b: admin1 (fixed) + urban (fixed) + admin2 (BYM2) ----
mod1b <- INLA::inla(Y ~ -1 + factor(strata) + factor(admin1) + 
                      f(admin2, graph = (admin2.mat.nested),model = "bym2",hyper=hyper.bym2, 
                        constr=T, scale.model = T, adjust.for.con.comp = T), #this combination forces sum-to-zero constraints on islands
                    data=dat_t, family='nbinomial', E=total,
                    control.predictor = list(compute = T, link = 1),
                    control.family = list(link = 'log'),
                    control.compute = list(return.marginals=T, 
                                           return.marginals.predictor=T,
                                           config=T))
nsim <- 1000
# sample from full posterior
{
  cs <- mod1b$misc$configs$contents$tag
  cs <- cs[cs != "Predictor"]
  select <- list()
  for (i in 1:length(cs)) {
    select[[i]] <- 0
    names(select)[i] <- cs[i]
  }
  
  sampFull <- INLA::inla.posterior.sample(n = nsim, result = mod1b, intern = TRUE, selection = select)
  sampFull.draws <- matrix(NA,nsim,length(sampFull[[1]]$latent))
  for(i in 1:nsim){
    sampFull.draws[i,] <- sampFull[[i]]$latent
  }
  colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)
  #which columns pertain to which factor
  admin2.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'admin2'))
  admin1.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'(admin1)'))
  strata.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'strata'))
  
  ## get admin2 stratified estimates
  outFull <- expand.grid(urban = unique(dat_t$urban),admin2.char = unique(dat_t$admin2.char))
  outFull <- merge(outFull,admin_key)
  outFull$admin2 <- as.numeric(str_remove(outFull$admin2.char,'admin2_'))
  outFull$admin1 <- as.numeric(str_remove(outFull$admin1.char,'admin1_'))
  outFull <- outFull[order(outFull$admin2),]
  # binary matrix which specifies which columns to index for a given combination of factors
  AA.loc <- matrix(0,nrow(outFull),ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin2)){
    AA.loc[,admin2.cols[k]] <- I(outFull$admin2==k)
  }
  for(k in 2:max(dat_t$admin1)){
    AA.loc[,admin1.cols[k-1]] <- I(outFull$admin1==k)
  }
  AA.loc[,strata.cols[2]] <- I(outFull$urban=='urban')
  AA.loc[,strata.cols[1]] <- 1 - I(outFull$urban=='urban')
  
  AA <- sampFull.draws %*% t(AA.loc)
  
  outFull$log.mean <- colMeans(AA)
  outFull$mean <- colMeans(exp(AA))
  outFull$median <- colMedians(exp(AA))
  outFull$variance <- colVars(exp(AA))
  outFull$lower <- apply(exp(AA),2,quantile,0.05)
  outFull$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get admin 1 stratified estimates
  outFullStrat <- expand.grid(urban = unique(dat_t$urban),admin1.char = unique(dat_t$admin1.char))
  outFullStrat$admin1 <- as.numeric(str_remove(outFullStrat$admin1.char,'admin1_'))
  outFullStrat <- outFullStrat[order(outFullStrat$admin1),]
  # matrix which specifies which columns to include and how to weight them
  AA.wt <- matrix(0,nrow(outFullStrat),ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin1)){
    #include admin1 FE
    if(k>1){
      AA.wt[,admin1.cols[k-1]] <- I(outFullStrat$admin1==k)
    }
    #include weights for admin2s -- really convoluted but it works
    adm2_weights_t <- adm2_to_adm1_weights_t[adm2_to_adm1_weights_t$admin1.char==paste0('admin1_',k),]
    for(m in unique(dat_t$admin2)){
      name_t <- paste0('admin2_',m)
      if(name_t %in% adm2_weights_t$region)
        AA.wt[,admin2.cols[m]] <- adm2_weights_t[adm2_weights_t$region==name_t,]$proportion*I(outFullStrat$admin1==k)
    }
  }
  AA.wt[,strata.cols[2]] <- I(outFullStrat$urban=='urban')
  AA.wt[,strata.cols[1]] <- 1 - I(outFullStrat$urban=='urban')
  
  AA <- sampFull.draws %*% t(AA.wt)
  outFullStrat$log.mean <- colMeans(AA)
  outFullStrat$mean <- colMeans(exp(AA))
  outFullStrat$median <- colMedians(exp(AA))
  outFullStrat$variance <- colVars(exp(AA))
  outFullStrat$lower <- apply(exp(AA),2,quantile,0.05)
  outFullStrat$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get admin 1 overall estimates
  outFullOverall <- data.frame(admin1 = 1:max(dat_t$admin1))
  outFullOverall$admin1.char <- paste0('admin1_', outFullOverall$admin1)
  # matrix which specifies which columns to include and how to weight them
  AA.wt <- matrix(0,nrow(outFullOverall),ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin1)){
    #include admin1 FE
    if(k>1){
      AA.wt[,admin1.cols[k-1]] <- I(outFullOverall$admin1==k)
    }
    #include weights for admin2s -- really convoluted but it works
    adm2_weights_t <- adm2_to_adm1_weights_t[adm2_to_adm1_weights_t$admin1.char==paste0('admin1_',k),]
    for(m in unique(dat_t$admin2)){
      name_t <- paste0('admin2_',m)
      if(name_t %in% adm2_weights_t$region)
        AA.wt[,admin2.cols[m]] <- adm2_weights_t[adm2_weights_t$region==name_t,]$proportion*I(outFullOverall$admin1==k)
    }
  }
  #include weights for urban/rural
  AA.wt[,strata.cols[2]] <- adm1_UR_weights_t$urban
  AA.wt[,strata.cols[1]] <- adm1_UR_weights_t$rural
  
  AA <- sampFull.draws %*% t(AA.wt)
  outFullOverall$log.mean <- colMeans(AA)
  outFullOverall$mean <- colMeans(exp(AA))
  outFullOverall$median <- colMedians(exp(AA))
  outFullOverall$variance <- colVars(exp(AA))
  outFullOverall$lower <- apply(exp(AA),2,quantile,0.05)
  outFullOverall$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get national estimates
  AA.wt <- rep(0,ncol(sampFull.draws))
  for(k in 2:max(dat_t$admin1)){
    #include admin1 FE
    AA.wt[admin1.cols[k-1]] <- adm1_weights_t[adm1_weights_t$region==paste0('admin1_',k),]$proportion
  }
  for(k in 1:max(dat_t$admin2)){
    #include weights for admin2s
    AA.wt[admin2.cols[k]] <- adm2_weights[adm2_weights$years==ref_year & adm2_weights$region==paste0('admin2_',k),]$proportion
  }
  #include weights for urban/rural
  AA.wt[strata.cols[2]] <- sum(adm1_UR_weights_t$urban*adm1_weights_t$proportion)
  AA.wt[strata.cols[1]] <- sum(adm1_UR_weights_t$rural*adm1_weights_t$proportion)
  
  AA <- sampFull.draws %*% AA.wt
  outFullNatl <- data.frame(mean = colMeans(exp(AA)))
  outFullNatl$log.mean <- colMeans(AA)
  outFullNatl$median <- colMedians(exp(AA))
  outFullNatl$variance <- colVars(exp(AA))
  outFullNatl$lower <- apply(exp(AA),2,quantile,0.05)
  outFullNatl$upper <- apply(exp(AA),2,quantile,0.95)
  
}

mod1b.adm2.strat <- outFull
mod1b.adm1.strat <- outFullStrat
mod1b.adm1 <- outFullOverall
mod1b.natl <- outFullNatl

# sample from marginal fixed effects posterior
{
  sampMarg <- matrix(unlist(lapply(mod1b$marginals.fixed,inla.rmarginal,n=nsim)),nrow=nsim,byrow = F)
  colnames(sampMarg) <- names(mod1b$marginals.fixed)
  
  ## get stratified estimates
  outMargStrat <- expand.grid(urban=unique(dat_t$urban), admin1.char = unique(dat_t$admin1.char))
  outMargStrat$admin1 <- as.numeric(str_remove(outMargStrat$admin1.char,'admin1_'))
  outMargStrat <- dat_t %>% select(urban,admin1.char,admin1) %>% unique() %>% arrange(admin1.char,urban)
  outMargStrat <- outMargStrat[order(outMargStrat$admin1),]
  
  AA.loc <- matrix(NA,nrow(outMargStrat),ncol(sampMarg))
  AA.loc[,2] <- I(outMargStrat$urban=='urban')
  AA.loc[,1] <- 1 - I(outMargStrat$urban=='urban')
  for(k in 2:max(dat_t$admin1)){
    AA.loc[,(k+1)] <- I(outMargStrat$admin1==k)
  }
  
  AA <- sampMarg %*% t(AA.loc)
  outMargStrat$log.mean <- colMeans(AA)
  outMargStrat$mean <- colMeans(exp(AA))
  outMargStrat$median <- colMedians(exp(AA))
  outMargStrat$variance <- colVars(exp(AA))
  outMargStrat$lower <- apply(exp(AA),2,quantile,0.05)
  outMargStrat$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get overall estimates
  outMargOverall <- data.frame(admin1 = 1:max(dat_t$admin1))
  outMargOverall$admin1.char <- paste0('admin1_', outMargOverall$admin1)
  
  AA.wt <- matrix(NA,nrow(outMargOverall),ncol(sampMarg))
  AA.wt[,2] <- adm1_UR_weights_t$urban
  AA.wt[,1] <- adm1_UR_weights_t$rural
  for(k in 2:max(dat_t$admin1)){
    AA.wt[,(k+1)] <- I(adm1_UR_weights_t$region==paste0('admin1_',k))
  }
  
  AA <- sampMarg %*% t(AA.wt)
  outMargOverall$log.mean <- colMeans(AA)
  outMargOverall$mean <- colMeans(exp(AA))
  outMargOverall$median <- colMedians(exp(AA))
  outMargOverall$variance <- colVars(exp(AA))
  outMargOverall$lower <- apply(exp(AA),2,quantile,0.05)
  outMargOverall$upper <- apply(exp(AA),2,quantile,0.95)
  
  mod1b.outMargStrat <- outMargStrat
  mod1b.outMargOverall <- outMargOverall
  
}

# mod2: urban*admin1 (fixed) -----
mod2 <- INLA::inla(Y ~ -1 + factor(admin1.strata),
                   data=dat_t, family='nbinomial', E=total,
                    control.predictor = list(compute = T, link = 1),
                    control.family = list(link = 'log'),
                    control.compute = list(config=T))

# sample from full posterior
mod2.est <- getSmoothedFullPost(mod2,dat_t,admin1.weights = adm1_weights_t,admin1.UR.weights = adm1_UR_weights_t)
mod2.adm1 <- mod2.est$adm1.est
mod2.natl <- mod2.est$natl.est

# mod2b: admin1 x urban (fixed) + admin2 (BYM2) ------
mod2b <- INLA::inla(Y ~ -1 + factor(admin1.strata) + 
                      f(admin2, graph = (admin2.mat.nested),model = "bym2",hyper=hyper.bym2, 
                        constr=T, scale.model = T, adjust.for.con.comp = T), #this combination forces sum-to-zero constraints on islands
                    data=dat_t, family='nbinomial', E=total,
                    control.predictor = list(compute = T, link = 1),
                    control.family = list(link = 'log'),
                    control.compute = list(return.marginals=T, 
                                           return.marginals.predictor=T,
                                           config=T))

# sample from full posterior -- estimates are very close to INLA summary output (but very slightly shrunken on average)
{
  cs <- mod2b$misc$configs$contents$tag
  cs <- cs[cs != "Predictor"]
  select <- list()
  for (i in 1:length(cs)) {
    select[[i]] <- 0
    names(select)[i] <- cs[i]
  }
  
  sampFull <- INLA::inla.posterior.sample(n = nsim, result = mod2b, intern = TRUE, selection = select)
  sampFull.draws <- matrix(NA,nsim,length(sampFull[[1]]$latent))
  for(i in 1:nsim){
    sampFull.draws[i,] <- sampFull[[i]]$latent
  }
  colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)
  
  ## get admin2 stratified estimates
  # could generalize by defining all potential variables (interactions)
  outFull <- expand.grid(urban = unique(dat_t$urban),admin2.char = unique(dat_t$admin2.char))
  outFull <- merge(outFull,admin_key)
  outFull$admin2 <- as.numeric(str_remove(outFull$admin2.char,'admin2_'))
  outFull$admin1 <- as.numeric(str_remove(outFull$admin1.char,'admin1_'))
  outFull <- outFull[order(outFull$admin2),]
  outFull$strata <- as.numeric(I(outFull$urban)=='urban')
  outFull$admin1.strata <- outFull$admin1 + (1-outFull$strata)*max(outFull$admin1)
  
  # binary matrix which specifies which columns to index for a given combination of factors
  # just changes based on whatever terms included in model (need to be careful if one factor includes intercept)
  admin2.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'admin2'))
  admin1.strata.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'admin1.strata'))
  AA.loc <- matrix(0,nrow(outFull),ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin2)){
    AA.loc[,admin2.cols[k]] <- I(outFull$admin2==k)
  }
  for(k in 1:max(dat_t$admin1.strata)){
    AA.loc[,admin1.strata.cols[k]] <- I(outFull$admin1.strata==k)
  }
  
  AA <- sampFull.draws %*% t(AA.loc)
  
  outFull$log.mean <- colMeans(AA)
  outFull$mean <- colMeans(exp(AA))
  outFull$median <- colMedians(exp(AA))
  outFull$variance <- colVars(exp(AA))
  outFull$lower <- apply(exp(AA),2,quantile,0.05)
  outFull$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get admin 1 stratified estimates
  outFullStrat <- expand.grid(urban = unique(dat_t$urban),admin1.char = unique(dat_t$admin1.char))
  outFullStrat$admin1 <- as.numeric(str_remove(outFullStrat$admin1.char,'admin1_'))
  outFullStrat <- outFullStrat[order(outFullStrat$admin1),]
  outFullStrat$strata <- as.numeric(I(outFullStrat$urban)=='urban')
  outFullStrat$admin1.strata <- outFullStrat$admin1 + (1-outFullStrat$strata)*max(outFullStrat$admin1)
  
  # matrix which specifies which columns to include and how to weight them
  AA.wt <- matrix(0,nrow(outFullStrat),ncol(sampFull.draws))
  #include admin1 FE
  for(k in 1:max(dat_t$admin1.strata)){
    AA.wt[,admin1.strata.cols[k]] <- I(outFullStrat$admin1.strata==k)
  }
  #include weights for admin2s
  for(k in 1:max(dat_t$admin1)){
    adm2_weights_t <- adm2_to_adm1_weights_t[adm2_to_adm1_weights_t$admin1.char==paste0('admin1_',k),]
    for(m in adm2_weights_t$region){
      AA.wt[,admin2.cols[as.numeric(str_remove(m,'admin2_'))]] <- adm2_weights_t[adm2_weights_t$region==m,]$proportion*I(outFullStrat$admin1==k)
    }
  }
  
  AA <- sampFull.draws %*% t(AA.wt)
  outFullStrat$log.mean <- colMeans(AA)
  outFullStrat$mean <- colMeans(exp(AA))
  outFullStrat$median <- colMedians(exp(AA))
  outFullStrat$variance <- colVars(exp(AA))
  outFullStrat$lower <- apply(exp(AA),2,quantile,0.05)
  outFullStrat$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get admin 1 overall estimates
  outFullOverall <- data.frame(admin1 = 1:max(dat_t$admin1))
  outFullOverall$admin1.char <- paste0('admin1_', outFullOverall$admin1)
  # matrix which specifies which columns to include and how to weight them
  AA.wt <- matrix(0,nrow(outFullOverall),ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin1)){
    #include weights for admin2s
    adm2_weights_t <- adm2_to_adm1_weights_t[adm2_to_adm1_weights_t$admin1.char==paste0('admin1_',k),]
    for(m in unique(dat_t$admin2)){
      name_t <- paste0('admin2_',m)
      if(name_t %in% adm2_weights_t$region)
        AA.wt[,admin2.cols[m]] <- adm2_weights_t[adm2_weights_t$region==name_t,]$proportion*I(outFullOverall$admin1==k)
    }
    #include urban rural weights
    strata_weights_t <- adm1_UR_weights_t[adm1_UR_weights_t$region==paste0('admin1_',k),]
    AA.wt[,admin1.strata.cols[k]] <- strata_weights_t$urban*I(outFullOverall$admin1==k) #urban
    AA.wt[,admin1.strata.cols[k + max(dat_t$admin1)]] <- strata_weights_t$rural*I(outFullOverall$admin1==k) #rural
  }
  
  AA <- sampFull.draws %*% t(AA.wt)
  outFullOverall$log.mean <- colMeans(AA)
  outFullOverall$mean <- colMeans(exp(AA))
  outFullOverall$median <- colMedians(exp(AA))
  outFullOverall$variance <- colVars(exp(AA))
  outFullOverall$lower <- apply(exp(AA),2,quantile,0.05)
  outFullOverall$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get national estimates
  AA.wt <- rep(0,ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin1)){
    #include admin1 x strata FE
    AA.wt[admin1.strata.cols[k]] <- adm1_weights_t[adm1_weights_t$region==paste0('admin1_',k),]$proportion*adm1_UR_weights_t[adm1_weights_t$region==paste0('admin1_',k),]$urban
    AA.wt[admin1.strata.cols[k + max(dat_t$admin1)]] <- adm1_weights_t[adm1_weights_t$region==paste0('admin1_',k),]$proportion*adm1_UR_weights_t[adm1_weights_t$region==paste0('admin1_',k),]$rural
  }
  for(k in 1:max(dat_t$admin2)){
    #include weights for admin2s
    AA.wt[admin2.cols[k]] <- adm2_weights[adm2_weights$years==ref_year & adm2_weights$region==paste0('admin2_',k),]$proportion
  }
  
  AA <- sampFull.draws %*% AA.wt
  outFullNatl <- data.frame(mean = colMeans(exp(AA)))
  outFullNatl$log.mean <- colMeans(AA)
  outFullNatl$median <- colMedians(exp(AA))
  outFullNatl$variance <- colVars(exp(AA))
  outFullNatl$lower <- apply(exp(AA),2,quantile,0.05)
  outFullNatl$upper <- apply(exp(AA),2,quantile,0.95)
}
mod2b.adm2.strat <- outFull
mod2b.adm1.strat <- outFullStrat
mod2b.adm1 <- outFullOverall
mod2b.natl <- outFullNatl

# sample from marginal fixed effects posterior
{
  sampMarg <- matrix(unlist(lapply(mod2b$marginals.fixed,inla.rmarginal,n=nsim)),nrow=nsim,byrow = F)
  colnames(sampMarg) <- names(mod2b$marginals.fixed)
  
  ## get stratified estimates (trivial because each grouping has its own parameter)
  outMargStrat <- expand.grid(urban=unique(dat_t$urban), admin1.char = unique(dat_t$admin1.char))
  outMargStrat$admin1 <- as.numeric(str_remove(outMargStrat$admin1.char,'admin1_'))
  outMargStrat <- dat_t %>% select(urban,admin1.char,admin1) %>% unique() %>% arrange(admin1.char,urban)
  outMargStrat <- outMargStrat[order(outMargStrat$admin1),]
  outMargStrat$strata <- as.numeric(I(outMargStrat$urban=='urban'))
  outMargStrat$admin1.strata <- outMargStrat$admin1 + (1-outMargStrat$strata)*max(outMargStrat$admin1)
  
  AA <- sampMarg
  outMargStrat$log.mean <- colMeans(AA)
  outMargStrat$mean <- colMeans(exp(AA))
  outMargStrat$median <- colMedians(exp(AA))
  outMargStrat$variance <- colVars(exp(AA))
  outMargStrat$lower <- apply(exp(AA),2,quantile,0.05)
  outMargStrat$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get overall estimates
  outMargOverall <- data.frame(admin1 = 1:max(dat_t$admin1))
  outMargOverall$admin1.char <- paste0('admin1_', outMargOverall$admin1)
  
  AA.wt <- matrix(NA,nrow(outMargOverall),ncol(sampMarg))
  for(k in 1:max(dat_t$admin1)){
    weights_t <- adm1_UR_weights_t[adm1_UR_weights_t$region==paste0('admin1_',k),]
    AA.wt[,k] <- weights_t$urban*I(outMargOverall$admin1==k) # urban weight
    AA.wt[,(k + max(dat_t$admin1))] <- weights_t$rural*I(outMargOverall$admin1==k) # rural weight
  }
  
  AA <- sampMarg %*% t(AA.wt)
  outMargOverall$log.mean <- colMeans(AA)
  outMargOverall$mean <- colMeans(exp(AA))
  outMargOverall$median <- colMedians(exp(AA))
  outMargOverall$variance <- colVars(exp(AA))
  outMargOverall$lower <- apply(exp(AA),2,quantile,0.05)
  outMargOverall$upper <- apply(exp(AA),2,quantile,0.95)
  
  mod2b.outMargStrat <- outMargStrat
  mod2b.outMargOverall <- outMargOverall
}

# mod3: urban (fixed) + admin2 (BYM2) ---------
mod3 <- INLA::inla(Y ~ -1 + factor(strata) +
                     f(admin2, graph = (admin2.mat),model = "bym2",hyper=hyper.bym2, 
                       constr=T, scale.model = T, adjust.for.con.comp = T), #this combination forces sum-to-zero constraints on islands
                   data=dat_t, family='nbinomial', E=total,
                   control.predictor = list(compute = T, link = 1),
                   control.family = list(link = 'log'),
                   control.compute = list(return.marginals=T, 
                                          return.marginals.predictor=T,
                                          config=T))

# sample from full posterior
{
  cs <- mod3$misc$configs$contents$tag
  cs <- cs[cs != "Predictor"]
  select <- list()
  for (i in 1:length(cs)) {
    select[[i]] <- 0
    names(select)[i] <- cs[i]
  }
  
  sampFull <- INLA::inla.posterior.sample(n = nsim, result = mod3, intern = TRUE, selection = select)
  sampFull.draws <- matrix(NA,nsim,length(sampFull[[1]]$latent))
  for(i in 1:nsim){
    sampFull.draws[i,] <- sampFull[[i]]$latent
  }
  colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)
  #which columns pertain to which factor
  admin2.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'admin2'))
  strata.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'strata'))
  
  ## get admin2 stratified estimates
  outFull <- expand.grid(urban = unique(dat_t$urban),admin2.char = unique(dat_t$admin2.char))
  outFull$admin2 <- as.numeric(str_remove(outFull$admin2.char,'admin2_'))
  outFull <- outFull[order(outFull$admin2),]
  # binary matrix which specifies which columns to index for a given combination of factors
  AA.loc <- matrix(0,nrow(outFull),ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin2)){
    AA.loc[,admin2.cols[k]] <- I(outFull$admin2==k)
  }
  AA.loc[,strata.cols[2]] <- I(outFull$urban=='urban')
  AA.loc[,strata.cols[1]] <- 1 - I(outFull$urban=='urban')
  
  AA <- sampFull.draws %*% t(AA.loc)
  
  outFull$log.mean <- colMeans(AA)
  outFull$mean <- colMeans(exp(AA))
  outFull$median <- colMedians(exp(AA))
  outFull$variance <- colVars(exp(AA))
  outFull$lower <- apply(exp(AA),2,quantile,0.05)
  outFull$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get admin 1 stratified estimates
  outFullStrat <- expand.grid(urban = unique(dat_t$urban),admin1.char = unique(dat_t$admin1.char))
  outFullStrat$admin1 <- as.numeric(str_remove(outFullStrat$admin1.char,'admin1_'))
  outFullStrat <- outFullStrat[order(outFullStrat$admin1),]
  # matrix which specifies which columns to include and how to weight them
  AA.wt <- matrix(0,nrow(outFullStrat),ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin1)){
    #include weights for admin2s
    adm2_weights_t <- adm2_to_adm1_weights_t[adm2_to_adm1_weights_t$admin1.char==paste0('admin1_',k),]
    for(m in adm2_weights_t$region){
      AA.wt[,admin2.cols[as.numeric(str_remove(m,'admin2_'))]] <- adm2_weights_t[adm2_weights_t$region==m,]$proportion*I(outFullStrat$admin1==k)
    }
  }
  AA.wt[,strata.cols[2]] <- I(outFullStrat$urban=='urban')
  AA.wt[,strata.cols[1]] <- 1 - I(outFullStrat$urban=='urban')
  
  AA <- sampFull.draws %*% t(AA.wt)
  outFullStrat$log.mean <- colMeans(AA)
  outFullStrat$mean <- colMeans(exp(AA))
  outFullStrat$median <- colMedians(exp(AA))
  outFullStrat$variance <- colVars(exp(AA))
  outFullStrat$lower <- apply(exp(AA),2,quantile,0.05)
  outFullStrat$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get admin 1 overall estimates
  outFullOverall <- data.frame(admin1 = 1:max(dat_t$admin1))
  outFullOverall$admin1.char <- paste0('admin1_', outFullOverall$admin1)
  # matrix which specifies which columns to include and how to weight them
  AA.wt <- matrix(0,nrow(outFullOverall),ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin1)){
    #include weights for admin2s
    adm2_weights_t <- adm2_to_adm1_weights_t[adm2_to_adm1_weights_t$admin1.char==paste0('admin1_',k),]
    for(m in adm2_weights_t$region){
      AA.wt[,admin2.cols[as.numeric(str_remove(m,'admin2_'))]] <- adm2_weights_t[adm2_weights_t$region==m,]$proportion*I(outFullOverall$admin1==k)
    }
  }
  #include weights for urban/rural
  AA.wt[,strata.cols[2]] <- adm1_UR_weights_t$urban
  AA.wt[,strata.cols[1]] <- adm1_UR_weights_t$rural
  
  AA <- sampFull.draws %*% t(AA.wt)
  outFullOverall$log.mean <- colMeans(AA)
  outFullOverall$mean <- colMeans(exp(AA))
  outFullOverall$median <- colMedians(exp(AA))
  outFullOverall$variance <- colVars(exp(AA))
  outFullOverall$lower <- apply(exp(AA),2,quantile,0.05)
  outFullOverall$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get national estimates
  AA.wt <- rep(0,ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin2)){
    #include weights for admin2s
    AA.wt[admin2.cols[k]] <- adm2_weights[adm2_weights$years==ref_year & adm2_weights$region==paste0('admin2_',k),]$proportion
  }
  #include weights for urban/rural
  AA.wt[strata.cols[2]] <- sum(adm1_UR_weights_t$urban*adm1_weights_t$proportion)
  AA.wt[strata.cols[1]] <- sum(adm1_UR_weights_t$rural*adm1_weights_t$proportion)
  
  AA <- sampFull.draws %*% AA.wt
  outFullNatl <- data.frame(mean = colMeans(exp(AA)))
  outFullNatl$log.mean <- colMeans(AA)
  outFullNatl$median <- colMedians(exp(AA))
  outFullNatl$variance <- colVars(exp(AA))
  outFullNatl$lower <- apply(exp(AA),2,quantile,0.05)
  outFullNatl$upper <- apply(exp(AA),2,quantile,0.95)
  
}
mod3.adm2.strat <- outFull
mod3.adm1.strat <- outFullStrat
mod3.adm1 <- outFullOverall
mod3.natl <- outFullNatl

# don't think there is a logical way to sample marginal without an admin1 fixed effect

# mod4: urban(fixed) + admin2 (fixed) -----
mod4 <- INLA::inla(Y ~ factor(strata) + factor(admin2) -1,
                   data=dat_t, family='nbinomial', E=total,
                   control.predictor = list(compute = T, link = 1),
                   control.family = list(link = 'log'),
                   control.compute = list(config=T))

# sample from full posterior
{
  cs <- mod4$misc$configs$contents$tag
  cs <- cs[cs != "Predictor"]
  select <- list()
  for (i in 1:length(cs)) {
    select[[i]] <- 0
    names(select)[i] <- cs[i]
  }
  
  sampFull <- INLA::inla.posterior.sample(n = nsim, result = mod4, intern = TRUE, selection = select)
  sampFull.draws <- matrix(NA,nsim,length(sampFull[[1]]$latent))
  for(i in 1:nsim){
    sampFull.draws[i,] <- sampFull[[i]]$latent
  }
  colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)
  #which columns pertain to which factor
  admin2.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'admin2'))
  strata.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'strata'))
  
  ## get admin 1 overall estimates
  outFullOverall <- data.frame(admin1 = 1:max(dat_t$admin1))
  outFullOverall$admin1.char <- paste0('admin1_', outFullOverall$admin1)
  # matrix which specifies which columns to include and how to weight them
  AA.wt <- matrix(0,nrow(outFullOverall),ncol(sampFull.draws))
  for(k in 1:max(dat_t$admin1)){
    #include weights for admin2s
    adm2_weights_t <- adm2_to_adm1_weights_t[adm2_to_adm1_weights_t$admin1.char==paste0('admin1_',k),]
    for(m in adm2_weights_t$region){
      if(m!='admin2_1'){
        AA.wt[,admin2.cols[as.numeric(str_remove(m,'admin2_'))-1]] <- adm2_weights_t[adm2_weights_t$region==m,]$proportion*I(outFullOverall$admin1==k)
      }
    }
  }
  #include weights for urban/rural
  AA.wt[,strata.cols[2]] <- adm1_UR_weights_t$urban
  AA.wt[,strata.cols[1]] <- adm1_UR_weights_t$rural
  
  AA <- sampFull.draws %*% t(AA.wt)
  outFullOverall$log.mean <- colMeans(AA)
  outFullOverall$mean <- colMeans(exp(AA))
  outFullOverall$median <- colMedians(exp(AA))
  outFullOverall$variance <- colVars(exp(AA))
  outFullOverall$lower <- apply(exp(AA),2,quantile,0.05)
  outFullOverall$upper <- apply(exp(AA),2,quantile,0.95)
  
  ## get national estimates
  AA.wt <- rep(0,ncol(sampFull.draws))
  for(k in 2:max(dat_t$admin2)){
    #include weights for admin2s
    AA.wt[admin2.cols[k-1]] <- adm2_weights[adm2_weights$years==ref_year & adm2_weights$region==paste0('admin2_',k),]$proportion
  }
  #include weights for urban/rural
  AA.wt[strata.cols[2]] <- sum(adm1_UR_weights_t$urban*adm1_weights_t$proportion)
  AA.wt[strata.cols[1]] <- sum(adm1_UR_weights_t$rural*adm1_weights_t$proportion)
  
  AA <- sampFull.draws %*% AA.wt
  outFullNatl <- data.frame(mean = colMeans(exp(AA)))
  outFullNatl$log.mean <- colMeans(AA)
  outFullNatl$median <- colMedians(exp(AA))
  outFullNatl$variance <- colVars(exp(AA))
  outFullNatl$lower <- apply(exp(AA),2,quantile,0.05)
  outFullNatl$upper <- apply(exp(AA),2,quantile,0.95)
  
}

mod4.adm1 <- outFullOverall
mod4.natl <- outFullNatl

# compare non-temporal models ----
direct.adm1$mod <- direct.natl$mod <- 'Direct'
sd.adm1$mod <- 'SD'
mod1.adm1.strat$mod <- mod1.adm1$mod <- mod1.natl$mod <- 'Urban + Admin1'
mod1b.adm1.strat$mod <- mod1b.adm1$mod <- mod1b.natl$mod <- 'Urban + Admin1 + Admin2 (nested BYM2)'
mod2.adm1$mod <- mod2.natl$mod <- 'Urban*Admin1'
mod2b.adm1.strat$mod <- mod2b.adm1$mod <- mod2b.natl$mod <- 'Urban*Admin1 + Admin2 (nested BYM2)'
mod3.adm2.strat$mod <- mod3.adm1.strat$mod <- mod3.adm1$mod <- mod3.natl$mod <- 'Urban + Admin2 (BYM2)'
mod4.adm1$mod <- mod4.natl$mod <- 'Urban + Admin2 (fixed)'

# compare Admin1 estimates
bind_vars <- c('admin1.char','mean','median','lower','upper','mod')
adm1.ests_t <- rbind(direct.adm1[,bind_vars],
                    #sd.adm1[,bind_vars],
                    mod1.adm1[,bind_vars],
                    mod1b.adm1[,bind_vars],
                    mod2.adm1[,bind_vars],
                    mod2b.adm1[,bind_vars],
                    mod3.adm1[,bind_vars],
                    mod4.adm1[,bind_vars])
adm1.ests_t <- adm1.ests_t %>% mutate(level = ifelse(mod %in% c('Direct','Urban + Admin1'),'admin1','admin2'))

# adm1.ests_t %>% ggplot(aes(x=admin1.char,y=mean,group=mod,color=mod)) + geom_point(aes(pch=level),size=3) + 
#   scale_color_discrete(name='Model') + scale_shape_discrete(name='Level of model',solid=F)

# compare National estimates
bind_vars <- c('mean','median','lower','upper','mod')
natl.ests_t <- rbind(direct.natl[,bind_vars],
                     mod1.natl[,bind_vars],
                     mod1b.natl[,bind_vars],
                     mod2.natl[,bind_vars],
                     mod2b.natl[,bind_vars],
                     mod3.natl[,bind_vars],
                     mod4.natl[,bind_vars])
natl.ests_t <- natl.ests_t %>% mutate(level = ifelse(mod =='Direct','national', ifelse(mod =='Urban + Admin1','admin1','admin2')))
natl.ests_t$admin1.char <- 'National'
# natl.ests_t %>% ggplot(aes(x=NA,y=mean,group=mod,color=mod)) + geom_point(aes(pch=level),size=3) + 
#   scale_color_discrete(name='Model') + scale_shape_discrete(name='Level of model',solid=F)

all.ests_t <- rbind(adm1.ests_t,natl.ests_t)
all.ests_t$period <- paste0(beg.years[period],'-',end.years[period])
all.ests_t$country <- country_t

all.ests <- rbind(all.ests,all.ests_t)
}


pdf("/Users/alanamcgovern/Desktop/Research/New Benchmarking/Nested Model Approach/Guinea Model Comparison.pdf")

g <- all.ests %>% ggplot(aes(x=admin1.char,y=median,group=mod,color=mod)) + 
  geom_jitter(aes(pch=mod,alpha=mod),size=3,width=0.2) +
  facet_wrap(~period) + scale_color_discrete(name='Model') + 
  scale_shape_manual(name='Model', values = c(8,8,17,17,17,8,17)) + 
  scale_alpha_manual(name = 'Model',values = c(1,1,0.6,0.6,0.6,1,0.6)) +
  ggtitle('Guinea, 2018 survey') + xlab('Admin 1') + ylab('Median NMR') +
  theme(legend.position = 'bottom',legend.box = 'vertical',legend.margin = margin(),
        axis.text.x = element_text(angle = 45,vjust = 0.75)) + 
  guides(color=guide_legend(nrow=3,byrow=T,title.position = 'top',title.hjust = 0.5))
print(g)

g <- all.ests %>% filter(!(mod %in% c('Urban*Admin1','Urban*Admin1 + Admin2 (nested BYM2)','Urban + Admin2 (BYM2)'))) %>% 
  ggplot(aes(x=admin1.char,y=median,group=mod,color=mod)) + 
  geom_jitter(aes(pch=mod,alpha=mod),size=3,width=0.2) +
  facet_wrap(~period) + scale_color_discrete(name='Model') + 
  scale_shape_manual(name='Model', values = c(8,8,17,17,17,8,17)) + 
  scale_alpha_manual(name = 'Model',values = c(1,1,0.6,0.6,0.6,1,0.6)) +
  ggtitle('Guinea, 2018 survey') + xlab('Admin 1') + ylab('Median NMR') +
  theme(legend.position = 'bottom',legend.box = 'vertical',legend.margin = margin(),
        axis.text.x = element_text(angle = 45,vjust = 0.75)) + 
  guides(color=guide_legend(nrow=3,byrow=T,title.position = 'top',title.hjust = 0.5))
print(g)

dev.off()


# fit models with time component -----
hyper.rw2 <-  list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
hyper.ar1 = list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                 theta2 = list(prior = "pc.cor1",param = c(0.7, 0.9)))

plot_list <- list()
for(survey_t in c(2010,2014,2015,2020)){
  #data
  dat_t <- dat[dat$survey==survey_t,]

  # admin 1 level models ---------
  #direct estimates and smoothed ------
adm1_dir_est <- getDirect(births = dat_t, years = unique(dat_t$years), regionVar = 'admin1.char',
                          timeVar = 'years', clusterVar = '~cluster', Ntrials = 'total')

natl_dir_est <- adm1_dir_est[adm1_dir_est$region=='All',]
adm1_dir_est <- adm1_dir_est[adm1_dir_est$region!='All',]

adm1_sd_fit <- smoothDirect(adm1_dir_est, Amat = admin1.mat, time.model = 'rw2',
                            year_label = as.character(2000:2015), type.st = 1,
                            year_range = 2000:2015, is.yearly = F)

adm1_sd_est <- getSmoothed(adm1_sd_fit, Amat = admin1.mat, CI=0.95,
                           year_label = as.character(2000:2015),
                           year_range = 2000:2015, control.compute=list(cpo=T))

adm1_sd_est$mean <- (adm1_sd_est$upper - adm1_sd_est$lower)/2 + adm1_sd_est$lower

  #mod1: just admin1 (fixed effect) -- may not have enough data if only 1 year of data is included/many admin1 areas -----
  mod1 <- INLA::inla(Y ~ urban + factor(admin1) +
                            f(years.int, model='rw2', hyper = hyper.rw2, scale.model = T,
                              constr = T, extraconstr = NULL, values = 1:max(years.int)) -1,
                        data=dat_t, family='nbinomial', E=total,
                        control.predictor = list(link = 1),
                        control.compute = list(cpo=T))

  mod1_adm1_fe <- rbind(rep(0,3),mod1$summary.fixed[3:nrow(mod1$summary.fixed),c(1,3,5)])
  mod1_strata_fe <- mod1$summary.fixed[1:2,c(1,3,5)]
  mod1_time_re <- mod1$summary.random$years.int[,c(1,2,4,6)]

  mod1_summary <- expand.grid(admin1.char = unique(dat_t$admin1.char),years = unique(dat_t$years),strata = unique(dat_t$urban))
  mod1_summary <- mod1_summary %>% mutate(years.int = years - min(mod1_summary$years) +1,
                                        strata.num = as.numeric(strata),
                                        admin1 = as.numeric(str_remove(admin1.char,'admin1_')),
                                        model = 'Adm1 FE')
  mod1_summary$logit.mean <- mod1_adm1_fe$mean[mod1_summary$admin1] + mod1_strata_fe$mean[mod1_summary$strata.num] + mod1_time_re$mean[mod1_summary$years.int]
  mod1_summary$mean <- expit(mod1_summary$logit.mean)
  mod1_summary  <- mod1_summary %>% select(model, admin1.char, years,strata, mean) %>% rename(region=admin1.char)

  # organize admin1 results ------

 adm1_dir_est <- data.frame(model='Direct',region = adm1_dir_est$region, years = adm1_dir_est$years,
                            mean = adm1_dir_est$mean)
 adm1_sd_est <- data.frame(model='SD',region = adm1_sd_est$region, years = adm1_sd_est$years.num,
                           mean = adm1_sd_est$mean)

  inla_adm1_ests <- mod1_summary
  inla_adm1_ests <- merge(adm1_UR_weights,inla_adm1_ests,by=c('region','strata','years'))
  inla_adm1_ests <- inla_adm1_ests %>% group_by(model,region,years) %>% summarise(mean =(sum(mean*proportion))) 
                                                                          

  # mod3: both admin1 (fixed effect) and admin2 (nested) -----
  mod3<- INLA::inla(Y ~ urban + factor(admin1) -1 +
                            f(admin2, graph = (admin2.mat.nested),model = "bym2",hyper=hyper.bym2, 
                              constr=T, scale.model = T, adjust.for.con.comp = T) + #this combination forces sum-to-zero constraints on islands
                            f(years.int, model='rw2', hyper = hyper.rw2, scale.model = T,
                              constr = T, extraconstr = NULL, values = 1:max(years.int)),
                          data=dat_t, family='nbinomial', E=total,
                          control.predictor = list(compute = F, link = 1),
                          control.compute = list(cpo=T))

  mod3_adm1_fe <- rbind(rep(0,3),mod3$summary.fixed[3:nrow(mod3$summary.fixed),c(1,3,5)])
  mod3_strata_fe <- mod3$summary.fixed[1:2,c(1,3,5)]
  mod3_adm2_re <- mod3$summary.random$admin2[,c(1,2,4,6)]
  mod3_time_re <- mod3$summary.random$years.int[,c(1,2,4,6)]

  mod3_summary <- expand.grid(admin2.char = unique(dat_t$admin2.char),years = unique(dat_t$years),strata = unique(dat_t$urban))
  mod3_summary <- merge(mod3_summary,admin_key,by='admin2.char')
  mod3_summary <- mod3_summary %>% mutate(years.int = years - min(mod1_summary$years) +1,
                                        strata.num = as.numeric(strata),
                                        admin1 = as.numeric(str_remove(admin1.char,'admin1_')),
                                        admin2 = as.numeric(str_remove(admin2.char,'admin2_')),
                                        model = 'Adm1 FE + Adm2')

  mod3_summary$logit.mean <- mod3_adm1_fe$mean[mod3_summary$admin1] + mod3_adm2_re$mean[mod3_summary$admin2] + mod3_strata_fe$mean[mod3_summary$strata.num] + mod3_time_re$mean[mod3_summary$years.int]
  mod3_summary$mean <- expit(mod3_summary$logit.mean)
  mod3_summary  <- mod3_summary %>% select(model, admin2.char, admin1.char, years,strata, mean) %>% rename(region=admin2.char)
  
  # mod4: just admin2 -----
  mod4<- INLA::inla(Y ~ urban -1 +
                      f(admin2, graph = (admin2.mat),model = "bym2",hyper=hyper.bym2, 
                        constr=T, scale.model = T, adjust.for.con.comp = T) + 
                      f(years.int, model='rw2', hyper = hyper.rw2, scale.model = T,
                        constr = T, extraconstr = NULL, values = 1:max(years.int)),
                    data=dat_t, family='nbinomial', E=total,
                    control.predictor = list(compute = F, link = 1),
                    control.compute = list(cpo=T))
  
  mod4_strata_fe <- mod4$summary.fixed[1:2,c(1,3,5)]
  mod4_adm2_re <- mod4$summary.random$admin2[,c(1,2,4,6)]
  mod4_time_re <- mod4$summary.random$years.int[,c(1,2,4,6)]
  
  mod4_summary <- expand.grid(admin2.char = unique(dat_t$admin2.char),years = unique(dat_t$years),strata = unique(dat_t$urban))
  mod4_summary <- mod4_summary %>% mutate(years.int = years - min(mod1_summary$years) +1,
                                          strata.num = as.numeric(strata),
                                          admin2 = as.numeric(str_remove(admin2.char,'admin2_')),
                                          model = 'Adm2')
  
  mod4_summary$logit.mean <- mod4_adm2_re$mean[mod4_summary$admin2] + mod4_strata_fe$mean[mod4_summary$strata.num] + mod4_time_re$mean[mod4_summary$years.int]
  mod4_summary$mean <- expit(mod4_summary$logit.mean)
  mod4_summary <- merge(mod4_summary,admin_key,by='admin2.char')
  
  mod4_summary  <- mod4_summary %>% select(model, admin2.char, admin1.char, years,strata, mean) %>% rename(region=admin2.char)

  # organize admin 2 results ----

  inla_adm2_ests <- rbind(mod3_summary,mod4_summary)
  #merge from U/R to admin2
  inla_adm2_ests <- merge(adm2_UR_weights,inla_adm2_ests,by=c('region','strata','years'))
  inla_adm2_ests <- inla_adm2_ests %>% group_by(model,admin1.char,region,years) %>% summarise(mean = (sum(mean*proportion)))
  #merge from admin2 to admin1
  inla_adm2_ests <- merge(adm2_to_adm1_weights,inla_adm2_ests,by=c('region','admin1.char','years'))
  inla_adm2_to_adm1_ests <- inla_adm2_ests %>% group_by(model,admin1.char,years) %>% summarise(mean = (sum(mean*proportion))) %>% rename(region=admin1.char)

  # combine results for admin ------
  adm1_ests_inla_all <- rbind(adm1_dir_est,
                            adm1_sd_est,
                            inla_adm1_ests,inla_adm2_to_adm1_ests)
  
  plot_list[[which(survey_t==c(2010,2014,2015,2020))]] <- adm1_ests_inla_all %>% ggplot() + geom_line(aes(x=years,y=mean,color=region,lty=model)) +ggtitle(paste0('Malawi, ',survey_t,' Survey'))
}

pdf("/Users/alanamcgovern/Desktop/Research/New Benchmarking/Nested Model Approach/Malawi Agg Admin1 Estimates.pdf")
  ggarrange(plotlist = plot_list, common.legend = T)
dev.off()

# aggregate to national ------
adm1_ests <- merge(adm1_ests_inla_all,weight.adm1.u1,by=c('region','years'))
natl_ests <- cbind(adm1_ests %>% group_by(model,years) %>% summarise(mean=sum(mean*proportion),lower95=sum(lower95*proportion), upper95=sum(upper95*proportion)))

natl_ests %>% ggplot() + geom_line(aes(x=years,y=mean,lty=model))
# compare CV diagnostics --------
## first check that they are reliable
sum(inla_adm1fe$cpo$failure + inla_adm1re$cpo$failure + 
      inla_bothfe$cpo$failure + inla_bothre$cpo$failure + inla_adm2re$cpo$failure)==0
sum(inla_adm1re$cpo$failure) + 
      sum(inla_bothfe$cpo$failure) + sum(inla_bothre$cpo$failure) + sum(inla_adm2re$cpo$failure)==0

# very similar, but are in agreement with what we see from comparing aggregated national estimates
cpo_results <- data.frame(model=c('Adm1 FE','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2'),
           CPO=c(-mean(log(inla_adm1fe$cpo$cpo)), -mean(log(inla_adm1re$cpo$cpo)),
                 -mean(log(inla_bothfe$cpo$cpo)), -mean(log(inla_bothre$cpo$cpo)),
                 -mean(log(inla_adm2re$cpo$cpo))))

cpo_results <- data.frame(model=c('Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2'),
                          CPO=c( -mean(log(inla_adm1re$cpo$cpo)),
                                -mean(log(inla_bothfe$cpo$cpo)), -mean(log(inla_bothre$cpo$cpo)),
                                -mean(log(inla_adm2re$cpo$cpo))))

pdf('Nested Model Approach/Malawi 2010 survey=2015 230814.pdf')
{
  plot(seq(-4,-2.8,.01),dnorm(seq(-4,-2.8,.01),natl_dir_est$logit.est,sqrt(1/natl_dir_est$logit.prec)),type='l',main='National Estimates')
  for(k in 1:nrow(natl_ests)){
    abline(v=logit(natl_ests$mean[k]),col=k+1)
  }
  legend('topright',legend=c('SD','Adm1 FE','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2'),
  #legend('topright',legend=c('SD','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2'),
         col = 2:(nrow(natl_ests)+1),lwd=1.5)
  
  g1 <- ggplot() + geom_point(data=adm1_ests_inla_all, aes(x=region,y=logit(mean),col=model),pch=3) +
    geom_point(data=adm1_dir_est,aes(x=region,y=logit.est),pch=2) +
    scale_color_discrete(labels=c('Adm1 FE','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2','SD')) +
    #scale_color_discrete(labels=c('Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2','SD')) +
    ggtitle('Admin1 Estimates')
  print(g1)
  
  #plot distribution of direct estimates
  for(area in 1:nrow(admin1.mat)){
    x_seq <- seq(adm1_dir_est$logit.est[area] - 3*sqrt(adm1_dir_est$logit.var.est[area]),
                 adm1_dir_est$logit.est[area] + 3*sqrt(adm1_dir_est$logit.var.est[area]), 0.01)
    plot(x_seq,dnorm(x_seq,adm1_dir_est$logit.est[area],sqrt(adm1_dir_est$logit.var.est[area])),type='l',
         main=paste0('Distribution of direct estimate for', adm1_dir_est$region[area]),ylab='',xlab='')
    for(k in 1:length(unique(adm1_ests_inla_all$model))){
      dat_tmp <- adm1_ests_inla_all[adm1_ests_inla_all$model == unique(adm1_ests_inla_all$model)[k],]
      abline(v=logit(dat_tmp[area,]$mean),col=k+1)
    }
    abline(v=natl_dir_est$logit.est,lty=2)
    legend('topright',legend=c('SD','Adm1 FE','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2',"Natl direct"),
    #legend('topright',legend=c('SD','Adm1','Adm1 FE + Adm2','Adm1 + Adm2','Adm2',"Natl direct"),
           col = c(2:(length(unique(adm1_ests_inla_all$model))+1),1),lwd=1.5,
           lty=c(rep(1,length(unique(adm1_ests_inla_all$model))),2))
    Sys.sleep(0.5)
  }

}
dev.off()
