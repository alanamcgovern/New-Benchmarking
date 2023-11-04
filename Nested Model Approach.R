library(INLA)
library(SUMMER)
library(tidyverse)
library(ggpubr)
library(Rfast)

# load data ------
load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Malawi/Malawi_cluster_dat.rda")
dat <- mod.dat[mod.dat$age==0,]
dat$died <- dat$Y
dat$years.int <- as.integer(dat$years)
dat$years <- as.numeric(as.character(dat$years))

dat$v005 <- dat$v005/1e6
admin_key <- dat %>% dplyr::select(admin1.char,admin2.char) %>% unique()
admin_key_num <- dat %>% dplyr::select(admin1,admin2) %>% unique()

#adjacency matrix
load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Malawi/shapeFiles_gadm/Malawi_Amat.rda")
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
load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Malawi/shapeFiles_gadm/Malawi_Amat_names.rda")

# get weights for aggregating from admin2 to admin 1
load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Malawi/worldpop/adm1_weights_u1.rda")
adm1_weights <- weight.adm1.u1
load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Malawi/worldpop/adm2_weights_u1.rda")
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
adm1_UR_weights <- readRDS("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Results/Malawi/UR/U1_fraction/admin1_u1_urban_weights.rds")
adm1_UR_weights <- gather(adm1_UR_weights,strata,proportion,urban:rural)
adm2_UR_weights <- readRDS("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Results/Malawi/UR/U1_fraction/admin2_u1_urban_weights.rds")
adm2_UR_weights <- gather(adm2_UR_weights,strata,proportion,urban:rural)

# INLA hyperpriors ------
hyper.rw2 <-  list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
hyper.ar1 = list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                 theta2 = list(prior = "pc.cor1",param = c(0.7, 0.9)))
hyper.bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                   phi = list(prior = "pc",  param = c(0.5, 2/3)))

# fit models as if 1 year of data ----
dat_t <- dat[dat$survey==2015 & dat$years%in%2009:2012,]
adm1_UR_weights_t <- adm1_UR_weights[adm1_UR_weights$years==2011,]
adm1_UR_weights_t <- spread(adm1_UR_weights_t,strata,proportion) %>% select(-years)
adm2_to_adm1_weights_t <- adm2_to_adm1_weights[adm2_to_adm1_weights$years==2011,]

#mod1: just admin1 (fixed effect) -- may not have enough data if only 1 year of data is included/many admin1 areas
mod1 <- INLA::inla(Y ~ urban + factor(admin1) -1,
                   data=dat_t, family='nbinomial', E=total,
                   control.predictor = list(compute = T, link = 1),
                   control.family = list(link = 'log'),
                   control.compute = list(return.marginals=T, 
                                          return.marginals.predictor=T,config=T))

# sample from whole joint posterior and from marginal fixed effects posterior gives same thing (which it should because this is just a fixed effects model)


mod3 <- INLA::inla(Y ~ -1 + urban + factor(admin1) + 
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
  admin1.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'admin1'))
  strata.cols <- which(str_detect(row.names(sampFull[[1]]$latent),'urban'))
  
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
  AA.loc[,strata.cols[1]] <- I(outFull$urban=='urban')
  AA.loc[,strata.cols[2]] <- 1 - I(outFull$urban=='urban')
  
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
  AA.wt[,strata.cols[1]] <- I(outFullStrat$urban=='urban')
  AA.wt[,strata.cols[2]] <- 1 - I(outFullStrat$urban=='urban')

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
  AA.wt[,strata.cols[1]] <- adm1_UR_weights_t$urban
  AA.wt[,strata.cols[2]] <- adm1_UR_weights_t$rural
  
  AA <- sampFull.draws %*% t(AA.wt)
  outFullOverall$log.mean <- colMeans(AA)
  outFullOverall$mean <- colMeans(exp(AA))
  outFullOverall$median <- colMedians(exp(AA))
  outFullOverall$variance <- colVars(exp(AA))
  outFullOverall$lower <- apply(exp(AA),2,quantile,0.05)
  outFullOverall$upper <- apply(exp(AA),2,quantile,0.95)
}

# sample from marginal fixed effects posterior
{
  sampMarg <- matrix(unlist(lapply(mod3$marginals.fixed,inla.rmarginal,n=nsim)),nrow=nsim,byrow = F)
  colnames(sampMarg) <- names(mod3$marginals.fixed)
  
  ## get stratified estimates
  outMargStrat <- expand.grid(urban=unique(dat_t$urban), admin1.char = unique(dat_t$admin1.char))
  outMargStrat$admin1 <- as.numeric(str_remove(outMargStrat$admin1.char,'admin1_'))
  outMargStrat <- dat_t %>% select(urban,admin1.char,admin1) %>% unique() %>% arrange(admin1.char,urban)
  outMargStrat <- outMargStrat[order(outMargStrat$admin1),]
  
  AA.loc <- matrix(NA,nrow(outMargStrat),ncol(sampMarg))
  AA.loc[,1] <- I(outMargStrat$urban=='urban')
  AA.loc[,2] <- 1 - I(outMargStrat$urban=='urban')
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
  AA.wt[,1] <- adm1_UR_weights_t$urban
  AA.wt[,2] <- adm1_UR_weights_t$rural
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
  
}

# compare outMargOverall and outFullOverall


# fit models with time component -----
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
