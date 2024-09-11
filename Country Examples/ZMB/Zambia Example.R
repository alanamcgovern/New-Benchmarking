library(sf)
library(surveyPrev)
library(SUMMER)
library(tidyverse)
library(data.table)
library(LaplacesDemon)
library(VGAM)
library(Matrix)
library(Rfast)
library(MASS)
library(locfit)
library(openxlsx)
library(ggpubr)
library(cmdstanr)
library(surveyPrev)
library(scales)
library(ggrepel)

options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

load(file = paste0('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Info/Zambia_general_info.Rdata')) # load the country info

# load data ----
# adjacency matrix
poly.path <- paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_ZMB_shp")

poly.adm1 <- st_read(dsn = poly.path, layer = poly.layer.adm1, options = "ENCODING=UTF-8")
poly.adm1$admin1 <- 1:nrow(poly.adm1)
poly.adm1$admin1.char <- paste0('admin1_',1:nrow(poly.adm1))

poly.adm2 <- st_read(dsn = poly.path, layer = poly.layer.adm2, options = "ENCODING=UTF-8")
poly.adm2$admin2 <- 1:nrow(poly.adm2)
poly.adm2$admin2.char <- paste0('admin2_',1:nrow(poly.adm2))
  
admin1.mat <- spdep::poly2nb(sp::SpatialPolygons(sf::as_Spatial(poly.adm1)@polygons))
admin1.mat <- spdep::nb2mat(admin1.mat, zero.policy = TRUE)
row.names(admin1.mat) <- colnames(admin1.mat) <- poly.adm1$admin1.char
n_admin1 <- nrow(admin1.mat)
  
admin2.mat <- spdep::poly2nb(sp::SpatialPolygons(sf::as_Spatial(poly.adm2)@polygons))
admin2.mat <- spdep::nb2mat(admin2.mat, zero.policy = TRUE)
row.names(admin2.mat) <- colnames(admin2.mat) <- poly.adm2$admin2.char
n_admin2 <- nrow(admin2.mat)

admin.key <- as.data.frame(poly.adm2[,c('NAME_1','NAME_2','admin2','admin2.char')])[,1:4]
admin.key <- merge(admin.key,poly.adm1[,c('NAME_1','admin1','admin1.char')])[,1:6]
admin.key <- admin.key[order(admin.key$admin2),] %>% rename(admin1.name=NAME_1,admin2.name=NAME_2)

#make admin2 mat nested
admin2.mat.nested <- admin2.mat
for(i in 1:nrow(admin2.mat.nested)){
  admin2.mat.nested[i,which(admin.key$admin1!=admin.key$admin1[i])] <- 0
  if(sum(admin2.mat.nested[i,])>0){
    admin2.mat.nested[i,] <- admin2.mat.nested[i,]/sum(admin2.mat.nested[i,])
  }
}

# load NMR data 
load('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Zambia/2009_2013_Zambia_cluster_dat.rda')


# load urban rural weights 
adm1_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_ZMB/admin1_u1_2011_urban_frac.rds"))
adm1_UR_weights <- adm1_UR_weights[order(as.numeric(str_remove(adm1_UR_weights$adm_idx,'admin1_'))),]
adm2_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_ZMB/admin2_u1_2011_urban_frac.rds"))
adm2_UR_weights <- adm2_UR_weights[order(as.numeric(str_remove(adm2_UR_weights$adm_idx,'admin2_'))),]

# admin1 and admin2
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Zambia/worldpop/adm1_weights_u1.rda"))
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Zambia/worldpop/adm2_weights_u1.rda"))

# get generalized inverse of scaled Q for BYM2 ----

Q.admin2 <- -admin2.mat
Q.admin2 <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})*Q.admin2
diag(Q.admin2) <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})
if(sum(diag(Q.admin2)==0)>0){
  stop('Admin 2 Island')
}
Q_scaled <- INLA::inla.scale.model(Q.admin2, constr=list(A=matrix(1,nrow=1,ncol=n_admin2), e=0))
Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))

Q.admin2.nested <- -admin2.mat.nested
Q.admin2.nested <- sapply(1:nrow(Q.admin2.nested),function(i){sum(I(Q.admin2.nested[i,]!=0))})*Q.admin2.nested
diag(Q.admin2.nested) <- sapply(1:nrow(Q.admin2.nested),function(i){sum(I(Q.admin2.nested[i,]!=0))})
constraints <- matrix(0,nrow=n_admin1,ncol=n_admin2) 
for(a in 1:n_admin2){constraints[admin.key$admin1[admin.key$admin2==a],a] <- 1}
Q_nested_scaled <- INLA::inla.scale.model(Q.admin2.nested, constr=list(A=constraints, e=rep(0,n_admin1)))
Q_nested_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_nested_scaled))

# direct estimates -----
data_for_direct <- nmr.dat
data_for_direct$age <- 0
data_for_direct$years <- 'All'
data_for_direct$died <- data_for_direct$Y
dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                             regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                             Ntrials = "total",weightsVar = "v005")

 sd.fit <- SUMMER::smoothDirect(dir.est, Amat=admin1.mat, time.model = NULL)
 sd.est <- SUMMER::getSmoothed(sd.fit)
 sd.est$admin1 <- as.numeric(unlist(lapply(str_split(sd.est$region,'_'),function(x){x[2]})))

natl.dir <- dir.est[dir.est$region=='All',]
admin1.dir <- dir.est[dir.est$region!='All',]
admin1.dir <- admin1.dir %>% dplyr::select(region,mean,lower,upper,logit.est,var.est) %>% rename(admin1_mean = mean, admin1_lower = lower,admin1_upper=upper,admin1.char = region)
admin1.dir$admin1 <- as.numeric(sapply(1:nrow(admin1.dir),function(i){str_split(admin1.dir$admin1.char[i],'_')[[1]][2]}))
admin1.dir <- admin1.dir[order(admin1.dir$admin1),]

# input strata level DHS survey info ----
dhs_dat <- data.frame(admin1 = rep(sort(unique(admin.key$admin1)),each=2),U = rep(c(1,0),max(admin.key$admin1)))
dhs_dat <- left_join(dhs_dat,unique.array(admin.key[,c('admin1','admin1.char','admin1.name')])) %>% rename(A1=admin1)
# total number of HH in frame (from 2018 DHS report)
dhs_dat$M <- c(76002,198744,336672,90217,47371,295534,44254,199656,422029,92051,26585,127665,31460,110464,
               44296,196260,79551,206791,27196,163099)
# total number of clusters in strata (from 2018 DHS report)
dhs_dat$n_clusters <- c(593,2234,2351,924,245,3279,334,1890,2767,833,188,1470,198,984,297,2207,541,2307,214,1775)
# average number of HH per cluster
dhs_dat$M_bar <- round(dhs_dat$M/dhs_dat$n_clusters)
dhs_dat$n_HH_samp <- 25

# estimate cluster size -----
frame.info <- merge(dhs_dat,nmr.dat %>% mutate(U=ifelse(urban=='urban',1,0)) %>%  group_by(admin1.char,U) %>% summarise(n_clusters_samp = n())) %>% 
  dplyr::select(A1,admin1.char,U,M_bar,n_clusters,n_clusters_samp,n_HH_samp)
sd_M <- 0.15*frame.info$M_bar

n_iter <- 1000
M_obs_bar <- rep(NA,n_iter)
frame.info$M_bar_est <- NA

for(i in 1:nrow(frame.info)){
  for(k in 1:n_iter){
    M = round(rnorm(frame.info$n_clusters[i],frame.info$M_bar[i],sd_M[i]))
    M_obs <- M[sample(1:frame.info$n_clusters[i],frame.info$n_clusters_samp[i],prob = M)]
    M_obs_bar[k] <- mean(M_obs)
  }
  frame.info$M_bar_est[i] <- mean(M_obs_bar)
}

# organize population and survey data for algorithm ----
alg_dat <- nmr.dat %>% mutate(U=ifelse(urban=='urban',1,0)) %>% 
  dplyr::select(cluster,total,Y,U,admin1,admin1.char,admin1.name,admin2,admin2.char,admin2.name) %>% 
  rename(Z=Y,n=total,A1=admin1, A2=admin2)
# add estimated number of total births per cluster
alg_dat <- merge(alg_dat,frame.info,by=c('admin1.char','A1','U')) %>% dplyr::select(-M_bar,-n_clusters,-n_clusters_samp) %>%
  mutate(N=round(M_bar_est*n/n_HH_samp))

#figure out 'best' way to get estimates of total number of births per strata --- use worldpop U1 and then rescale using admin1 births from census
pop_fracs <- NULL 
for(year in 2009:2013){
  adm2_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_ZMB/admin2_u1_",year,"_urban_frac.rds"))
  pop_fracs_tmp <- data.frame(admin2.char = rep(admin.key$admin2.char,each=2),U=rep(c(0,1),max(admin.key$admin2)))
  pop_fracs_tmp <- merge(pop_fracs_tmp,weight.adm2.u1[weight.adm2.u1$year==year,c(1,4)])
  pop_fracs_tmp <- merge(pop_fracs_tmp,adm2_UR_weights[,2:3],by.y='adm_idx',by.x='admin2.char')
  pop_fracs_tmp <- pop_fracs_tmp %>% mutate(full_prop=ifelse(U==1,proportion*urb_frac,proportion*(1-urb_frac)))
  pop_fracs_tmp$year <- year
  
  pop_fracs <- rbind(pop_fracs,pop_fracs_tmp[,c(1,2,5,6)])
}
pop_fracs <- merge(admin.key,pop_fracs)

# https://zambia.opendataforafrica.org/kvfbjpd/annual-number-of-births 
subnational_births <- data.frame(admin1.char = rep(unique(admin.key$admin1.char),5),
                                 admin1.name = rep(unique(admin.key[order(admin.key$admin1),]$admin1.name),5),
                                 year = rep(2009:2013,each=n_admin1),
                                 births = c(rep(NA,2*n_admin1),
                                   60731,81678,77200,51237,94465,38110,35111, 55726,73979,39839, #2011
                                            61775, 83311, 78501, 51766, 98514, 39198, 35545, 56814, 75475, 39834, #2012
                                            62828, 84980, 79776, 52284, 102443, 40233, 35978, 57888, 76902, 39840)) #2013

# estimate subnational births in 2009 and 2010 (no data)
for(area in unique(subnational_births$admin1.char)){
  change <- (subnational_births[subnational_births$admin1.char==area & subnational_births$year==2013,]$births - 
    subnational_births[subnational_births$admin1.char==area & subnational_births$year==2011,]$births)/2
  subnational_births[subnational_births$admin1.char==area & subnational_births$year==2010,]$births <- 
    subnational_births[subnational_births$admin1.char==area & subnational_births$year==2011,]$births - change
  subnational_births[subnational_births$admin1.char==area & subnational_births$year==2009,]$births <- 
    subnational_births[subnational_births$admin1.char==area & subnational_births$year==2011,]$births - 2*change
}

pop_dat <- NULL
for(area in unique(admin.key$admin1.char)){
  for(year_t in 2009:2013){
    pop_dat_est <- pop_fracs %>% filter(admin1.char==area & year==year_t)
    pop_dat_est$N <- 
      pop_dat_est$full_prop*subnational_births[subnational_births$admin1.char==area & subnational_births$year==year_t,]$births/sum(pop_dat_est$full_prop)
    pop_dat <- rbind(pop_dat,pop_dat_est)
  }
}

pop_strata_dat <- pop_dat %>% group_by(admin2.char,admin1.char,admin1,admin2,U) %>% summarise(N = round(sum(N)))  %>% rename(A1=admin1,A2=admin2)
pop_strata_dat <- pop_strata_dat[order(pop_strata_dat$A2,pop_strata_dat$U),]

# some of these don't match up -- urban/rural fraction suggests the area is 100% rural, but there are actually some urban clusters
# checked map to see if any look like they could have been jittered over admin2 lines, but does not appear to be the case
# one way to interpret this is that it means this areas very recently became urban, so we will change them to rural 
alg_dat %>% filter(admin2.char %in% c(pop_strata_dat[pop_strata_dat$N==0& pop_strata_dat$U==1,]$admin2.char) & U==1)
alg_dat <- alg_dat %>% mutate(U = ifelse(admin2.char %in% c(pop_strata_dat[pop_strata_dat$N==0,]$admin2.char),0,U))

# INLA model (for scaling matrix) ------

inla.fit <- INLA::inla(Y ~ urban + factor(admin1) + f(admin2, model="bym2", graph = admin2.mat.nested, scale.model=T, constr=T,
                                                          adjust.for.con.comp=T, 
                                                          hyper=list(phi=list(prior="pc", param=c(0.5, 2/3), initial=1),
                                                                     prec=list(prior="pc.prec",param=c(1,0.01),initial = 5))) -1,
                           family = 'nbinomial', data = nmr.dat, 
                           control.family=list(link='log'),
                           control.fixed = list(mean =list("factor(U)0"=-3.5, "factor(U)1"=-3.5, default = 0),
                                            prec=list(default = 1/9)),
                           control.predictor = list(compute = FALSE, link = 1), 
                           control.compute = list(config=T),
                           E = nmr.dat$total)

# sample from posterior
cs <- inla.fit$misc$configs$contents$tag
cs <- cs[cs != "Predictor"]
select <- list()
for (i in 1:length(cs)) {
  select[[i]] <- 0
  names(select)[i] <- cs[i]
}
sampFull <- INLA::inla.posterior.sample(n = 1000, result = inla.fit, intern = TRUE, selection = select)
#get estimate of covariance matrix
sampFull.draws <- matrix(NA,1000,4+max(admin.key$admin1)+max(admin.key$admin2))
for(i in 1:1000){
  sampFull.draws[i,] <- c(sampFull[[i]]$hyperpar[2:3],
                          -sampFull[[i]]$hyperpar[1], #kind of comparable to log(d)
                          sampFull[[i]]$latent[(max(2*admin.key$admin2)+1):(max(2*admin.key$admin2)+2)], # urban rural intercepts
                          tail(sampFull[[i]]$latent,max(admin.key$admin1)-1),
                          sampFull[[i]]$latent[c(1:max(admin.key$admin2))]) # admin2 spatial effect
}
Sigma <- cov(sampFull.draws)

# sampFull.draws <- matrix(NA,1000,length(sampFull[[1]]$latent))
# for(i in 1:1000){
#   sampFull.draws[i,] <- sampFull[[i]]$latent
# }
# colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)
# sampFull.draws <- sampFull.draws[,-c((n_admin2+1):(2*n_admin2))] # extra remove BYM2 components
# fields <- colnames(sampFull.draws)
# 
# eta_samples <- matrix(0,ncol=2*n_admin2,nrow=1000)
# for(area in 1:n_admin2){
#   adm1_area_tmp <- admin.key[admin.key$admin2==area,]$admin1
#   # urban rate
#   eta_samples[,area] <- sampFull.draws[,n_admin2+2] + sampFull.draws[,area]
#   # rural rate
#   eta_samples[,area + n_admin2] <- sampFull.draws[,n_admin2+1] + sampFull.draws[,area]
#   # add admin1 fixed effect if not in admin1_1
#   if(adm1_area_tmp!=1){
#     eta_samples[,area] <-  eta_samples[,area] + sampFull.draws[,n_admin2+1+adm1_area_tmp]
#     eta_samples[,area + n_admin2] <- eta_samples[,area + n_admin2] + sampFull.draws[,n_admin2+1+adm1_area_tmp]
#   }
# }
# 
# admin2_strat_samples <- exp(eta_samples)
# admin2_samples <- admin2_strat_samples%*%adm2_wt_mat
# 
# admin2.inla <- data.frame(A2 = 1:n_admin2, admin2.char = paste0('admin2_',1:n_admin2),
#                           inla.mean = apply(admin2_samples,2,mean),
#                           inla.median = apply(admin2_samples,2,median),
#                           inla.sd = apply(admin2_samples,2,sd),
#                           inla.lower = apply(admin2_samples,2,quantile,prob = 0.05),
#                           inla.upper = apply(admin2_samples,2,quantile,prob = 0.95))
# admin2.inla <- merge(admin.key,admin2.inla)
# 
# admin1_samples <- admin2_samples%*%adm1_wt_mat
# admin1.inla <- data.frame(admin1 = 1:n_admin1, admin1.char = paste0('admin1_',1:n_admin1),
#                           inla.mean = apply(admin1_samples,2,mean),
#                           inla.median = apply(admin1_samples,2,median),
#                           inla.sd = apply(admin1_samples,2,sd),
#                           inla.lower = apply(admin1_samples,2,quantile,prob = 0.05),
#                           inla.upper = apply(admin1_samples,2,quantile,prob = 0.95))
# 

# Stan models ------

standardmod.nofe <- cmdstan_model('Stan models/mod3.stan')

list_dat <- list(lenA2=n_admin2,
                 lenC=nrow(alg_dat),
                 urban_id=alg_dat$U,
                 admin2_id=alg_dat$A2,
                 Y=alg_dat$Z,
                 N=alg_dat$n,
                 Q_scaled_inv = Q_scaled_inv)

stan.fit1 <- standardmod.nofe$sample(
  data = list_dat,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh=500)

out1 <- stan.fit1$draws(format='df',variables=c('alpha','b'))
stan.nofe.samples <- exp(cbind(out1$`alpha[1]` + out1[,3:117],out1$`alpha[2]` + out1[,3:117]))

standardmod <- cmdstan_model('Simulation_Study/Model_compare.stan')

list_dat <- list(lenA1=n_admin1, lenA2=n_admin2,
                 lenC=nrow(alg_dat),
                 admin1_id=alg_dat$A1,
                 admin2_id=alg_dat$A2,
                 urban_id=alg_dat$U,
                 Y=alg_dat$Z,
                 N=alg_dat$n,
                 Q_scaled_inv = Q_nested_scaled_inv)

stan.fit2 <- standardmod$sample(
  data = list_dat,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh=500)

out2 <- as.data.frame(stan.fit2$draws(format='df',variables=c('alpha','beta','b')))
eta_postsamp <- cbind(out2$`alpha[1]` + out2[,(2+n_admin1):(1+n_admin1+n_admin2)],
                      out2$`alpha[2]` + out2[,(2+n_admin1):(1+n_admin1+n_admin2)])
for(i in 1:nrow(admin.key)){
  if(admin.key$admin1[i]>1){
    eta_postsamp[,i] <- eta_postsamp[,i] + out2[,admin.key$admin1[i]+1] 
    eta_postsamp[,n_admin2 + i] <- eta_postsamp[,n_admin2 + i] + out2[,admin.key$admin1[i]+1] 
  }
}
stan.samples <- exp(eta_postsamp)

adm2_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_ZMB/admin2_u1_2011_urban_frac.rds"))
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Zambia/worldpop/adm2_weights_u1.rda"))
adm2_weights <- weight.adm2.u1[weight.adm2.u1$year==2011,]
adm2_UR_weights <- adm2_UR_weights[order(as.numeric(str_remove(adm2_UR_weights$adm_idx,'admin2_'))),]
adm2_UR_weights$A2 <- 1:nrow(adm2_UR_weights)
adm2_weights <- adm2_weights[order(as.numeric(str_remove(adm2_weights$admin2.char,'admin2_'))),]
adm2_weights$A2 <- 1:nrow(adm2_weights)

adm2_wt_mat <- matrix(0,2*n_admin2,n_admin2)
for(area in 1:n_admin2){
  adm2_wt_mat[area,area] <- adm2_UR_weights[adm2_UR_weights$A2==area,]$urb_frac #urban
  adm2_wt_mat[area+n_admin2,area] <- 1-adm2_UR_weights[adm2_UR_weights$A2==area,]$urb_frac # rural
}

adm1_wt_mat <- matrix(0,n_admin2,n_admin1)
for(area in 1:n_admin1){
  adm2_areas_tmp <- admin.key[admin.key$admin1==area,]$admin2
  adm1_wt_mat[adm2_areas_tmp,area] <- adm2_weights[adm2_weights$A2 %in% adm2_areas_tmp,]$proportion/sum(adm2_weights[adm2_weights$A2 %in% adm2_areas_tmp,]$proportion)
}

admin2.samples <- as.matrix(stan.nofe.samples)%*%adm2_wt_mat

admin2.res.stan.nofe <- data.frame(admin2 = 1:n_admin2, 
                              stan.nofe.mean = apply(admin2.samples,2,mean),
                              stan.nofe.median = apply(admin2.samples,2,median),
                              stan.nofe.sd = apply(admin2.samples,2,sd),
                              stan.nofe.lower = apply(admin2.samples,2,quantile,prob = 0.05),
                              stan.nofe.upper = apply(admin2.samples,2,quantile,prob = 0.95))
admin2.res.stan.nofe <- merge(admin.key,admin2.res.stan.nofe)

admin1.samples <- admin2.samples%*%adm1_wt_mat
admin1.res.stan.nofe <- data.frame(admin1 = 1:n_admin1, 
                              stan.nofe.mean = apply(admin1.samples,2,mean),
                              stan.nofe.median = apply(admin1.samples,2,median),
                              stan.nofe.sd = apply(admin1.samples,2,sd),
                              stan.nofe.lower = apply(admin1.samples,2,quantile,prob = 0.05),
                              stan.nofe.upper = apply(admin1.samples,2,quantile,prob = 0.95))

admin2.samples <- as.matrix(stan.samples)%*%adm2_wt_mat
admin2.res.stan <- data.frame(admin2 = 1:n_admin2, 
                             stan.mean = apply(admin2.samples,2,mean),
                             stan.median = apply(admin2.samples,2,median),
                             stan.sd = apply(admin2.samples,2,sd),
                             stan.lower = apply(admin2.samples,2,quantile,prob = 0.05),
                             stan.upper = apply(admin2.samples,2,quantile,prob = 0.95))
admin2.res.stan <- merge(admin.key,admin2.res.stan)

admin1.samples <- admin2.samples%*%adm1_wt_mat
admin1.res.stan <- data.frame(admin1 = 1:n_admin1, 
                             stan.mean = apply(admin1.samples,2,mean),
                             stan.median = apply(admin1.samples,2,median),
                             stan.sd = apply(admin1.samples,2,sd),
                             stan.lower = apply(admin1.samples,2,quantile,prob = 0.05),
                             stan.upper = apply(admin1.samples,2,quantile,prob = 0.95))

# shrinkage
plot(admin2.res.stan$stan.median,admin2.res.stan.nofe$stan.nofe.median,ylim=range(admin2.res.stan$stan.median))
plot(admin1.res.stan$stan.median,admin1.res.stan.nofe$stan.nofe.median,ylim=range(admin1.res.stan$stan.median))
# including admin1 FE helps with shrinkage, but doesn't necessarily bring closer to direct estimate
plot(admin1.dir$admin1_mean,admin1.res.stan$stan.median)
###

## run model -------

ZMB_data <- list()
ZMB_data[[1]] <- list(M_init = solve(Sigma),
                      logit_r_hat = admin1.dir$logit.est,
                      logit_r_hat_var = admin1.dir$var.est,
                      data_list =  list(obs_dat = alg_dat,
                                        pop_strata_dat = as.data.frame(pop_strata_dat),
                                        Q_scaled_inv=Q_nested_scaled_inv))

save(ZMB_data,file = '/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Country Examples/ZMB_data.rda')

# NUTS_5.4 <- postsamp_Alg5.4_NUTS_MCMC_onechain(eps0 = 0.04, 
#                                         M_init = solve(Sigma),
#                                         logit_r_hat = admin1.dir$logit.est,
#                                         logit_r_hat_var = admin1.dir$var.est,
#                                         data_list =  list(obs_dat = alg_dat,
#                                                           pop_strata_dat = as.data.frame(pop_strata_dat),
#                                                           Q_scaled_inv=Q_nested_scaled_inv), 
#                                         n_iter = 2000)

#save(NUTS_5.4,file="Handcoded MCMC Simulations/Alg 5.4 240416, ZMB 2009-2013, 4 chains of 3k iterations.rda")
#load("Handcoded MCMC Simulations/Alg 5.4 240405, ZMB 2009-2013, 4 chains of 2k iterations.rda")

ZMB_res <- NULL
for(chain in 1:4){
  res.tmp <- read.xlsx(paste0('Country Examples/ZMB/Results',chain,'.xlsx'))
  res.tmp$iter <- 1:nrow(res.tmp)
  ZMB_res <- rbind(ZMB_res, res.tmp)
}

## Diagnostic plots ------

plot(ZMB_res[ZMB_res$iter>1000,]$alphaU, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$alphaR, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$logtau, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$logitphi, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$logd, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$beta2, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$beta6, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$beta9, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$b1, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$b17, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$Yplus1, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$Yplus2, type='l')
plot(ZMB_res[ZMB_res$iter>1000,]$Yplus9, type='l')

each_chain <- 3000
par(mfrow=c(1,3))
for(k in 1:n_admin1){
  plot(cumsum(ZMB_res[1:each_chain,360 + k])/(1:each_chain),type='l',col=1,ylab='',
       ylim=c(min(ZMB_res[,360 + k]),max(ZMB_res[,360+k])),main = paste0('Running mean of Yplus_',k))
  for(c in 2:4){
    lines(cumsum(ZMB_res[((c-1)*each_chain+1):(c*each_chain),360+k])/(1:each_chain),col=c)
  }
}

# extract NUTS results ------

raw.rates <- ZMB_res %>% filter(iter>1000 & iter < 2001) %>% dplyr::select(rU1:rR115)

admin2.samples <- as.matrix(raw.rates)%*%adm2_wt_mat

admin2.res.bench <- data.frame(A2 = 1:n_admin2, 
                             bench.mean = apply(admin2.samples,2,mean),
                             bench.median = apply(admin2.samples,2,median),
                             bench.sd = apply(admin2.samples,2,sd),
                             bench.lower = apply(admin2.samples,2,quantile,prob = 0.05),
                             bench.upper = apply(admin2.samples,2,quantile,prob = 0.95))
admin2.res.bench$admin2 <- paste0('admin2_',admin2.res.bench$A2)

admin1.samples <- admin2.samples%*%adm1_wt_mat
admin1.res.bench <- data.frame(admin1 = 1:n_admin1, 
                             bench.mean = apply(admin1.samples,2,mean),
                             bench.median = apply(admin1.samples,2,median),
                             bench.sd = apply(admin1.samples,2,sd),
                             bench.lower = apply(admin1.samples,2,quantile,prob = 0.05),
                             bench.upper = apply(admin1.samples,2,quantile,prob = 0.95))
admin1.res.bench$admin1.char <- paste0('admin1_',admin1.res.bench$admin1)

## plots -----

area.order <- data.frame(admin1.ordered = 1:n_admin1,
                         admin1 = admin1.dir$admin1[order(admin1.dir$admin1_mean)])
area.order <- merge(area.order,admin.key %>% dplyr::select(admin1.name,admin1) %>% unique(),by='admin1')
area.order <- area.order[order(area.order$admin1.ordered),]

admin1.dir <- merge(admin1.dir,area.order)
admin1.res <- merge(admin1.dir,admin1.res.bench)
#admin1.res <- merge(admin1.res,admin1.res.inla)
admin1.res <- merge(admin1.res,admin1.res.stan)

admin2.res <- merge(admin.key,admin2.res.bench,by.x='admin2',by.y = 'A2')
admin2.res <- merge(area.order, admin2.res)
#admin2.res <- merge(admin2.res,admin2.res.inla)
admin2.res <- merge(admin2.res,admin2.res.stan)

admin2.res.stan <- merge(admin2.res.stan,area.order)
admin1.res.stan <- merge(admin1.res.stan,area.order)
admin2.res.stan.nofe <- merge(admin2.res.stan.nofe,area.order)
admin1.res.stan.nofe <- merge(admin1.res.stan.nofe,area.order)

colors <- hue_pal()(4)

# Map plots summarizing areas/data ====

ggplot() +  geom_sf(data=poly.adm2,aes(fill=NAME_1)) + theme_minimal() +
  #geom_sf_text(data=poly.adm1,aes(label=NAME_1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.position = 'bottom') +
  scale_fill_discrete(name='First Adminstrative Area')

# plot of data
poly.adm1 <- merge(poly.adm1,nmr.dat %>% group_by(admin1.char) %>% reframe(deaths = sum(Y),births = sum(total)))
poly.adm2 <- merge(poly.adm2,nmr.dat %>% group_by(admin2.char) %>% reframe(deaths = sum(Y),births = sum(total)))

b1 <- ggplot() +  geom_sf(data=poly.adm1,aes(fill=births)) + theme_bw() +
  geom_sf_label(data=poly.adm1,aes(label=NAME_1),alpha=0.75, 
                label.size = NA,label.padding = unit(0.15, "lines"), size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.position = 'top') +
  scale_fill_gradientn(name='Number of observed births',
                       limits = c(5,1600),
                       breaks = c(5,500,1000,1500),
                       values=rescale(c(5,1000,1600)),
                       colors = c('red','white','blue')) +
  ggtitle('First administrative area')

b2 <- ggplot() +  geom_sf(data=poly.adm2,aes(fill=births)) + theme_bw() +
  geom_sf(fill = "transparent", color = "grey30", lwd=0.75, data = poly.adm2 %>% group_by(NAME_1) %>% summarise()) +
  geom_sf_label(data=poly.adm1,aes(label=NAME_1),alpha=0.75, 
                label.size = NA,label.padding = unit(0.15, "lines"), size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.position = 'top') +
  scale_fill_gradientn(name='Number of observed births',
                       limits = c(5,1600),
                       breaks = c(5,500,1000,1500),
                       values=rescale(c(5,1000,1600)),
                       colors = c('red','white','blue')) +
  ggtitle('Second administrative area')

ggarrange(plotlist = list(b1,b2),common.legend = T)

d1 <- ggplot() +  geom_sf(data=poly.adm1,aes(fill=deaths)) + theme_bw() +
  geom_sf_label(data=poly.adm1,aes(label=NAME_1),alpha=0.75, 
                label.size = NA,label.padding = unit(0.15, "lines"), size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.position = 'top') +
  scale_fill_gradientn(name='Number of observed deaths',
                        limits = c(0,55),
                        breaks = c(0,10,20,30,40,50),
                       values=rescale(c(0,20,55)),
                       colors = c('red','white','blue')) +
  ggtitle('First administrative area')

d2 <- ggplot() +  geom_sf(data=poly.adm2,aes(fill=deaths)) + theme_bw() +
  geom_sf(fill = "transparent", color = "grey30", lwd=0.75, data = poly.adm2 %>% group_by(NAME_1) %>% summarise()) +
  geom_sf_label(data=poly.adm1,aes(label=NAME_1),alpha=0.75, 
                label.size = NA,label.padding = unit(0.15, "lines"), size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.position = 'top')+
  scale_fill_gradientn(limits = c(0,55),
                       values=rescale(c(0,20,55)),
                       colors = c('red','white','blue')) + 
  ggtitle('Second administrative area')

ggarrange(plotlist = list(d1,d2),common.legend = T)

# Aggregation plots ======

# 1 = direct
# 2 = No FE admin2
# 3 = No FE admin1
# 4 = FE admin2
# 5 = FE admin1
# 6 = bench admin2
# 7 = bench admin1

## motivating example
ggplot() + geom_point(data=admin2.res.stan.nofe,aes(1000*stan.nofe.median,admin1.ordered-0.15,col='2',pch='2'),size=2) +
  geom_point(data=admin1.res.stan.nofe,aes(1000*stan.nofe.median,admin1.ordered + 0.01,col='3',pch='3'),size=4,) +
  geom_point(data=admin1.dir,aes(1000*admin1_mean,admin1.ordered,col='1',pch='1'),size=4) + # admin1 direct
  geom_point(data=admin2.res.stan,aes(1000*stan.median,admin1.ordered-0.25,col='4',pch='4'),size=2) +
  geom_point(data=admin1.res.stan,aes(1000*stan.median,admin1.ordered - 0.01,col='5',pch='5'),size=4) +
 # ggtitle('Zambia 2009-2013') + 
  theme_bw() +
  ylab('') + xlab('Neo-natal deaths per 1000 live births') + theme(legend.position = 'bottom',
                                 axis.text = element_text(size=12),
                                 legend.text = element_text(size=10)) +
  scale_y_continuous(breaks=1:n_admin1,labels = area.order$admin1.name) +
  scale_colour_manual(name = '', values =c('1'='black',
                                           '5'=colors[1],'3'=colors[4],
                                           '4'=colors[1],'2'=colors[4]),
                      labels = c('1' = 'Direct admin 1',
                                 '4' = 'Nested admin 2', '5' = 'Nested (aggregated to admin1)',
                                 '2' = 'Non-nested admin 2', '3' = 'Non-nested (aggregated to admin1)')) +
  scale_shape_manual(name = '', values =c('1'=5,'5'=5,'3'=5,'4'=1,'2'=1),
                     labels = c('1' = 'Direct admin 1',
                                '4' = 'Nested admin 2', '5' = 'Nested (aggregated to admin1)',
                                '2' = 'Non-nested admin 2', '3' = 'Non-nested (aggregated to admin1)')) +
  guides(color=guide_legend(nrow=2,reverse = T),pch=guide_legend(nrow=2,reverse = T))

# with DABUL
ggplot() + geom_point(data=admin1.res,aes(1000*admin1_mean,admin1.ordered,col='1',pch='1'),size=4) + # admin1 direct
  geom_point(data=admin2.res,aes(1000*bench.median,admin1.ordered-0.15,col='6',pch='6'),size=2) +
  geom_point(data=admin1.res,aes(1000*bench.median,admin1.ordered + 0.01,col='7',pch='7'),size=4,) +
  geom_point(data=admin2.res,aes(1000*stan.median,admin1.ordered-0.25,col='4',pch='4'),size=2) +
  geom_point(data=admin1.res,aes(1000*stan.median,admin1.ordered - 0.01,col='5',pch='5'),size=4) + # admin1 aggregated INLA
  #ggtitle('Zambia 2009-2013') + 
  theme_bw() +
  ylab('') + xlab('Neo-natal deaths per 1000 live births') + theme(legend.position = 'bottom',
                                 axis.text = element_text(size=12),
                                 legend.text = element_text(size=10)) +
  scale_y_continuous(breaks=1:n_admin1,labels = area.order$admin1.name) +
  scale_colour_manual(name = '', values =c('1'='black','4'=colors[1],'5'=colors[1],'6'=colors[3],'7'=colors[3]),
                      labels = c('1' = 'Direct admin 1',
                                 '4' = 'Standard unit-level nested admin 2', '5' = 'Standard unit-level nested (aggregated to admin 1)',
                                 '6' = 'DABUL admin 2','7' = 'DABUL (aggregated to admin 1)')) +
  scale_shape_manual(name = '', values =c('1'=5,'4'=1,'5'=5,'6'=1,'7'=5),
                     labels = c('1' = 'Direct admin 1',
                               '4' = 'Standard unit-level nested admin 2', '5' = 'Standard unit-level nested (aggregated to admin 1)',
                               '6' = 'DABUL admin 2','7' = 'DABUL (aggregated to admin 1)')) +
  guides(color=guide_legend(nrow=2,reverse = T),pch=guide_legend(nrow=2,reverse = T))

# Scatter plots =====

# Compare medians

admin2.res %>% ggplot() + 
 geom_point(aes(1000*bench.median,1000*stan.median,color=factor(admin1.name)),pch=1,size=3) +
 geom_abline(slope=1,intercept = 0) +
 theme_bw() +
 ylab('Standard unit-level model') + xlab('DABUL model') +
 theme(legend.position = 'bottom',  
       axis.title = element_text(size=12),
       axis.text = element_text(size=10), 
       legend.title = element_text(size=12),
       legend.text = element_text(size=10)) +
  scale_color_discrete(name = 'First administrative areas') +
  scale_y_continuous(limits = 1000*range(admin2.res$bench.median,admin2.res$stan.median) + c(-1,1),
                     breaks = c(20,25,30)) +
  scale_x_continuous(limits = 1000*range(admin2.res$bench.median,admin2.res$stan.median)+ c(-1,1),
                     breaks = c(20,25,30))

# Comparing precision
admin2.order <- data.frame(admin2.ordered = 1:n_admin2,
                         admin2 = admin2.res$admin2[order(admin2.res$bench.median)])
admin2.res <- merge(admin2.res,admin2.order)

admin2.res %>% ggplot() + 
  geom_line(aes(admin2.ordered,bench.sd/bench.median,col='Bench')) +
  geom_line(aes(admin2.ordered,stan.sd/stan.median,col='No Bench')) +
  theme_bw() +
  ylab('Coefficient of variation') + xlab('Second administrative area (in order of NMR estimate)') +
  theme(legend.position = 'bottom', axis.text.x = element_blank(), 
        axis.title = element_text(size=12),
        axis.text.y = element_text(size=12), legend.text = element_text(size=10)) +
  scale_y_continuous(breaks = c(0.15,0.2,0.25,0.3,0.35)) +
  scale_colour_manual(name = '', values =c('No Bench'=colors[1],'Bench'=colors[3]),
                      labels=c('No Bench'='Standard unit-level nested','Bench'='DABUL'))

# Compare lower bound of confidence interval
ggplot() + geom_point(aes(1000*admin2.res.stan$stan.lower,1000*admin2.res.bench$bench.lower), pch=3) +
  theme_bw() +
  xlab('90% lower bound on standard model estimates') +
  ylab('90% lower bound on DABUL estimates') +
  geom_abline(slope=1,intercept=0) +
  geom_hline(yintercept = 12, col='blue') + geom_vline(xintercept = 12, col='blue') +
  scale_x_continuous(breaks=c(10,12,15,20,25)) +
  scale_y_continuous(breaks=c(10,12,15,20,25)) +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=10), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=10))


# Map plots of NMR =====
poly.adm2$admin2 <- 1:nrow(poly.adm2)
poly.adm2$admin2.char <- paste0('admin2_',poly.adm2$admin2)
poly.adm2 <- merge(poly.adm2,admin2.res)
poly.adm2 <- merge(poly.adm2,admin1.dir,by.x='NAME_1',by.y='admin1.name')
poly.adm2 <- merge(poly.adm2,admin2.res.stan.nofe)
poly.adm2 <- merge(poly.adm2,alg_dat %>% group_by(admin2.char) %>% reframe(births = sum(n),deaths = sum(Z)))

map_template <- ggplot() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size=14),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.position = 'bottom') +
  scale_fill_gradientn(name='Neo-natal deaths per 1000 live births',colors=c(colors[3],colors[4],colors[1]),
                       limits=1000*range(poly.adm2$stan.median,poly.adm2$bench.median),
                       breaks = c(20,25,30))

# map plots of NMR estimates
unbench.map.nofe <-  map_template +  geom_sf(data=poly.adm2,aes(fill=1000*stan.nofe.median),color='grey80',) + 
  geom_sf(fill = "transparent", size=1, color = "grey50", lwd=0.75, data = poly.adm2 %>% group_by(NAME_1) %>% summarise()) +
  #geom_sf_label(data=poly.adm1,aes(label=NAME_1),alpha=0.5, label.size = NA,label.padding = unit(0.15, "lines")) +
  geom_label_repel(data=poly.adm1,aes(label=NAME_1,geometry=geometry),stat="sf_coordinates",
                   alpha=0.75,label.size=NA,label.padding = unit(0.15, "lines")) +
  ggtitle('Non-nested model')

unbench.map <- map_template +  geom_sf(data=poly.adm2,aes(fill=1000*stan.median),color='grey80') + 
  geom_sf(fill = "transparent", size=1, color = "grey50", lwd=0.75, data = poly.adm2 %>% group_by(NAME_1) %>% summarise()) +
  #geom_sf_label(data=poly.adm1,aes(label=NAME_1),alpha=0.5, label.size = NA,label.padding = unit(0.15, "lines")) +
  geom_label_repel(data=poly.adm1,aes(label=NAME_1,geometry=geometry),stat="sf_coordinates",
                   alpha=0.75,label.size=NA,label.padding = unit(0.15, "lines")) +
  #ggtitle('Nested model')
  ggtitle('Standard unit-level nested model')

dir.map <- map_template + geom_sf(data=poly.adm2,aes(fill=1000*admin1_mean),color='grey80',lwd=0.001) + 
  geom_sf(fill = "transparent", size=1, color = "grey50", lwd=0.75, data = poly.adm2 %>% group_by(NAME_1) %>% summarise()) +
  #geom_sf_label(data=poly.adm1,aes(label=NAME_1),alpha=0.5, label.size = NA,label.padding = unit(0.15, "lines")) +
  geom_label_repel(data=poly.adm1,aes(label=NAME_1,geometry=geometry),stat="sf_coordinates",
                   alpha=0.75,label.size=NA,label.padding = unit(0.15, "lines")) +
  ggtitle('Direct estimates')

bench.map <- map_template +  geom_sf(data=poly.adm2,aes(fill=1000*bench.median),color='grey80') + 
  geom_sf(fill = "transparent", color = "grey50", lwd=0.75, data = poly.adm2 %>% group_by(NAME_1) %>% summarise()) +
  geom_label_repel(data=poly.adm1,aes(label=NAME_1,geometry=geometry),stat="sf_coordinates",
                   alpha=0.75,label.size=NA,label.padding = unit(0.15, "lines")) +
  ggtitle('DABUL model')


ggarrange(plotlist=list(unbench.map.nofe,unbench.map,dir.map),common.legend = T,legend='bottom',nrow=1)
ggarrange(plotlist=list(unbench.map,bench.map),common.legend = T,legend='bottom')




# testing ======

