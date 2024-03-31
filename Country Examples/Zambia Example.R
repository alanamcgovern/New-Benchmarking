library(sf)
library(surveyPrev)
library(SUMMER)
library(tidyverse)
library(data.table)
library(parallel)
library(LaplacesDemon)
library(VGAM)
library(Matrix)
library(Rfast)
library(MASS)
library(locfit)

load(file = paste0('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Info/Zambia_general_info.Rdata')) # load the country info

# load data ----
# adjacency matrix
poly.path <- paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_ZMB_shp")

poly.adm1 <- st_read(dsn = poly.path, layer = poly.layer.adm1, options = "ENCODING=UTF-8")
poly.adm2 <- st_read(dsn = poly.path, layer = poly.layer.adm2, options = "ENCODING=UTF-8")
  
admin1.info <- adminInfo(geo = poly.adm1, admin = 1)
admin1.info$admin.info$admin1 <- 1:nrow(admin1.info$admin.info)
admin1.info$admin.info$admin1.char <- paste0('admin1_',admin1.info$admin.info$admin1)
admin1.mat <- admin1.info$admin.mat
row.names(admin1.mat) <- colnames(admin1.mat) <- admin1.info$admin.info$admin1.char
  
admin2.info <- adminInfo(geo = poly.adm2, admin = 2)
admin2.info$admin.info$admin2 <- 1:nrow(admin2.info$admin.info)
admin2.info$admin.info$admin2.char <- paste0('admin2_',admin2.info$admin.info$admin2)
admin2.mat <- admin2.info$admin.mat
row.names(admin2.mat) <- colnames(admin2.mat) <- admin2.info$admin.info$admin2.char
  
admin.key <- merge(admin1.info$admin.info[,c(1,4,5)],admin2.info$admin.info[,c(1,2,6,7)])
admin.key <- admin.key[order(admin.key$admin2),]
  
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
Q_scaled <- INLA::inla.scale.model(Q.admin2, constr=list(A=t(eigen(Q.admin2)$vectors[,eigen(Q.admin2)$values<1e-10]), e=rep(0,sum(eigen(Q.admin2)$values<1e-10))))
Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))

# direct estimates -----
data_for_direct <- nmr.dat
data_for_direct$age <- 0
data_for_direct$years <- 'All'
data_for_direct$died <- data_for_direct$Y
dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                             regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                             Ntrials = "total",weightsVar = "v005")

# sd.fit <- SUMMER::smoothDirect(dir.est, Amat=admin1.mat, time.model = NULL)
# sd.est <- SUMMER::getSmoothed(sd.fit)

natl.dir <- dir.est[dir.est$region=='All',]
admin1.dir <- dir.est[dir.est$region!='All',]
admin1.dir <- admin1.dir %>% dplyr::select(region,mean,lower,upper,logit.est,var.est) %>% rename(admin1_mean = mean, admin1_lower = lower,admin1_upper=upper,admin1.char = region)
admin1.dir$admin1 <- as.numeric(sapply(1:nrow(admin1.dir),function(i){str_split(admin1.dir$admin1.char[i],'_')[[1]][2]}))
admin1.dir <- admin1.dir[order(admin1.dir$admin1),]
admin1.dir$my.var.est <- (admin1.dir$admin1_upper - admin1.dir$admin1_lower)/(2*1.96)

# input strata level DHS survey info ----
dhs_dat <- data.frame(admin1 = rep(sort(unique(admin.key$admin1)),each=2),U = rep(c(1,0),max(admin.key$admin1)))
dhs_dat <- left_join(dhs_dat,unique.array(admin.key[,c(1:3)])) %>% rename(A1=admin1)
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

#need to figure out 'best' way to get these estimates --- use worldpop U1 and then rescale using admin1 births from census
births_dat <- data.frame(admin1.char = paste0('admin1_',sort(unique(admin.key$admin1))),
                         births = 5*c(60731,81678,77200,51237,94465,38110,35111, 55726,73979,39839))
pop_fracs <- NULL 
for(year in 2009:2013){
  adm2_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_ZMB/admin2_u1_",year,"_urban_frac.rds"))
  adm2_UR_weights <- adm2_UR_weights[order(as.numeric(str_remove(adm2_UR_weights$adm_idx,'admin2_'))),]
  pop_fracs_tmp <- data.frame(admin2.char = rep(admin.key$admin2.char,each=2),U=rep(c(0,1),max(admin.key$admin2)))
  pop_fracs_tmp <- merge(pop_fracs_tmp,weight.adm2.u1[weight.adm2.u1$year==year,c(1,4)])
  pop_fracs_tmp <- merge(pop_fracs_tmp,adm2_UR_weights[,2:3],by.y='adm_idx',by.x='admin2.char')
  pop_fracs_tmp <- pop_fracs_tmp %>% mutate(full_prop=ifelse(U==1,proportion*urb_frac,proportion*(1-urb_frac)))
  pop_fracs_tmp$year <- year
  
  pop_fracs <- rbind(pop_fracs,pop_fracs_tmp[,c(1,2,5,6)])
}
pop_fracs <- merge(admin.key,pop_fracs)

pop_dat <- NULL
for(area in sort(unique(admin.key$admin1))){
  pop_dat_est <- pop_fracs %>% filter(admin1==area)
  pop_dat_est$N <- pop_dat_est$full_prop*births_dat$births[area]/sum(pop_dat_est$full_prop)
  pop_dat <- rbind(pop_dat,pop_dat_est)
}

pop_strata_dat <- pop_dat %>% group_by(admin2.char,U) %>% summarise(N = round(sum(N)))
pop_strata_dat <- merge(admin.key,pop_strata_dat) %>% rename(A1=admin1,A2=admin2)
pop_strata_dat <- pop_strata_dat[order(pop_strata_dat$A2,pop_strata_dat$U),]

# some of these don't match up -- TEMPORARY FIX (this might be because of jittering that they are put into a different admin2 area, will check later)
alg_dat %>% filter(admin2.char %in% c(pop_strata_dat[pop_strata_dat$N==0& pop_strata_dat$U==1,]$admin2.char) & U==1)
alg_dat <- alg_dat %>% mutate(U = ifelse(admin2.char %in% c(pop_strata_dat[pop_strata_dat$N==0,]$admin2.char),0,U))

data_list <- list(obs_dat = alg_dat,
                  pop_strata_dat = as.data.frame(pop_strata_dat),
                  Q_scaled_inv = Q_scaled_inv)


# INLA model ------
inla.fit <- INLA::inla(Y ~ urban + factor(admin1) + f(admin2, model="bym2", graph = admin2.mat.nested, scale.model=T, constr=T,
                                                          adjust.for.con.comp=T, 
                                                          hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1))) -1,
                           family = 'nbinomial', data = nmr.dat, 
                           control.family=list(link='log'),
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
sampFull.draws <- matrix(NA,1000,length(sampFull[[1]]$latent))
for(i in 1:1000){
  sampFull.draws[i,] <- sampFull[[i]]$latent
}
fields <- colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)

admin2.inla <- matrix(0,ncol=(max(admin.key$admin1)+2*max(admin.key$admin2)+1),nrow=max(admin.key$admin2))
for(area in 1:max(admin.key$admin2)){
  admin2.inla[area,area] <- 1 #choose admin1 area
  if(admin.key[admin.key$admin2==area,]$admin1>1){
    admin2.inla[area,2*max(admin.key$admin2)+1+admin.key[admin.key$admin2==area,]$admin1] <- 1 #choose admin1 area
  }
}
admin2.inla[,2*max(admin.key$admin2)+1] <- adm2_UR_weights$urb_frac
admin2.inla[,2*max(admin.key$admin2)+2] <- 1 - adm2_UR_weights$urb_frac

admin2.inla <- admin2.inla%*%t(sampFull.draws)
admin2.res.inla <- data.frame(admin2 = 1:nrow(admin.key))
admin2.res.inla$admin2.char <- paste0('admin2_', admin2.res.inla$admin2)  

admin2.res.inla$adm2_median <- (apply(admin2.inla,1,function(x){median(exp(x))}))
admin2.res.inla$adm2_lower95 <- (apply(admin2.inla,1,function(x){quantile(exp(x),0.025)}))
admin2.res.inla$adm2_upper95 <- (apply(admin2.inla,1,function(x){quantile(exp(x),0.975)}))
admin2.res.inla <- merge(admin.key,admin2.res.inla)

admin1.res.inla <- data.frame(admin1 = rep(sort(unique(admin.key$admin1)),each=2),
                           U = rep(c(0,1),max(admin.key$admin1)))
admin1.res.inla$est <- rep(c(0,inla.fit$summary.fixed$`0.5quant`[3:nrow(inla.fit$summary.fixed)]),each=2) + 
  rep(inla.fit$summary.fixed$`0.5quant`[c(2,1)],max(admin.key$admin1))
admin1.res.inla$admin1.char <- paste0('admin1_',admin1.res.inla$admin1)
admin1.res.inla <- merge(admin1.res.inla,adm1_UR_weights,by.x='admin1.char',by.y='adm_idx')
admin1.res.inla[admin1.res.inla$U==0,]$urb_frac <- 1-admin1.res.inla[admin1.res.inla$U==0,]$urb_frac
admin1.res.inla <- admin1.res.inla %>% group_by(admin1.char,admin1) %>% summarise(adm1_median=exp(sum(est*urb_frac)))

g <- ggplot() + geom_point(data=admin1.dir, aes(admin1_mean,admin1.char,col='1',pch='1'),size=3) +
  geom_point(data=admin1.dir,aes(admin1_lower,admin1.char,col='1'),pch=3,size=3) +
  geom_point(data=admin1.dir,aes(admin1_upper,admin1.char,col='1'),pch=3,size=3) +
  geom_point(data=admin2.res.inla,aes(adm2_median,admin1.char,col='2',pch='2'),size=3) + 
  geom_point(data=admin1.res.inla,aes(adm1_median,admin1.char,col='3',pch='3'),size=3) +
  ggtitle('Zambia 2009-2013') + theme_minimal() +
  ylab('') + xlab('NMR') + theme(legend.position = 'bottom') +
  scale_colour_manual(name = '', values =c('1'='blue','2'='orange','3'='darkgreen'),
                      labels = c('Admin 1 Direct','Admin 2 BYM2 (with Admin 1 FE)','Admin 1 FE')) +
  scale_shape_manual(name = '', values =c('1'=17,'2'=1,'3'=2),
                     labels = c('Admin 1 Direct','Admin 2 BYM2 (with Admin 1 FE)','Admin 1 FE'))

## run model -------

#run model to get estimate of covariance matrix
inla.fit <- INLA::inla(Y ~ urban + f(admin2, model="bym2", graph = admin2.mat.nested, scale.model=T, adjust.for.con.comp=T, 
                                                      hyper=list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1))) -1,
                       family = 'nbinomial', data = nmr.dat, 
                       control.family=list(link='log'),
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
sampFull.draws <- matrix(NA,1000,5+max(admin.key$admin2))
for(i in 1:1000){
  sampFull.draws[i,] <- c(sampFull[[i]]$hyperpar[2:3],
                          -sampFull[[i]]$hyperpar[1], #kind of comparable to log(d)
                          tail(sampFull[[i]]$latent,2),sampFull[[i]]$latent[c(1:max(admin.key$admin2))])
}
Sigma <- cov(sampFull.draws)

start.time <- Sys.time()
NUTS_5.3a <- postsamp_Alg5.3a_NUTS_MCMC(eps0 = 0.1, 
                                        M_init = solve(Sigma),
                                        logit_r_hat = admin1.dir$logit.est,
                                        logit_r_hat_var = admin1.dir$var.est/1000,
                                        data_list = data_list, 
                                        n_iter = 2000, chains=4)
Sys.time() - start.time
NUTS_5.3a$tree_depth
save(NUTS_5.3a,file="Handcoded MCMC Simulations/Alg 5.3a BYM2 NUTS 240322, ZMB 2009-2013, 4 chains of 2k iterations posterior dist.rda")


iter.include <- c(500:1000,1500:2000,2500:3000,3500:4000)
iter.include <- c(500:2000,2500:4000,4500:6000,6500:8000)
iter.include <- c(50:200,250:400,450:600,650:800)
iter.include <- 1:400

plot(NUTS_5.3a$logtau[iter.include],type='l')
plot(NUTS_5.3a$logitphi[iter.include],type='l')
plot(NUTS_5.3a$logd[iter.include],type='l')
plot(NUTS_5.3a$alphaU[iter.include],type='l')
plot(NUTS_5.3a$alphaR[iter.include],type='l')
plot(NUTS_5.3a$b[iter.include,1],type='l')
plot(NUTS_5.3a$b[iter.include,17],type='l')
plot(NUTS_5.3a$Yplus[iter.include,1],type='l')
plot(NUTS_5.3a$Yplus[iter.include,2],type='l')
plot(NUTS_5.3a$Yplus[iter.include,3],type='l')

each_chain <- nrow(NUTS_5.3a$reg_params)/4
par(mfrow=c(1,3))
for(k in 1:ncol(NUTS_5.3a$Yplus)){
  plot(cumsum(NUTS_5.3a$Yplus[1:each_chain,k])/(1:each_chain),type='l',col=1,ylab='',
       ylim=c(min(NUTS_5.3a$Yplus[,k]),max(NUTS_5.3a$Yplus[,k])),main = paste0('Running mean of Yplus_',k))
  for(c in 2:chains){
    lines(cumsum(NUTS_5.3a$Yplus[((c-1)*each_chain+1):(c*each_chain),k])/(1:each_chain),col=c)
  }
}


adm2_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_ZMB/admin2_u1_2011_urban_frac.rds"))
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Zambia/worldpop/adm2_weights_u1.rda"))
adm2_weights <- weight.adm2.u1[weight.adm2.u1$years==2011,]

strat.res.dat <- data_list$pop_strata_dat[order(data_list$pop_strata_dat$U,decreasing=T),]
strat.res.dat$mean <- colMeans(NUTS_5.3a$r[iter.include,],na.rm = T)
strat.res.dat$lower <- apply(NUTS_5.3a$r[iter.include,],2,quantile,probs=0.025)
strat.res.dat$upper <- apply(NUTS_5.3a$r[iter.include,],2,quantile,probs=0.975)
strat.res.dat$median <- colMedians(NUTS_5.3a$r[iter.include,],na.rm = T)

# fix later to get bounds
adm2.res.dat <- merge(strat.res.dat,adm2_UR_weights[,2:3],by.x='admin2.char',by.y='adm_idx') %>% mutate(urb_frac = ifelse(U==0,1-urb_frac,urb_frac)) %>% 
  group_by(admin2.char,admin1.char) %>% summarise(mean = sum(urb_frac*mean),median = sum(urb_frac*median))


ggplot() + geom_point(data=adm2.res.dat,aes(mean,admin1.char,col='2',pch='2'),size=4) +
  geom_point(data=admin1.dir,aes(admin1_mean,admin1.char,col='1',pch='1'),size=4) + # admin1 direct
  geom_point(data=admin1.dir,aes(admin1_lower,admin1.char,col='1'),pch=3,size=3) + # admin1 direct lower
  geom_point(data=admin1.dir,aes(admin1_upper,admin1.char,col='1'),pch=3,size=3) + # admin1 direct upper
  geom_point(data=admin1.res.inla,aes(adm1_median,admin1.char,col='3',pch='3'),size=3) + # admin1 FE INLA
  #geom_point(data=sd.est,aes(median,region)) +
 # geom_point(data=sd.est,aes(lower,region),pch=3) +
 # geom_point(data=sd.est,aes(upper,region),pch=3) +
 # geom_point(data=all_res,aes(adm1_median_agg,admin1.char,col='4',pch='4'),size=4) +
  ggtitle('Zambia 2009-2013') + theme_minimal() +
  ylab('') + xlab('NMR') + theme(legend.position = 'bottom') +
  scale_colour_manual(name = '', values =c('1'='blue','3'='red','2'='orange','4'='red'),
                      labels = c('Admin 1 Direct','Admin 2','Admin 2, Aggregated','Admin 2 INLA, Aggregated')) +
  scale_shape_manual(name = '', values =c('1'=17,'2'=1,'3'=2,'4'=17),
                     labels = c('Admin 1 Direct','Admin 2','Admin 2, Aggregated','Admin 2 INLA, Aggregated'))

