library(tidyverse)
library(locfit)
library(rdhs)
library(surveyPrev)
library(SUMMER)
library(labelled)
library(sf)
#library(ggpubr)

load(file ='/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Info/Kenya_general_info.Rdata') # load the country info

# load shapefiles -----
poly.path <- "/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_KEN_shp"
poly.adm1 <- st_read(dsn = poly.path, layer = poly.layer.adm1, options = "ENCODING=UTF-8")

geo <- getDHSgeo(country = 'Kenya', year = 2022)

admin1.info <- adminInfo(geo = poly.adm1, admin = 1)
admin1.info$admin.info$admin1 <- 1:nrow(admin1.info$admin.info)
admin1.info$admin.info$admin1.char <- paste0('admin1_',admin1.info$admin.info$admin1)
row.names(admin1.info$admin.mat) <- colnames(admin1.info$admin.mat) <- admin1.info$admin.info$admin1
# re-structure admin1 matrix for INLA
admin1.mat <- matrix(0,ncol=nrow(admin1.info$admin.mat),nrow = nrow(admin1.info$admin.mat))
for(i in 1:nrow(admin1.mat)){
  admin1.mat[i,] <- -1*I(admin1.info$admin.mat[i,]!=0)
}
diag(admin1.mat) <- -colSums(admin1.mat)
save(admin1.mat, file='/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Kenya/Kenya_Amat.rda')
load('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Kenya/Kenya_Amat.rda')

# load data ------
load(paste0('/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Kenya/Kenya_cluster_dat.rda')) 

adm1_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_ken/admin1_u1_2020_urban_frac.rds"))
adm1_UR_weights <- adm1_UR_weights[order(as.numeric(str_remove(adm1_UR_weights$adm_idx,'admin1_'))),]

# admin1 weights 
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Kenya/worldpop/adm1_weights_u1.rda"))
adm1_weights <- weight.adm1.u1[weight.adm1.u1$years==2020,]

# direct estimates -----
data_for_direct <- nmr.dat
data_for_direct$age <- 0
data_for_direct$years <- 'All'
data_for_direct$died <- data_for_direct$Y
dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                             regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                             Ntrials = "total",weightsVar = "v005")

natl.dir <- dir.est[dir.est$region=='All',]
admin1.dir <- dir.est[dir.est$region!='All',]
admin1.dir <- admin1.dir %>% dplyr::select(region,mean,lower,upper) %>% rename(admin1_mean = mean, admin1_lower = lower,admin1_upper=upper) %>% rename(admin1.char = region)
admin1.dir$admin1 <- as.numeric(sapply(1:nrow(admin1.dir),function(i){str_split(admin1.dir$admin1.char[i],'_')[[1]][2]}))

# admin 1 random effect model-----
admin1RE.fit <- INLA::inla(Y ~ urban -1 +  f(admin1, model="bym2", graph = admin1.mat, scale.model=T, constr=T,
                                             adjust.for.con.comp=T, 
                                             hyper=list(phi=list(prior="pc", param=c(0.5, 0.5)),
                                                        prec=list(prior='pc.prec', param=c(1, 0.01)))),
                           family = 'binomial', 
                           data = nmr.dat, 
                           control.family=list(variant=1,link='log'),
                           control.predictor = list(compute = FALSE, link = 1), 
                           control.compute = list(config=T),
                           Ntrials = nmr.dat$total)

# sample from posterior
cs <- admin1RE.fit$misc$configs$contents$tag
cs <- cs[cs != "Predictor"]
select <- list()
for (i in 1:length(cs)) {
  select[[i]] <- 0
  names(select)[i] <- cs[i]
}
sampFull <- INLA::inla.posterior.sample(n = 1000, result = admin1RE.fit, intern = TRUE, selection = select)
sampFull.draws <- matrix(NA,1000,length(sampFull[[1]]$latent))
for(i in 1:1000){
  sampFull.draws[i,] <- sampFull[[i]]$latent
}
fields <- colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)

admin1RE.inla <- matrix(0,ncol=(2*max(nmr.dat$admin1)+2),nrow=max(nmr.dat$admin1))
for(area in 1:max(nmr.dat$admin1)){
  admin1RE.inla[area,area] <- 1 #choose admin1 area
}
admin1RE.inla[,2*max(nmr.dat$admin1)+1] <- adm1_UR_weights$urb_frac
admin1RE.inla[,2*max(nmr.dat$admin1)+2] <- 1 - adm1_UR_weights$urb_frac

admin1RE.inla <- admin1RE.inla%*%t(sampFull.draws)
res.admin1RE <- data.frame(admin1 = unique(nmr.dat$admin1),admin1.char = unique(nmr.dat$admin1.char))

res.admin1RE$adm1_median <- (apply(admin1RE.inla,1,function(x){median(exp(x))}))
res.admin1RE$adm1_lower95 <- (apply(admin1RE.inla,1,function(x){quantile(exp(x),0.025)}))
res.admin1RE$adm1_upper95 <- (apply(admin1RE.inla,1,function(x){quantile(exp(x),0.975)}))

natl.inla <- rep(0,2*max(nmr.dat$admin1)+2)
for(area in 1:max(nmr.dat$admin1)){
  natl.inla[area] <- weight.adm1.u1$proportion[area] #choose admin1 area
}
natl.inla[2*max(nmr.dat$admin1)+1] <- sum(adm1_UR_weights$urb_frac)/nrow(adm1_UR_weights)
natl.inla[2*max(nmr.dat$admin1)+2] <- 1 - sum(adm1_UR_weights$urb_frac)/nrow(adm1_UR_weights)
natl.inla <- natl.inla%*%t(sampFull.draws)

res.natl <- data.frame(median = median(exp(natl.inla)),
                       lower = quantile(exp(natl.inla),0.025),
                       upper = quantile(exp(natl.inla),0.975))

#plot -----
plot(res.admin1RE$adm1_median)
abline(h=natl.dir$mean,col='red')
abline(h=natl.dir$lower,col='red',lty=2)
abline(h=natl.dir$upper,col='red',lty=2)
abline(h=res.natl$median,col='blue')
abline(h=res.natl$lower,col='blue',lty=2)
abline(h=res.natl$upper,col='blue',lty=2)

# get generalized inverse of scaled Q for BYM2 ----
Q.admin1 <- admin1.mat
Q_scaled <- INLA::inla.scale.model(Q.admin1, constr=list(A=t(eigen(Q.admin1)$vectors[,eigen(Q.admin1)$values<1e-10]), e=rep(0,sum(eigen(Q.admin1)$values<1e-10))))
Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))

# input strata level DHS survey info ----
dhs_dat <- data.frame(admin1 = rep(sort(unique(admin.key$admin1)),each=2),U = rep(c(1,0),max(admin.key$admin1)))
dhs_dat <- left_join(dhs_dat,unique.array(admin.key[,c(1:3)])) %>% rename(A1=admin1)
# total number of HH in frame
dhs_dat$M 
# total number of clusters in strata 
dhs_dat$n_clusters 
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
  dplyr::select(cluster,total,Y,U,admin1,admin1.char,admin1.name) %>% 
  rename(Z=Y,n=total,A1=admin1)
# add estimated number of total births per cluster
alg_dat <- merge(alg_dat,frame.info,by=c('admin1.char','A1','U')) %>% dplyr::select(-M_bar,-n_clusters,-n_clusters_samp) %>%
  mutate(N=round(M_bar_est*n/n_HH_samp))

# use worldpop U1 and then rescale using admin1 births from census
births_dat <- data.frame(admin1.char = paste0('admin1_',sort(unique(admin.key$admin1))),
                         births = 5*c())
pop_fracs <- NULL 
for(year in 2018:2022){
  adm1_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_KEN/admin1_u1_",year,"_urban_frac.rds"))
  adm1_UR_weights <- adm1_UR_weights[order(as.numeric(str_remove(adm1_UR_weights$adm_idx,'admin1_'))),]
  pop_fracs_tmp <- data.frame(admin1.char = rep(admin.key$admin1.char,each=2),U=rep(c(0,1),max(admin.key$admin1)))
  pop_fracs_tmp <- merge(pop_fracs_tmp,weight.adm1.u1[weight.adm1.u1$year==year,c(1,4)])
  pop_fracs_tmp <- merge(pop_fracs_tmp,adm1_UR_weights[,2:3],by.y='adm_idx',by.x='admin1.char')
  pop_fracs_tmp <- pop_fracs_tmp %>% mutate(full_prop=ifelse(U==1,proportion*urb_frac,proportion*(1-urb_frac)))
  pop_fracs_tmp$year <- year
  
  pop_fracs <- rbind(pop_fracs,pop_fracs_tmp[,c(1,2,5,6)])
}
pop_fracs <- merge(admin.key,pop_fracs)

pop_fracs$N <- pop_dat_est$full_prop*births_dat$births/sum(pop_dat_est$full_prop)

pop_strata_dat <- pop_fracs %>% group_by(admin1.char,U) %>% summarise(N = round(sum(N)))
pop_strata_dat <- merge(admin.key,pop_strata_dat) %>% rename(A1=admin1)
pop_strata_dat <- pop_strata_dat[order(pop_strata_dat$A1,pop_strata_dat$U),]

# check if these match up
alg_dat %>% filter(admin1.char %in% c(pop_strata_dat[pop_strata_dat$N==0& pop_strata_dat$U==1,]$admin1.char) & U==1)

data_list <- list(obs_dat = alg_dat,
                  pop_strata_dat = as.data.frame(pop_strata_dat),
                  Q_scaled_inv = Q_scaled_inv)


