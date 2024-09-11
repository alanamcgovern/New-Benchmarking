library(sf)
library(SUMMER)
library(tidyverse)
library(spdep)

# setting_number <- 1
# 
# load(file=paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim_setting_",setting_number,".rda"))
# #load(file=paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim_setting_1.rda"))
# 
# all_dat <- sim_setting$all_dat
# n_clusters_urban <- (all_dat %>% filter(U==1) %>% group_by(A1) %>% summarise(n=n()))$n
# n_clusters_rural <- (all_dat %>% filter(U==0) %>% group_by(A1) %>% summarise(n=n()))$n

# choose some country to use as population surface: Angola ----
poly.adm1 <- read_sf(dsn = file.path("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_AGO_shp/gadm41_AGO_1.shp"))
poly.adm1$A1 <-1:nrow(poly.adm1)
n_admin1 <- nrow(poly.adm1)

poly.adm2 <- read_sf(dsn = file.path("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_AGO_shp/gadm41_AGO_2.shp"))
admin2.mat <- nb2mat(poly2nb(poly.adm2), zero.policy = TRUE)
colnames(admin2.mat) <- rownames(admin2.mat) <- paste0("admin2_", 1:dim(admin2.mat)[1])
poly.adm2$A2 <- 1:nrow(admin2.mat)
n_admin2 <- nrow(poly.adm2)

admin.key <- data.frame(NAME_1 = poly.adm2$NAME_1, A2 = poly.adm2$A2)
admin.key <- merge(admin.key,data.frame(NAME_1 = poly.adm1$NAME_1, A1 = poly.adm1$A1))
# adjust adjacency matrix -----

{
  admin.key$A1 <- ifelse(admin.key$NAME_1 %in% c('Zaire','Uíge'),1,
                                    ifelse(admin.key$NAME_1 %in% c('Bengo','Luanda'),2,
                                           ifelse(admin.key$NAME_1 %in% c('Cuanza Norte','Cuanza Sul'),3,
                                                  ifelse(admin.key$NAME_1 %in% c('Benguela','Namibe','Huíla'),4,
                                                         ifelse(admin.key$NAME_1 %in% c('Cunene','Cuando Cubango'),5,
                                                                ifelse(admin.key$NAME_1 %in% c('Huambo','Bié'),6,
                                                                       ifelse(admin.key$NAME_1 %in% c('Lunda Norte','Malanje'),7,
                                                                              ifelse(admin.key$NAME_1 %in% c('Lunda Sul','Moxico'),8, NA))))))))
  poly.adm2$A1 <- ifelse(poly.adm2$NAME_1 %in% c('Zaire','Uíge'),1,
                         ifelse(poly.adm2$NAME_1 %in% c('Bengo','Luanda'),2,
                                ifelse(poly.adm2$NAME_1 %in% c('Cuanza Norte','Cuanza Sul'),3,
                                       ifelse(poly.adm2$NAME_1 %in% c('Benguela','Namibe','Huíla'),4,
                                              ifelse(poly.adm2$NAME_1 %in% c('Cunene','Cuando Cubango'),5,
                                                     ifelse(poly.adm2$NAME_1 %in% c('Huambo','Bié'),6,
                                                            ifelse(poly.adm2$NAME_1 %in% c('Lunda Norte','Malanje'),7,
                                                                   ifelse(poly.adm2$NAME_1 %in% c('Lunda Sul','Moxico'),8, NA))))))))
  admin.key <- admin.key[!is.na(admin.key$A1),]
  n_admin1 <- length(unique(admin.key$A1))
}

{
  # n_admin1 <- 8
  # admin.key <- admin.key[admin.key$A1 %in% 1:n_admin1,]
}

# map plot
poly.adm2 %>% filter(!is.na(A1)) %>% ggplot() +geom_sf(aes(fill=factor(A1))) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.position = 'bottom') +
  scale_fill_discrete(name='First administrative area')

admin2.mat <- admin2.mat[admin.key$A2,admin.key$A2]
admin.key$admin1.char <- paste0('admin1_',admin.key$A1)
n_admin2 <- length(unique(admin.key$A2))
admin.key$A2 <- 1:nrow(admin.key)
admin.key$admin2.char <- paste0('admin2_',admin.key$A2)
colnames(admin2.mat) <- rownames(admin2.mat) <- admin.key$admin2.char

admin2.mat.nested <- admin2.mat
for(i in 1:nrow(admin2.mat.nested)){
  admin2.mat.nested[i,which(admin.key$A1!=admin.key$A1[i])] <- 0
  if(sum(admin2.mat.nested[i,])>0){
    admin2.mat.nested[i,] <- admin2.mat.nested[i,]/sum(admin2.mat.nested[i,])
  }
}

Q.admin2 <- -admin2.mat
Q.admin2 <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})*Q.admin2
diag(Q.admin2) <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})
diag(Q.admin2)[diag(Q.admin2)==0] <- 1
Q_scaled <- INLA::inla.scale.model(Q.admin2, constr=list(A=matrix(1,nrow=1,ncol=n_admin2), e=0))
Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))

Q.admin2.nested <- -admin2.mat.nested
Q.admin2.nested <- sapply(1:nrow(Q.admin2.nested),function(i){sum(I(Q.admin2.nested[i,]!=0))})*Q.admin2.nested
diag(Q.admin2.nested) <- sapply(1:nrow(Q.admin2.nested),function(i){sum(I(Q.admin2.nested[i,]!=0))})
constraints <- matrix(0,nrow=n_admin1,ncol=n_admin2) #n_admin1 rows, n_admin2 columns
for(a in 1:n_admin2){constraints[admin.key$A1[admin.key$A2==a],a] <- 1}
Q_nested_scaled <- INLA::inla.scale.model(Q.admin2.nested, constr=list(A=constraints, e=rep(0,n_admin1)))
Q_nested_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_nested_scaled))

# population surface parameters -----

n_clusters_urban <- pmax(10,round(rnorm(n_admin1,700,100))) #number of urban clusters in each admin1 area
n_clusters_rural <- pmax(10,round(rnorm(n_admin1,800,100))) #number of rural clusters in each admin1 area

n_births_urban <- 100 #average number of births per urban cluster
sd_n_births_urban <- 10
n_births_rural <-  125 #average number of births per rural cluster
sd_n_births_rural <- 10

# simulate population surface -----

all_dat <- data.frame(
  # (true) number of births in each cluster
  N = round(c(rnorm(sum(n_clusters_urban),n_births_urban,sd_n_births_urban),
              rnorm(sum(n_clusters_rural),n_births_rural,sd_n_births_rural))),
  # urban or rural strata (U=1 if urban, otherwise 0)
  U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
  # admin area of each cluster
  A1 = c(unlist(sapply(1:n_admin1,function(i){
    rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_admin1,function(i){rep(i,n_clusters_rural[i])}))))
#assign clusters
all_dat$cluster <- 1:nrow(all_dat)
#assign admin2 -- CAN CHANGE TO VARY SIZES OF ADMIN2 AREAS
all_dat$A2 <- sapply(1:nrow(all_dat),function(x){
  adm2.areas.tmp <- admin.key[admin.key$A1==all_dat$A1[x],]$A2
  sample(adm2.areas.tmp,1,prob = rep(1,length(adm2.areas.tmp)))
})

# sampling parameters ------
n_clusters_urban_samp <- round(0.05*n_clusters_urban) #round(0.1*n_clusters_urban) #number of urban clusters sampled in each admin2 area
n_clusters_rural_samp <- round(0.03*n_clusters_rural) #round(0.08*n_clusters_rural) #number of rural clusters sampled in each admin2 area
n_births_urban_samp <- 15 # average number of births sampled from each urban cluster 
sd_births_urban_samp <- 2
n_births_rural_samp <- 20 # average number of births sampled from each rural cluster
sd_n_births_rural_samp <- 2

# calculate sampling and population weights -----
pop_dat <- all_dat %>% group_by(A1,U) %>% summarise(N=sum(N))
pop_dat$num_clust_samp <- pop_dat$num_births_samp <- NA
pop_dat$num_clust_samp[pop_dat$U==0] <- n_clusters_rural_samp
pop_dat$num_clust_samp[pop_dat$U==1] <- n_clusters_urban_samp
pop_dat$num_births_samp[pop_dat$U==0] <- n_births_rural_samp 
pop_dat$num_births_samp[pop_dat$U==1] <- n_births_urban_samp
pop_dat$wt <- pop_dat$N/(pop_dat$num_births_samp*pop_dat$num_clust_samp)
pop_dat$wt <- pop_dat$wt/100
admin1_UR_weights <- all_dat %>% group_by(A1) %>% summarise(urban_prop = sum(U*N)/sum(N))
admin2_UR_weights <- all_dat %>% group_by(A2) %>% summarise(urban_prop = sum(U*N)/sum(N))
admin2_weights <- all_dat %>% group_by(A2,A1) %>% summarise(prop = sum(N)/sum(all_dat$N))
admin2_to_admin1_weights <- admin2_weights %>% group_by(A1) %>% mutate(prop = prop/sum(prop))

# risk surface parameters -----
setting_number <- 4
d <- 0.25
sigma <- c(0.15,0.05,0.05,0.15)[setting_number] #for BYM2
phi <- c(0.25,0.25,0.7,0.25)[setting_number]
alphaU <- log(0.02)
alphaR <- log(0.025)
beta <- c(-0.025,0.025,-0.05,0.05,-0.1,0.1,-.2,0.5)

#plot(exp(admin1_UR_weights$urban_prop*betaU + (1-admin1_UR_weights$urban_prop)*betaR))

# simulate  risk surface ------
b <- as.vector(Rfast::rmvnorm(1,rep(0,n_admin2),sigma^2*(diag((1-phi),n_admin2) + phi*Q_scaled_inv)))

#draw deaths
all_dat$Y <- sapply(1:nrow(all_dat),function(i){
  rnbinom(1,
          size = 1/d*all_dat$N[i]*exp(alphaU*all_dat$U[i] + alphaR*(1-all_dat$U[i]) + beta[all_dat$A1[i]] + b[all_dat$A2[i]]),
          prob = 1/(1+d))})

# simulations ------

nsim <- 500
sim_data <- list()
admin1.dir.res <- NULL
#admin1.res <- admin2.res <- NULL
start.time <- Sys.time()
for(k in 1:nsim){
  if(k%%10==0){
    message(paste0('starting simulation ', k))
  }
  # simulate data ------
  ### sample clusters
  obs_dat <- NULL
  for(area in 1:n_admin1){
    # randomly select clusters from urban strata based on size of cluster
    urban_strata_dat <- all_dat[all_dat$A1==area & all_dat$U==1,]
    rural_strata_dat <- all_dat[all_dat$A1==area & all_dat$U==0,]
    obs_dat <- rbind(obs_dat,urban_strata_dat[sample(1:n_clusters_urban[area],n_clusters_urban_samp[area],prob = urban_strata_dat$N),],
                     # randomly select cluster from rural strata  based on size of cluster
                     rural_strata_dat[sample(1:n_clusters_rural[area],n_clusters_rural_samp[area],prob = rural_strata_dat$N),])
  }
  
  # sample births
  obs_dat$n <- obs_dat$Z <- NA
  obs_dat[obs_dat$U==1,]$n <- round(rnorm(sum(obs_dat$U),n_births_urban_samp,sd_births_urban_samp)) 
  obs_dat[obs_dat$U==0,]$n <- round(rnorm(sum(1-obs_dat$U),n_births_rural_samp,sd_n_births_rural_samp))
  obs_dat[obs_dat$n>obs_dat$N,]$n <- obs_dat[obs_dat$n>obs_dat$N,]$N # make sure it less than total births
  if(length(obs_dat[obs_dat$n<1,]$n)>0)
    obs_dat[obs_dat$n<1,]$n <- 1 #make sure its larger than 0
  for(i in 1:nrow(obs_dat)){
    obs_dat$Z[i] <- sum(sample(c(rep(1,obs_dat$Y[i]),rep(0,(obs_dat$N[i]-obs_dat$Y[i]))), obs_dat$n[i]))
  }
  
  # add additional variables (sampling weights and things for SUMMER package)
  obs_dat <- merge(obs_dat,pop_dat[c('A1','U','wt')])
  #obs_dat <- obs_dat[order(obs_dat$A1),]
  obs_dat$admin1.char <- paste0('admin1_',obs_dat$A1)
  obs_dat$years <- 'All'
  obs_dat$age <- '0'
  obs_dat$died <- obs_dat$Z
  obs_dat$strata <- paste0(obs_dat$A1,'.',c('rural','urban')[obs_dat$U+1])
  
  # get direct estimates -----
  suppressMessages({
    admin1.dir <- SUMMER::getDirect(obs_dat, years='All',
                                 regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                                 Ntrials = "n", weightsVar = "wt", CI=0.9)
  })
  admin1.dir <- admin1.dir[admin1.dir$region!='All',]
  admin1.dir$admin1 <- as.numeric(str_remove(admin1.dir$region,'admin1_'))
  admin1.dir <- admin1.dir[order(admin1.dir$admin1),]
  
  dir.est <- admin1.dir %>% dplyr::select(admin1,mean,lower,upper) %>% 
    rename(dir.mean=mean,dir.lower=lower,dir.upper=upper)
  dir.est$dir.sd <- (dir.est$dir.upper - dir.est$dir.lower)/(2*1.645)
  
  # fit INLA model: admin1 fixed effect + admin2 BYM2 ------
  admin2.inla.fit <- INLA::inla(Z ~ factor(U) + factor(A1) + 
                                 f(A2, model="bym2", graph = admin2.mat.nested, scale.model=T, constr=T,
                                    adjust.for.con.comp=T,hyper=list(phi=list(prior="pc", param=c(0.5, 2/3)),
                                                                     prec=list(prior="pc.prec",param=c(1,0.01)))) -1,
                                family = 'nbinomial', 
                                data = obs_dat, 
                                control.family=list(link='log'),
                                control.fixed = list(mean =list("factor(U)0"=-3.5, "factor(U)1"=-3.5, default = 0),
                                                     prec=list(default = 1/9)),
                                control.predictor = list(compute = FALSE, link = 1), 
                                control.compute = list(config=T),
                                E = n)
  
  # sample from posterior
  cs <- admin2.inla.fit$misc$configs$contents$tag
  cs <- cs[cs != "Predictor"]
  select <- list()
  for (i in 1:length(cs)) {
    select[[i]] <- 0
    names(select)[i] <- cs[i]
  }
  sampFull <- INLA::inla.posterior.sample(n = 1000, result = admin2.inla.fit, intern = TRUE, selection = select)
  
  # save data for benchmarking model -----
  #get estimate of covaraiance matrix (to use as inverse of scaling matrix)
  sampFull.draws <- matrix(NA,1000,4+n_admin1+n_admin2)
  for(i in 1:1000){
    sampFull.draws[i,] <- c(sampFull[[i]]$hyperpar[2:3],
                            -sampFull[[i]]$hyperpar[1], #kind of comparable to log(d)
                            sampFull[[i]]$latent[(2*n_admin2+2)], # urban intercept
                            sampFull[[i]]$latent[(2*n_admin2+1)], # rural intercept
                            tail(sampFull[[i]]$latent,n_admin1-1),
                            sampFull[[i]]$latent[c(1:n_admin2)]) # admin2 spatial effect
  }
  Sigma <- cov(sampFull.draws)
  
  pop_strata_dat <- all_dat %>% group_by(A2,U,A1) %>% summarise(N=sum(N)) 
  
  sim_data[[k]] <- list(M_init = solve(Sigma),
                    logit_r_hat = admin1.dir$logit.est,
                    logit_r_hat_var = admin1.dir$var.est,
                    obs_dat = obs_dat,
                    pop_strata_dat = as.data.frame(pop_strata_dat),
                    Q_scaled_inv=Q_nested_scaled_inv)
  
  # get INLA estimates --------
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
  #   adm1_area_tmp <- admin.key[admin.key$A2==area,]$A1
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
  # 
  # adm2_wt_mat <- matrix(0,2*n_admin2,n_admin2)
  # for(area in 1:n_admin2){
  #   adm2_wt_mat[area,area] <- admin2_UR_weights[admin2_UR_weights$A2==area,]$urban_prop #urban
  #   adm2_wt_mat[area+n_admin2,area] <- 1-admin2_UR_weights[admin2_UR_weights$A2==area,]$urban_prop # rural
  # }
  # admin2_samples <- admin2_strat_samples%*%adm2_wt_mat
  # 
  # admin2.inla <- data.frame(A2 = 1:n_admin2, admin2.char = paste0('admin2_',1:n_admin2),
  #                           inla.mean = apply(admin2_samples,2,mean),
  #                           inla.median = apply(admin2_samples,2,median),
  #                           inla.sd = apply(admin2_samples,2,sd),
  #                           inla.lower = apply(admin2_samples,2,quantile,prob = 0.05),
  #                           inla.upper = apply(admin2_samples,2,quantile,prob = 0.95))
  # 
  # adm1_wt_mat <- matrix(0,n_admin2,n_admin1)
  # for(area in 1:n_admin1){
  #   adm2_areas_tmp <- admin.key[admin.key$A1==area,]$A2
  #   adm1_wt_mat[adm2_areas_tmp,area] <- admin2_to_admin1_weights[admin2_to_admin1_weights$A2 %in% adm2_areas_tmp,]$prop
  # }
  # admin1_samples <- admin2_samples%*%adm1_wt_mat
  # admin1.inla <- data.frame(admin1 = 1:n_admin1, admin1.char = paste0('admin1_',1:n_admin1),
  #                           inla.mean = apply(admin1_samples,2,mean),
  #                           inla.median = apply(admin1_samples,2,median),
  #                           inla.sd = apply(admin1_samples,2,sd),
  #                           inla.lower = apply(admin1_samples,2,quantile,prob = 0.05),
  #                           inla.upper = apply(admin1_samples,2,quantile,prob = 0.95))
  
  # what do we want to record? ----
  
  dir.est$iter.id <- k
  admin1.dir.res <- rbind(admin1.dir.res,dir.est)
  
  # admin1.res.tmp <- merge(admin1.inla,dir.est)
  # data.info <- obs_dat %>% group_by(A1) %>% reframe(n_births = sum(n), n_deaths = sum(Z), n_clusters = n(),
  #                                                   urban_births = sum(U*n), urban_deaths = sum(U*Z), urban_clusters = sum(U))
  # admin1.res.tmp <- merge(admin1.res.tmp,data.info,by.x='admin1',by.y='A1')
  # admin1.res.tmp$iter.id <- k
  # admin1.res <- rbind(admin1.res,admin1.res.tmp)
  #
  # data.info.admin2 <- obs_dat %>% group_by(A2) %>% reframe(n_births = sum(n), n_deaths = sum(Z), n_clusters = n())
  # missing.admin2 <- (1:n_admin2)[!(1:n_admin2 %in% data.info.admin2$A2)]
  # if(length(missing.admin2)>0){
  #   data.info.admin2 <- rbind(data.info.admin2, data.frame(A2 = missing.admin2, 
  #                                                          n_births = 0, n_deaths = NA, n_clusters = 0))
  # }
  # admin2.res.tmp <- merge(admin2.inla,data.info.admin2)
  # admin2.res.tmp$iter.id <- k
  # admin2.res <- rbind(admin2.res, admin2.res.tmp)
  
}
Sys.time() - start.time

# saving simulations ------

# save simulation settings
sim_setting <- list(all_dat=all_dat,
                    index = setting_number,
                    alphaU=alphaU, alphaR=alphaR, beta=beta, 
                    b=b, sigma=sigma, phi=phi, d=d)
save(sim_setting,file=paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim_setting_",setting_number,".rda"))

# save all simulated datasets (for cluster)
save(sim_data,file=paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_number,"_datasets.rda"))

# save all direct/INLA results
save(admin1.dir.res,file=paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_number,"_Dir.rda"))

# # save all direct/INLA results
# save(admin1.res,file=paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_number,"_DirINLA_admin1.rda"))
# save(admin2.res,file=paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_number,"_DirINLA_admin2.rda"))

# get summary of datasets ------

load(file=paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_number,"_datasets.rda"))

simdata_summary <- matrix(NA,500,5)
for(i in 1:500){
  dat_tmp <- sim_data[[i]]$obs_dat[,c('A1','A2','n','Z')]
  # add rows for missing admin2
  missing.adm2 <- which(!(1:n_admin2 %in% unique(dat_tmp$A2)))
  if(length(missing.adm2)>0)
    dat_tmp <- rbind(dat_tmp,data.frame(A1=1,A2=missing.adm2,n=0,Z=0))
  
  adm1_tmp <- dat_tmp %>% group_by(A1) %>% reframe(sum_Z=sum(Z),sum_n=sum(n)) %>% 
    reframe(avg_births_adm1 = mean(sum_n),avg_deaths_adm1 = mean(sum_Z))
  
  adm2_tmp <- dat_tmp %>% group_by(A2) %>% reframe(sum_Z=sum(Z),sum_n=sum(n)) %>% 
    reframe(avg_births_adm2 = mean(sum_n),avg_deaths_adm2 = mean(sum_Z),no_deaths_adm2 = mean(sum_Z==0))
  
  simdata_summary[i,] <- c(adm1_tmp$avg_births_adm1,adm1_tmp$avg_deaths_adm1,
                           adm2_tmp$avg_births_adm2,adm2_tmp$avg_deaths_adm2,adm2_tmp$no_deaths_adm2)
}
colnames(simdata_summary) <- c('avg_births_adm1','avg_deaths_adm1','avg_births_adm2','avg_deaths_adm2','no_deaths_adm2')
colMeans(simdata_summary)
