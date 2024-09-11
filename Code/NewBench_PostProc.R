## process results from cluster and combined with direct and INLA estimates

library(tidyverse)

# simulation setting number
setting_num <- 3
# number of simulated datasets
nsim <- 500

load(file = paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim_setting_",setting_num,".rda"))

# calculate population weights (same for all datasets simulated under the same seting)
admin1_UR_weights <- sim_setting$all_dat %>% group_by(A1) %>% summarise(urban_prop = sum(U*N)/sum(N))
admin2_UR_weights <- sim_setting$all_dat %>% group_by(A2) %>% summarise(urban_prop = sum(U*N)/sum(N))
admin2_weights <- sim_setting$all_dat %>% group_by(A2,A1) %>% summarise(adm1.wt = sum(N)/sum(sim_setting$all_dat$N))

admin.key <- admin2_weights[,c('A1','A2')]
n_admin1 <- max(admin.key$A1)
n_admin2 <- max(admin.key$A2)

# Process Benchmarked estimates ============
setwd(paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/Sim",setting_num))
admin1.res.bench <- admin2.res.bench <- NULL
for(i in 1:nsim){
  print(i)
  
  raw.res <- readRDS(paste0('Results',i,'.rds'))
  
  ## summarise estimates for simulation (admin2 and admin1 rate means, medians, SDs, quantiles)
  raw.res <- raw.res[,(n_admin2+n_admin1+6):(3*n_admin2+n_admin1+5)]
  
  adm2_wt_mat <- matrix(0,2*n_admin2,n_admin2)
  for(k in 1:n_admin2){
    adm2_wt_mat[k,k] <- admin2_UR_weights[admin2_UR_weights$A2==k,]$urban_prop #urban
    adm2_wt_mat[k+n_admin2,k] <- 1-admin2_UR_weights[admin2_UR_weights$A2==k,]$urban_prop # rural
  }
  admin2.samples <- as.matrix(raw.res)%*%adm2_wt_mat
  
  admin2.res.tmp <- data.frame(A2 = 1:n_admin2, iter.id = i,
                               bench.mean = apply(admin2.samples,2,mean),
                               bench.median = apply(admin2.samples,2,median),
                               bench.sd = apply(admin2.samples,2,sd),
                               bench.lower = apply(admin2.samples,2,quantile,prob = 0.05),
                               bench.upper = apply(admin2.samples,2,quantile,prob = 0.95))
  admin2.res.bench <- rbind(admin2.res.bench,admin2.res.tmp)
  
  adm1_wt_mat <- matrix(0,n_admin2,n_admin1)
  for(k in 1:n_admin1){
    adm2_areas_tmp <- admin.key[admin.key$A1==k,]$A2
    adm1_wt_mat[adm2_areas_tmp,k] <- admin2_weights[admin2_weights$A2 %in% adm2_areas_tmp,]$adm1.wt/sum(admin2_weights[admin2_weights$A2 %in% adm2_areas_tmp,]$adm1.wt)
  }
  admin1.samples <- admin2.samples%*%adm1_wt_mat
  admin1.res.tmp <- data.frame(admin1 = 1:n_admin1, iter.id = i,
                               bench.mean = apply(admin1.samples,2,mean),
                               bench.median = apply(admin1.samples,2,median),
                               bench.sd = apply(admin1.samples,2,sd),
                               bench.lower = apply(admin1.samples,2,quantile,prob = 0.05),
                               bench.upper = apply(admin1.samples,2,quantile,prob = 0.95))
  admin1.res.bench <- rbind(admin1.res.bench,admin1.res.tmp)
}

# merge with direct estimates and save
load(paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_num,"_Dir.rda"))

admin1.res <- merge(admin1.dir.res, admin1.res.bench)
admin2.res <- admin2.res.bench
admin2.res <- merge(admin.key,admin2.res)

# add population NMR
admin1.res <- merge(admin1.res,sim_setting$all_dat %>% group_by(A1) %>% summarise(pop_nmr=sum(Y)/sum(N)),by.x='admin1',by.y='A1')
admin2.res <- merge(admin2.res,sim_setting$all_dat %>% group_by(A2) %>% summarise(pop_nmr=sum(Y)/sum(N)))

save(admin1.res, file = paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_num,"results_admin1.rda"))
save(admin2.res, file = paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_num,"results_admin2.rda"))


## Add Stan estimates =============
setwd(paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking"))
for(setting_number in 1:4){
  
  admin1.res.stan <- readRDS(file=paste0('Simulation_Study/Stan_Admin1Summary_',setting_number,'.rds'))
  admin2.res.stan <- readRDS(file=paste0('Simulation_Study/Stan_Admin2Summary_',setting_number,'.rds'))
  
  load(paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_number,"results_admin1.rda"))
  load(paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_number,"results_admin2.rda"))
  
  admin1.res <- merge(admin1.res,admin1.res.stan,by=c('admin1','iter.id'))
  admin2.res <- merge(admin2.res,admin2.res.stan,by=c('A2','iter.id'))
  
  save(admin1.res, file = paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_number,"results+stan_admin1.rda"))
  save(admin2.res, file = paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_number,"results+stan_admin2.rda"))
  
}

## Process Exact benchmarking estimates ===========
setwd(paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/Sim",setting_num,"_Exact"))
admin1.res.exact <- admin2.res.exact <- NULL
for(i in 1:nsim){
  print(i)
  
  raw.res <- readRDS(paste0('Results',i,'.rds'))
  
  ## summarise estimates for simulation (admin2 and admin1 rate means, medians, SDs, quantiles)
  raw.res <- raw.res[,(n_admin2+n_admin1+6):(3*n_admin2+n_admin1+5)]
  
  adm2_wt_mat <- matrix(0,2*n_admin2,n_admin2)
  for(k in 1:n_admin2){
    adm2_wt_mat[k,k] <- admin2_UR_weights[admin2_UR_weights$A2==k,]$urban_prop #urban
    adm2_wt_mat[k+n_admin2,k] <- 1-admin2_UR_weights[admin2_UR_weights$A2==k,]$urban_prop # rural
  }
  admin2.samples <- as.matrix(raw.res)%*%adm2_wt_mat
  
  admin2.res.tmp <- data.frame(A2 = 1:n_admin2, iter.id = i,
                               exact.mean = apply(admin2.samples,2,mean),
                               exact.median = apply(admin2.samples,2,median),
                               exact.sd = apply(admin2.samples,2,sd),
                               exact.lower = apply(admin2.samples,2,quantile,prob = 0.05),
                               exact.upper = apply(admin2.samples,2,quantile,prob = 0.95))
  admin2.res.exact <- rbind(admin2.res.exact,admin2.res.tmp)
  
  adm1_wt_mat <- matrix(0,n_admin2,n_admin1)
  for(k in 1:n_admin1){
    adm2_areas_tmp <- admin.key[admin.key$A1==k,]$A2
    adm1_wt_mat[adm2_areas_tmp,k] <- admin2_weights[admin2_weights$A2 %in% adm2_areas_tmp,]$adm1.wt/sum(admin2_weights[admin2_weights$A2 %in% adm2_areas_tmp,]$adm1.wt)
  }
  admin1.samples <- admin2.samples%*%adm1_wt_mat
  admin1.res.tmp <- data.frame(admin1 = 1:n_admin1, iter.id = i,
                               exact.mean = apply(admin1.samples,2,mean),
                               exact.median = apply(admin1.samples,2,median),
                               exact.sd = apply(admin1.samples,2,sd),
                               exact.lower = apply(admin1.samples,2,quantile,prob = 0.05),
                               exact.upper = apply(admin1.samples,2,quantile,prob = 0.95))
  admin1.res.exact <- rbind(admin1.res.exact,admin1.res.tmp)
}

# merge with direct, Stan, and Benchmarked estimates and save
load(paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_num,"results+stan_admin1.rda"))
load(paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_num,"results+stan_admin2.rda"))

admin1.res <- merge(admin1.res, admin1.res.exact)
admin2.res <- merge(admin2.res, admin2.res.exact) 
admin2.res <- merge(admin.key,admin2.res)

save(admin1.res, file = paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_num,"results_all_admin1.rda"))
save(admin2.res, file = paste0("/Users/alanamcgovern/Desktop/Research/New_Benchmarking/Simulation_Study/sim",setting_num,"results_all_admin2.rda"))
