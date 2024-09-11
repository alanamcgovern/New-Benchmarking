suppressMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly=TRUE)
total_sims <- as.numeric(args[1])
setting_number <-  as.numeric(args[2])

load(paste0('sim_setting_',setting_number,'.rda'))
n_admin1 <- length(unique(sim_setting$all_dat$A1))
n_admin2 <- length(unique(sim_setting$all_dat$A2))
admin.key <- sim_setting$all_dat %>% dplyr::select(A1,A2) %>% unique() %>% arrange(A2)
admin2_UR_weights <- sim_setting$all_dat %>% group_by(A2) %>% summarise(urban_prop = sum(U*N)/sum(N))
admin2_weights <- sim_setting$all_dat %>% group_by(A2,A1) %>% summarise(adm1.wt = sum(N)/sum(sim_setting$all_dat$N))

admin2.res.stan <- admin1.res.stan <- NULL
for(k in 1:total_sims){
  out <- readRDS(file=paste0('Sim',setting_number,'_stan/Results',k,'.rds'))
  
  suppressWarnings({
    eta_postsamp <- matrix(NA,nrow(out$alpha),2*n_admin2)
  })
  
  for(i in 1:nrow(admin.key)){
    if(admin.key$A1[i]>1){
      eta_postsamp[,i] <- out$alpha[,1] + out$beta[,admin.key$A1[i]-1] + out$b[,admin.key$A2[i]]
      eta_postsamp[,n_admin2 + i] <- out$alpha[,2] + out$beta[,admin.key$A1[i]-1] + out$b[,admin.key$A2[i]]
    }else{
      eta_postsamp[,i] <- out$alpha[,1] + out$b[,admin.key$A2[i]]
      eta_postsamp[,n_admin2 + i] <- out$alpha[,2] + out$b[,admin.key$A2[i]]
    }
  }
  r_postsamp <- exp(eta_postsamp)
  
  adm2_wt_mat <- matrix(0,2*n_admin2,n_admin2)
  for(area in 1:n_admin2){
    adm2_wt_mat[area,area] <- admin2_UR_weights[admin2_UR_weights$A2==area,]$urban_prop #urban
    adm2_wt_mat[area+n_admin2,area] <- 1-admin2_UR_weights[admin2_UR_weights$A2==area,]$urban_prop # rural
  }
  admin2.samples <- r_postsamp%*%adm2_wt_mat
  
  admin2.res.tmp <- data.frame(A2 = 1:n_admin2, iter.id = k,
                               stan.mean = apply(admin2.samples,2,mean),
                               stan.median = apply(admin2.samples,2,median),
                               stan.sd = apply(admin2.samples,2,sd),
                               stan.lower = apply(admin2.samples,2,quantile,prob = 0.05),
                               stan.upper = apply(admin2.samples,2,quantile,prob = 0.95))
  admin2.res.stan <- rbind(admin2.res.stan,admin2.res.tmp)
  
  adm1_wt_mat <- matrix(0,n_admin2,n_admin1)
  for(area in 1:n_admin1){
    adm2_areas_tmp <- admin.key[admin.key$A1==area,]$A2
    adm1_wt_mat[adm2_areas_tmp,area] <- admin2_weights[admin2_weights$A2 %in% adm2_areas_tmp,]$adm1.wt/sum(admin2_weights[admin2_weights$A2 %in% adm2_areas_tmp,]$adm1.wt)
  }
  admin1.samples <- admin2.samples%*%adm1_wt_mat
  admin1.res.tmp <- data.frame(admin1 = 1:n_admin1, iter.id = k,
                               stan.mean = apply(admin1.samples,2,mean),
                               stan.median = apply(admin1.samples,2,median),
                               stan.sd = apply(admin1.samples,2,sd),
                               stan.lower = apply(admin1.samples,2,quantile,prob = 0.05),
                               stan.upper = apply(admin1.samples,2,quantile,prob = 0.95))
  admin1.res.stan <- rbind(admin1.res.stan,admin1.res.tmp)
}

saveRDS(admin1.res.stan,file=paste0('Stan_Admin1Summary_',setting_number,'.rds'))
saveRDS(admin2.res.stan,file=paste0('Stan_Admin2Summary_',setting_number,'.rds'))

