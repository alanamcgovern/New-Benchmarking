# required packages
suppressMessages({
 # library(tidyverse) # loaded in cluster
  library(data.table) # loaded in cluster
  library(LaplacesDemon) # loaded in cluster
  library(writexl) # loaded in cluster
})
source('Alg5.4_NUTSfn.R')

# step size
eps0 <- 0.05


args <- commandArgs(trailingOnly=TRUE)

data_number = as.numeric(args[1])
chain_number = as.numeric(args[2])
n_iter = as.numeric(args[3])
setting_number <- as.numeric(args[4])


if(!file.exists(paste0('Sim',setting_number,'/Results',data_number,'-',chain_number,'.xlsx')) & 
   !file.exists(paste0('Sim',setting_number,'/Results',data_number,'.rds'))){
  
  load(paste0('sim',setting_number,'_datasets.rda'))
  
  M_init = sim_data[[data_number]]$M_init
  logit_r_hat = sim_data[[data_number]]$logit_r_hat
  logit_r_hat_var = sim_data[[data_number]]$logit_r_hat_var
  data_list = list(obs_dat = sim_data[[data_number]]$obs_dat,
                   pop_strata_dat = sim_data[[data_number]]$pop_strata_dat,
                   Q_scaled_inv = sim_data[[data_number]]$Q_scaled_inv)
  
  
  result <- postsamp_Alg5.4_NUTS_MCMC_onechain(eps0 = eps0, M_init = M_init,
                                               logit_r_hat = logit_r_hat,
                                               logit_r_hat_var = logit_r_hat_var,
                                               data_list = data_list,
                                               n_iter = n_iter)

  result$chain <- chain_number
  result$dataset <- data_number

  write_xlsx(result,path=paste0('Sim',setting_number,'/Results',data_number,'-',chain_number,'.xlsx'))
}

