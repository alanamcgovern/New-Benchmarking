# required packages
suppressMessages({
  library(data.table) # loaded in cluster
  library(LaplacesDemon) # loaded in cluster
  library(writexl) # loaded in cluster
})

source('Alg5.4_NUTSfn.R')

# step size
eps0 <- 0.05

args <- commandArgs(trailingOnly=TRUE)

chain_number = as.numeric(args[1])
n_iter = as.numeric(args[2])

load(paste0('ZMB_data.rda'))

M_init = ZMB_data[[1]]$M_init
logit_r_hat = ZMB_data[[1]]$logit_r_hat
logit_r_hat_var = ZMB_data[[1]]$logit_r_hat_var
data_list = ZMB_data[[1]]$data_list


result <- postsamp_Alg5.4_NUTS_MCMC_onechain(eps0 = eps0, M_init = M_init,
                                             logit_r_hat = logit_r_hat,
                                             logit_r_hat_var = logit_r_hat_var,
                                             data_list = data_list,
                                             n_iter = n_iter)

result$chain <- chain_number

write_xlsx(result,path=paste0('ZMB/Results',chain_number,'.xlsx'))

