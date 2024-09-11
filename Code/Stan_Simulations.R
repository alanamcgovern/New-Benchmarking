suppressMessages({
  library(cmdstanr)
})

args <- commandArgs(trailingOnly=TRUE)

data_number = as.numeric(args[1])
setting_number <- as.numeric(args[2])

if(!file.exists(paste0('Sim',setting_number,'_stan/Results',data_number,'.rds'))){
  mod <- cmdstan_model('Model_compare.stan')
  
  load(paste0('sim_setting_',setting_number,'.rda'))
  load(paste0('sim',setting_number,'_datasets.rda'))
  n_admin1 <- length(unique(sim_setting$all_dat$A1))
  n_admin2 <- length(unique(sim_setting$all_dat$A2))
  
  list_dat <- list(lenA1=n_admin1, lenA2=n_admin2,
                   lenC=nrow(sim_data[[data_number]]$obs_dat),
                   admin1_id=sim_data[[data_number]]$obs_dat$A1,
                   admin2_id=sim_data[[data_number]]$obs_dat$A2,
                   urban_id=sim_data[[data_number]]$obs_dat$U,
                   Y=sim_data[[data_number]]$obs_dat$Z,
                   N=sim_data[[data_number]]$obs_dat$n,
                   Q_scaled_inv = sim_data[[data_number]]$Q_scaled_inv)
  
  stan.fit <- mod$sample(
    data = list_dat,
    chains = 4,
    parallel_chains = 4,
    refresh=1000
  )
  
  out <- list(hyper = as.data.frame(stan.fit$draws(format='df',variables=c('tau','phi','d'))),
              alpha = as.data.frame(stan.fit$draws(format='df',variables=c('alpha'))),
              beta = as.data.frame(stan.fit$draws(format='df',variables=c('beta'))),
              b = as.data.frame(stan.fit$draws(format='df',variables=c('b'))))
  
  saveRDS(out,file=paste0('Sim',setting_number,'_stan/Results',data_number,'.rds'))
}
