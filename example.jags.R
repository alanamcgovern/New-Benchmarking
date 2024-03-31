library(rjags)
library(R2jags)


dat <- list("n" = alg_dat$n, "Z" = alg_dat$Z, "U" = alg_dat$U + 1, # 1= rural, 2 = urban
            #"A1" = alg_dat$A1, 
            "A2" = alg_dat$A2,
            "prec_alpha" = diag(1/9,2),
            "prec_b" = diag(1/40,max(alg_dat$A2)),
            "n_cluster" = nrow(alg_dat))

jags.m <- jags.model( file = "JAGSmodel.txt", data=dat, n.chains=1, n.adapt=500 )
samps <- coda.samples(jags.m, variable.names = c('alpha') ,n.iter=1000 )
