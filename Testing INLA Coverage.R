
# a very simple model
log_r <- -3.5
r <- exp(log_r)
alpha <- 5
  
#sample size (number of clusters)
n <- 500
# average number of births per cluster
N <- 10

exdat <- data.frame(id=1:n)
exdat$N <- rpois(n,N)

nsim <- 100
inla.out.fv.cov <- inla.out.lp.cov <- inla.out.od.cov <- samp.out.fv.cov <- samp.out.lp.cov <- samp.out.od.cov <- rep(NA,nsim)
inla.out.fv.width <- inla.out.lp.width <- inla.out.od.width <- samp.out.fv.width <- samp.out.lp.width <- samp.out.od.width <- rep(NA,nsim)
for(i in 1:nsim){
  exdat$Y <- sapply(1:n, function(i){rnbinom(1,size = alpha, mu = exdat$N[i]*exp(log_r + rnorm(1,0,0.01)))})
  
  inla.mod1 <- inla(Y ~ 1,
                    data=exdat, family='nbinomial', E=N,
                    control.predictor = list(compute = T, link = 1),
                    control.family = list(control.link = list(model = 'log'),
                                          hyper = list(theta = list(prior="pc.mgamma", param=0.25))),
                    control.compute = list(config=T),
                    )
  
  # coverage of r from INLA (summary.fitted.values)
  inla.out.fv.cov[i] <- I(r>inla.mod1$summary.fitted.values$`0.025quant`[1] & r<inla.mod1$summary.fitted.values$`0.975quant`[1])
  inla.out.fv.width[i] <- inla.mod1$summary.fitted.values$`0.975quant`[1] - inla.mod1$summary.fitted.values$`0.025quant`[1]
  # coverage of linear predictor from INLA (summary.linear.predictor)
  inla.out.lp.cov[i] <- I(log_r>inla.mod1$summary.linear.predictor$`0.025quant`[1] & log_r<inla.mod1$summary.linear.predictor$`0.975quant`[1])
  inla.out.lp.width[i] <- inla.mod1$summary.linear.predictor$`0.975quant`[1] - inla.mod1$summary.linear.predictor$`0.025quant`[1]
  # coverage of overdispersion param from INLA (summary.hyperpar)
  inla.out.od.cov[i] <- I(alpha>inla.mod1$summary.hyperpar$`0.025quant` & alpha<inla.mod1$summary.hyperpar$`0.975quant`)
  inla.out.od.width[i] <- inla.mod1$summary.hyperpar$`0.975quant` - inla.mod1$summary.hyperpar$`0.025quant`
  
  # sample from posterior myself
  cs <- inla.mod1$misc$configs$contents$tag
  cs <- cs[cs != "Predictor"]
  select <- list()
  for (j in 1:length(cs)) {
    select[[j]] <- 0
    names(select)[j] <- cs[j]
  }
  
  sampFull <- INLA::inla.posterior.sample(n = 100, result = inla.mod1, intern = TRUE, selection = select)
  postSamples <- unlist(lapply(1:length(sampFull), function(x){sampFull[[x]]$latent}))
  # coverage of r from posterior draws
  samp.out.fv.cov[i] <-  I(r>quantile(exp(postSamples),0.025) & r<quantile(exp(postSamples),0.975))
  samp.out.fv.width[i] <- quantile(exp(postSamples),0.975) - quantile(exp(postSamples),0.025)
  # coverage of linear predictor from posterior draws
  samp.out.lp.cov[i] <- I(log_r>quantile(postSamples,0.025) & log_r<quantile(postSamples,0.975))
  samp.out.lp.width[i] <- quantile(postSamples,0.975) - quantile(postSamples,0.025)
  # coverage of overdispersion param from posterior draws
  hyperparSamples <- unlist(lapply(1:length(sampFull), function(x){sampFull[[x]]$hyperpar}))
  samp.out.od.cov[i] <- I(alpha>quantile(exp(hyperparSamples),0.025) & alpha<quantile(exp(hyperparSamples),0.975))
  samp.out.od.width[i] <- quantile(exp(hyperparSamples),0.975) - quantile(exp(hyperparSamples),0.025)
}

mean(inla.out.fv.cov)
mean(inla.out.lp.cov)
mean(samp.out.fv.cov)
mean(samp.out.lp.cov)

mean(inla.out.od.cov)
mean(samp.out.od.cov)
# 

