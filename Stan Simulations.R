library(rstan)
library(tidyverse)
library(locfit)
library(INLA)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]
setwd(paste(code.path.splitted[1: (length(code.path.splitted)-1)], collapse = "/"))
source('Data Generation.R')

#next steps: move into sampling context with SRS (don't know all births)

#Alg 2.1: MCMC to get P(r|Y) (fixed spatial effect) -------------
#simulate data w a fixed spatial effect
n_areas <- 5
n_clusters <- rep(15,n_areas)
b <- rnorm(n_areas,-2,0.25)

get.sims1 <- function(n_areas,n_clusters,b,nsims){
  b.res.stan <- list()
  b.res.inla <- list()
  r.res <- list()
  for(i in 1:nsims){
    message(paste0('Starting simulation ',i))
    dat1 <- simdat_fsp(n_areas,n_clusters,b)
    #sum over admin areas and put in list format
    list_dat1 <- list(lenA = n_areas,
                      N=(dat1 %>% group_by(A) %>% summarise(N=sum(N)))$N,
                      Y=(dat1 %>% group_by(A) %>% summarise(Y=sum(Y)))$Y)
    
    #run MCMC
    fit1.out <- stan(file = 'Alg2.1.stan', data = list_dat1, iter=2000)
    #sort results
    fit1 <- rstan::extract(fit1.out)
    #check against INLA
    inlafit1 <- INLA::inla(Y ~ as.factor(A) -1,
                           data=dat1, 
                           family='poisson', 
                           control.family = list(link = "log"), 
                           E = N)
    
    b.res.stan[[i]] <- fit1$b
    b.res.inla[[i]] <- inlafit1$summary.fixed[,1:2]
    r.res[[i]] <- fit1$r
  }
  return(list(b.res.stan,b.res.inla,r.res))
}

results1 <- get.sims1(n_areas,n_clusters,b,1)
results1.b <- results1[[1]]
results1.r <- results1[[3]]

#compare STAN and INLA
lapply(results1[[1]],colMeans)
results1[[2]]

#some traceplots
plot.ts(results1.b[[1]],type='l')
plot.ts(results1.b[[1]][1000:2000,],type='l')

#Alg 2.2: MCMC to get P(r|Y,Y+) (fixed spatial effect) -------------
#simulate data w a fixed spatial effect
n_areas <- 5
n_clusters <- rep(15,n_areas)
b <- rnorm(n_areas,-2,0.25)

get.sims2 <- function(n_areas,n_clusters,b,nsims){
  b.res.stan <- list()
  #b.res.inla <- list()
  r.res <- list()
  for(i in 1:nsims){
    message(paste0('Starting simulation ',i))
    dat2 <- simdat_fsp(n_areas,rep(15,n_areas),b)
    #sum over admin areas and put in list format
    list_dat2 <- list(lenA = n_areas,
                  N=(dat2 %>% group_by(A) %>% summarise(N=sum(N)))$N,
                  Y=(dat2 %>% group_by(A) %>% summarise(Y=sum(Y)))$Y,
                  Yplus=sum(dat2$Y))

    #run MCMC
    fit2.out <- stan(file = 'Alg2.2.stan', data = list_dat2,iter = 5000)
    fit2 <- rstan::extract(fit2.out)
    
    #is there a way to do this in INLA -- potentially, need to ask Serge? or Geir Arne?
    
    b.res.stan[[i]] <- fit2$b
    #b.res.inla[[i]] <- inlafit1$summary.fixed[,1:2]
    r.res[[i]] <- fit2$r
  }
  return(list(b.res.stan,r.res))
}

results2 <- get.sims2(n_areas,n_clusters,b,1)
lapply(results2[[1]],colMeans)


#Alg 2.3a: MCMC to get P(r,phi,tau|Y) (IID spatial effect) + Alg 2.4a: MCMC to get P(r, phi, tau|Y,Y+) (IID spatial effect) -------
## define ground truth 
{
n_areas <- 40
n_clusters <- rep(10,n_areas)
#spatial relationship
Q <- matrix(NA,n_areas,n_areas)
for(i in 1:(n_areas-1)){
  Q[i,(i+1):n_areas] <- -rbinom(n_areas-i,1,0.25)
  if(sum(Q[i,(i+1):n_areas])==0){
    Q[i,i+1] <- -1
  }
  Q[(i+1):n_areas,i] <- Q[i,(i+1):n_areas]
}
for(i in 1:n_areas){
  Q[i,i] <- -sum(Q[i,],na.rm = T)
}

#intercept (fixed effect) 
a <- rnorm(1,-6,0.25)

#draw spatial precision from hyperprior (same as INLA default)
tau <- rgamma(1,1,0.5)

# draw random effects
b <- Rfast::rmvnorm(1,rep(0,n_areas),diag(1/tau,n_areas))


}
## run simulations 
{
nsims <- 100
a.res.stan3 <- a.res.stan4 <- a.res.inla <- a.res.eb <- list()
b.res.stan3 <- b.res.stan4 <- list()
r.res.inla <- r.res.eb <- r.res.mle <- r.res.stan3 <- r.res.stan4 <- list()
tau.res.stan3 <- tau.res.stan4 <- tau.res.inla <- tau.res.eb <- tau.res.mle <- list()
for(i in 1:nsims){
  message(paste0('Starting simulation ',i))
  datiid.info <- simdat_iid(n_areas=n_areas,n_clusters=n_clusters,a=a,tau=tau,b=NULL)
  dat <- datiid.info$dat
  list_dat3 <- list(lenA=n_areas,
                    N=(dat %>% group_by(A) %>% summarise(N=sum(N)))$N,
                    Y=(dat %>% group_by(A) %>% summarise(Y=sum(Y)))$Y)
  list_dat4 <- list(lenA=n_areas,
                    N=(dat %>% group_by(A) %>% summarise(N=sum(N)))$N,
                    Y=(dat %>% group_by(A) %>% summarise(Y=sum(Y)))$Y,
                    Yplus=sum(dat$Y))
  
  #check against MLE
  r.res.mle[[i]] <- (dat %>% group_by(A) %>% summarise(mle=sum(Y)/sum(N)))$mle
  tau.res.mle[[i]] <- 1/var(as.vector(datiid.info$b))
  
  # check against empirical Bayes
  ebfit <- INLA::inla(Y ~ f(A, graph = abs(Q), model = "iid"),
                      data=dat, family='poisson', E=N,
                      control.inla = list(int.strategy = "eb"),
                      control.predictor = list(compute = FALSE, link = 1))
  # 
  # #some INLA results
  a.res.eb[[i]] <- ebfit$summary.fixed[,1]
  r.res.eb[[i]] <- unique(ebfit$summary.fitted.values[,1])
  tau.res.eb[[i]] <- ebfit$summary.hyperpar[1,1]
  
  # #check against INLA
  hyper.inla <- list(prec = list(prior = "loggamma", param = c(0.01,0.01)))
  
  inlafit3 <- INLA::inla(Y ~ f(A, graph = abs(Q), model = "iid", hyper=hyper.inla),
                         data=dat, family='poisson', E=N,
                         control.predictor = list(compute = FALSE, link = 1))
  # 
  # #some INLA results
  a.res.inla[[i]] <- inlafit3$summary.fixed[,1]
  r.res.inla[[i]] <- unique(inlafit3$summary.fitted.values[,1])
  tau.res.inla[[i]] <- inlafit3$summary.hyperpar[1,1]
  
  
  # #run MCMC Alg 2.3
  fit3.out <- stan(file = 'Alg2.3a.stan', data = list_dat3)
  fit3 <- rstan::extract(fit3.out)

  # #save STAN results
  a.res.stan3[[i]] <- fit3$a
  b.res.stan3[[i]] <- fit3$b
  r.res.stan3[[i]] <- fit3$r
  tau.res.stan3[[i]] <- fit3$tau


  
  #run MCMC Alg 2.4
  fit4.out <- stan(file = 'Alg2.4a.stan', data = list_dat4)
  fit4 <- rstan::extract(fit4.out)

  # #save STAN results
  a.res.stan4[[i]] <- fit4$a
  b.res.stan4[[i]] <- fit4$b
  r.res.stan4[[i]] <- fit4$r
  tau.res.stan4[[i]] <- fit4$tau
}
}
## organize results 
{
-------
# stan: for each parameter of interest we get a list whose elements are sublists (one for each iteration)
# inla: for each parameter of interest we get a list whose elements are the parameter estimates for each iteration
results <- (list(stan3_a = a.res.stan3, stan3_b = b.res.stan3, stan3_r = r.res.stan3,
                 stan3_tau = tau.res.stan3, 
                 
                 stan4_a = a.res.stan4, stan4_b = b.res.stan4, stan4_r = r.res.stan4,
                 stan4_tau = tau.res.stan4,
                 
                 mle_r = r.res.mle,
                 mle_tau = tau.res.mle,
                 
                 inla_a = a.res.inla,inla_r= r.res.inla,
                 inla_tau = tau.res.inla,
                 
                 eb_a = a.res.eb,eb_r= r.res.eb,
                 eb_tau = tau.res.eb))

summary_results <- list(res_a=rbind(data.frame(a=unlist(lapply(results$stan3_a,mean)),model='stan3'),
                                 data.frame(a=unlist(lapply(results$stan4_a,mean)),model='stan4'),
                                 data.frame(a=unlist(results$inla_a),model='inla'),
                                 data.frame(a=unlist(results$eb_a),model='eb')),
                     res_b=rbind(data.frame(matrix(unlist(lapply(results$stan3_b,colmeans)),ncol=n_areas,byrow = T),model='stan3'),
                                 data.frame(matrix(unlist(lapply(results$stan4_b,colmeans)),ncol=n_areas,byrow = T),model='stan4')),
                     res_r=rbind(pivot_longer(data.frame(matrix(unlist(lapply(results$stan3_r,colmeans)),ncol=n_areas,byrow = T),model='stan3'),cols=X1:X15,names_to = 'area',values_to = 'r'),
                                 pivot_longer(data.frame(matrix(unlist(lapply(results$stan4_r,colmeans)),ncol=n_areas,byrow = T),model='stan4'),cols=X1:X15,names_to = 'area',values_to = 'r'),
                                 pivot_longer(data.frame(matrix(unlist(results$inla_r),ncol=n_areas,byrow = T),model='inla'),cols=X1:X15,names_to = 'area',values_to = 'r'),
                                 pivot_longer(data.frame(matrix(unlist(results$eb_r),ncol=n_areas,byrow = T),model='eb'),cols=X1:X15,names_to = 'area',values_to = 'r'),
                                 pivot_longer(data.frame(matrix(unlist(results$mle_r),ncol=n_areas,byrow = T),model='mle'),cols=X1:X15,names_to = 'area',values_to = 'r')),
                     res_tau=rbind(data.frame(tau=unlist(lapply(results$stan3_tau,mean)),model='stan3'),
                                   data.frame(tau=unlist(lapply(results$stan4_tau,mean)),model='stan4'),
                                   data.frame(tau=unlist(results$inla_tau),model='inla'),
                                   data.frame(tau=unlist(results$eb_tau),model='eb'),
                                   data.frame(tau=unlist(results$mle_tau),model='mle')))

summary_results$res_r$area <- as.numeric(str_remove(summary_results$res_r$area,'X'))

}
## save results 
{
-------
 alg_res <- list(truth=list(a=a,tau=tau),results=summary_results)
 save(alg_res,file='alg_res230615_random b N=5000x15.rda')
 
 # alg_res <- list(truth=list(a=a,b=b,tau=tau),results=summary_results)
 # save(alg_res,file='alg_res230615_fixed b N=5000x15.rda')

 }
## visualize results
{

# pdf(file='Alg Results 230615 fixed b N=5000x15.pdf')
# {
# 
#   #compare r estimates (risk) for each area
#   n_pages <- ceiling(n_areas/9)
#   truth_dat <- data.frame(cbind(area=1:n_areas,r=exp(a+as.vector(b))))
#   for(page in 1:n_pages){
#     area_ids <- (9*(page-1)+1):(9*page)
#     area_ids <- area_ids[area_ids<=n_areas]
#     results_sub <- summary_results$res_r %>% filter(area %in% area_ids)
#     truth_dat <- data.frame(cbind(area=area_ids,r=exp(a+as.vector(b)[area_ids])))
#     g <- ggplot(results_sub) + geom_boxplot(aes(x=model,y=r)) + geom_hline(data=truth_dat, aes(yintercept=r),col='red') +
#       ylab('') + ylim(c(max(exp(a+b[i])-0.001,0)),exp(a+b[i])+0.001) + facet_wrap(~area)
#     print(g)
#   }
# 
# 
#   #compare tau estimates
#   ggplot(summary_results$res_tau) + geom_boxplot(aes(x=model,y=log(tau))) + geom_hline(yintercept=(log(tau)),col='red') +
#     ggtitle('Spatial precision, tau')
# 
# }
# dev.off()

pdf(file='Alg Results 230615 random b N=5000x15.pdf')
{

  #compare tau estimates
  ggplot(summary_results$res_tau) + geom_boxplot(aes(x=model,y=log(tau))) + geom_hline(yintercept=(log(tau)),col='red') +
    ggtitle('Spatial precision, tau')

}
dev.off()


}
#Alg 2.3: MCMC to get P(r,phi,tau|Y) (BYM2 spatial effect) + Alg 2.4: MCMC to get P(r, phi, tau|Y,Y+) (BYM2 spatial effect) -------------

## define ground truth 
{
n_areas <- 75
n_clusters <- rep(10,n_areas)
#spatial relationship
Q <- matrix(NA,n_areas,n_areas)
for(i in 1:(n_areas-1)){
  Q[i,(i+1):n_areas] <- -rbinom(n_areas-i,1,0.25)
  if(sum(Q[i,(i+1):n_areas])==0){
    Q[i,i+1] <- -1
  }
  Q[(i+1):n_areas,i] <- Q[i,(i+1):n_areas]
}
for(i in 1:n_areas){
  Q[i,i] <- -sum(Q[i,],na.rm = T)
}

#intercept (fixed effect) 
a <- rnorm(1,-6,0.25)

#draw spatial precision from hyperprior
log_tau_b <- rnorm(1,3,0.5)
tau_b <- exp(log_tau_b)

#draw spatial variation from beta hyperprior
phi <- rbeta(1,2,3)

}
## arguments for BYM2 
{
## reparametrized matrix for BYM2 parameters
Q_scaled <- inla.scale.model(Q, constr=list(A=t(eigen(Q)$vectors[,eigen(Q)$values<1e-10]), e=rep(0,sum(eigen(Q)$values<1e-10))))
Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))
Var_b <- (1/tau_b)*(diag((1-phi),n_areas) + phi*Q_scaled_inv)
b <- Rfast::rmvnorm(1,rep(0,n_areas),Var_b)

## find lambda to properly calibrate PC prior on phi
alpha <- 2/3
U <- 0.5

kld_phi <- function(phi, Q_inv){
  out <- phi*tr(Q_inv) - phi*nrow(Q_inv) - as.numeric(determinant(diag(1-phi,nrow(Q_inv)) + phi*Q_inv,logarithm = T)$modulus)
  return(out/2)
}

gamma <- eigen(Q_scaled)$values
gamma_til <- 1/gamma
gamma_til[gamma_til>1e+10] <- 0

pc_pdf_phi<- function(phi,lambda,Q_inv,gamma_til){
  out <- lambda/sqrt(8*kld_phi(phi,Q_inv))*exp(-lambda*sqrt(2*kld_phi(phi,Q_inv)))*
    (tr(Q_inv)-nrow(Q_inv)-sum((gamma_til-1)/(1+phi*(gamma_til-1))))
  return(out)
}

#use prior from INLA to calibrate to proper lambda for given alpha and u
pc_pdf_phi_inla <- INLA:::inla.pc.bym.phi(Q=Q_scaled,rankdef = 1,alpha = alpha,u=U)
prior_diff <- function(lambda,Q_scaled_inv,gamma_til){
  phi_seq <- seq(0.01,0.99,0.01)
  sum((sapply(phi_seq,function(x){exp(pc_pdf_phi_inla(x))}) - sapply(phi_seq,function(x){pc_pdf_phi(x,lambda,Q_scaled_inv,gamma_til)}))^2)
}

lambda_opt <- optim(1,prior_diff,Q_scaled_inv=Q_scaled_inv,gamma_til=gamma_til,method = 'Brent',lower = 0,upper=10)$par


}
## run simulations 
{
nsims <- 100
a.res.stan3 <- a.res.stan4 <- a.res.inla <-a.res.inla.pc <- a.res.eb <- list()
b.res.stan3 <- b.res.stan4 <- list()
r.res.inla <- r.res.inla.pc <- r.res.mle <- r.res.eb <- r.res.stan3 <- r.res.stan4 <- list()
tau.res.stan3 <- phi.res.stan3 <- tau.res.stan4 <- phi.res.stan4 <- tau.res.inla <- phi.res.inla <- tau.res.inla.pc <- phi.res.inla.pc <- tau.res.mle <- tau.res.eb <- phi.res.eb <- list()
  for(i in 1:nsims){
    message(paste0('Starting simulation ',i))
    datbym2.info <- simdat_bym2sp(n_areas=n_areas,n_clusters=n_clusters,n_births=500,
                                  a=a,Var_b=Var_b,b=b)
    dat <- datbym2.info$dat
    list_dat3 <- list(lenA=n_areas,
                      N=(dat %>% group_by(A) %>% summarise(N=sum(N)))$N,
                      Y=(dat %>% group_by(A) %>% summarise(Y=sum(Y)))$Y,
                      Q_scaled_inv = Q_scaled_inv,
                      eigen = gamma_til,
                      lambda = lambda_opt)
                      
    list_dat4 <- list(lenA=n_areas,
                      N=(dat %>% group_by(A) %>% summarise(N=sum(N)))$N,
                      Y=(dat %>% group_by(A) %>% summarise(Y=sum(Y)))$Y,
                      Yplus=sum(dat$Y),
                      Q_scaled_inv = Q_scaled_inv,
                      eigen = gamma_til,
                      lambda = lambda_opt)
                      
    #run MCMC Alg 2.3
    fit3.out <- stan(file = 'Alg2.3.stan', data = list_dat3)
    fit3 <- rstan::extract(fit3.out)

    # #save STAN results
    a.res.stan3[[i]] <- fit3$a
    b.res.stan3[[i]] <- fit3$b
    r.res.stan3[[i]] <- fit3$r
    tau.res.stan3[[i]] <- fit3$tau
    phi.res.stan3[[i]] <- fit3$phi

    # #run MCMC Alg 2.4
    fit4.out <- stan(file = 'Alg2.4.stan', data = list_dat4)
    fit4 <- rstan::extract(fit4.out)

    # #save STAN results
    a.res.stan4[[i]] <- fit4$a
    b.res.stan4[[i]] <- fit4$b
    r.res.stan4[[i]] <- fit4$r
    tau.res.stan4[[i]] <- fit4$tau
    phi.res.stan4[[i]] <- fit4$phi

    # check against empirical Bayes
      # ebfit <- INLA::inla(Y ~ f(A, graph = abs(Q), model = "bym2",
      #                           scale.model = TRUE,
      #                           adjust.for.con.comp = TRUE),
      #                     data=dat, family='poisson', E=N,
      #                     control.inla = list(int.strategy = "eb"),
      #                     control.predictor = list(compute = FALSE, link = 1))
      # 
      # #some INLA results
      # a.res.eb[[i]] <- ebfit$summary.fixed[,1]
      # r.res.eb[[i]] <- unique(ebfit$summary.fitted.values[,1])
      # tau.res.eb[[i]] <- ebfit$summary.hyperpar[1,1]
      # phi.res.eb[[i]] <- ebfit$summary.hyperpar[2,1]

    #check against INLA
    # hyper.inla <- list(prec = list(prior = "pc.prec", param = c(1,0.01), initial=5),
    #                    phi = list(prior= "logitbeta", param=c(2,3)))
    # 
    # inlafit <- INLA::inla(Y ~ f(A, graph = abs(Q),
    #                              model = "bym2",
    #                              hyper = hyper.inla,
    #                              scale.model = TRUE,
    #                              adjust.for.con.comp = TRUE),
    #                        data=dat, family='poisson', E=N,
    #                        control.predictor = list(compute = FALSE, link = 1))
    # 
    # #some INLA results
    # a.res.inla[[i]] <- inlafit$summary.fixed[,1]
    # r.res.inla[[i]] <- unique(inlafit$summary.fitted.values[,1])
    # tau.res.inla[[i]] <- inlafit$summary.hyperpar[1,1]
    # phi.res.inla[[i]] <- inlafit$summary.hyperpar[2,1]
    
    # check against INLA with all PC priors
    hyper.inla.pc <- list(prec = list(prior = "pc.prec", param = c(1,0.01), initial=5),
                       phi = list(prior= "pc", param=c(0.5,2/3)))

    inlafit.pc <- INLA::inla(Y ~ f(A, graph = abs(Q),
                                model = "bym2",
                                hyper = hyper.inla.pc,
                                scale.model = TRUE,
                                adjust.for.con.comp = TRUE),
                          data=dat, family='poisson', E=N,
                          control.predictor = list(compute = FALSE, link = 1))
    
    # #some INLA results
    a.res.inla.pc[[i]] <- inlafit.pc$summary.fixed[,1]
    r.res.inla.pc[[i]] <- unique(inlafit.pc$summary.fitted.values[,1])
    tau.res.inla.pc[[i]] <- inlafit.pc$summary.hyperpar[1,1]
    phi.res.inla.pc[[i]] <- inlafit.pc$summary.hyperpar[2,1]
    
    #check against MLE
    r.res.mle[[i]] <- (dat %>% group_by(A) %>% summarise(mle=sum(Y)/sum(N)))$mle
    # can only solve for known phi
    tau.res.mle[[i]] <- n_areas/(datbym2.info$b%*%solve(diag(1-phi,n_areas) + phi*Q_scaled_inv)%*%t(datbym2.info$b))
    
  }

}
## organize results 
{
# stan: for each parameter of interest we get a list whose elements are sublists (one for each iteration)
# inla: for each parameter of interest we get a list whose elements are the parameter estimates for each iteration
results <- (list(
                  stan3_a = a.res.stan3, stan3_b = b.res.stan3, stan3_r = r.res.stan3,
                  stan3_tau = tau.res.stan3, stan3_phi = phi.res.stan3,
                 
                 stan4_a = a.res.stan4, stan4_b = b.res.stan4, stan4_r = r.res.stan4,
                 stan4_tau = tau.res.stan4, stan4_phi = phi.res.stan4,
                 
                 mle_tau = tau.res.mle,
                 
                 # eb_a = a.res.eb,eb_r= r.res.eb,
                 # eb_tau = tau.res.eb, eb_phi = phi.res.eb,
                 
                 # inla_a = a.res.inla,inla_r= r.res.inla,
                 # inla_tau = tau.res.inla, inla_phi = phi.res.inla,
                  
                 inla_pc_a = a.res.inla.pc,inla_pc_r= r.res.inla.pc,
                 inla_pc_tau = tau.res.inla.pc, inla_pc_phi = phi.res.inla.pc))
              
summary_results <- list(res_a=rbind(data.frame(a=unlist(lapply(results$stan3_a,mean)),model='stan3'),
                                 data.frame(a=unlist(lapply(results$stan4_a,mean)),model='stan4'),
                                 #data.frame(a= unlist(results$eb_a),model='eb'),
                                 #data.frame(a= unlist(results$inla_a),model='inla'),
                                 data.frame(a= unlist(results$inla_pc_a),model='inla.pc')),
                     res_b=rbind(data.frame(matrix(unlist(lapply(results$stan3_b,colmeans)),ncol=n_areas,byrow = T),model='stan3'),
                                  data.frame(matrix(unlist(lapply(results$stan4_b,colmeans)),ncol=n_areas,byrow = T),model='stan4')),
                     res_r=rbind(data.frame(matrix(unlist(lapply(results$stan3_r,colmeans)),ncol=n_areas,byrow = T),model='stan3'),
                                 data.frame(matrix(unlist(lapply(results$stan4_r,colmeans)),ncol=n_areas,byrow = T),model='stan4'),
                                 #data.frame(matrix(unlist(results$eb_r),ncol=n_areas,byrow = T),model='eb'),
                                 #data.frame(matrix(unlist(results$inla_r),ncol=n_areas,byrow = T),model='inla'),
                                 data.frame(matrix(unlist(results$inla_pc_r),ncol=n_areas,byrow = T),model='inla.pc')),
                               #  data.frame(matrix(unlist(results$mle_r),ncol=n_areas,byrow = T),model='mle')),
                     res_tau=rbind(data.frame(tau=unlist(lapply(results$stan3_tau,mean)),model='stan3'),
                                   data.frame(tau=unlist(lapply(results$stan4_tau,mean)),model='stan4'),
                                 #  data.frame(tau=unlist(results$eb_tau),model='eb'),
                                 # data.frame(tau=unlist(results$inla_tau),model='inla'),
                                   data.frame(tau=unlist(results$inla_pc_tau),model='inla.pc'),
                                   data.frame(tau=unlist(results$mle_tau),model='mle')),
                     res_phi=rbind(data.frame(phi=unlist(lapply(results$stan3_phi,mean)),model='stan3'),
                                   data.frame(phi=unlist(lapply(results$stan4_phi,mean)),model='stan4'),
                                #   data.frame(phi=unlist(results$eb_phi),model='eb'),
                                #   data.frame(phi=unlist(results$inla_phi),model='inla'),
                                   data.frame(phi=unlist(results$inla_pc_phi),model='inla.pc')))
}
## save results 
{
   prior_info <- list(#inla.tau = 'PC(1,0.01)', inla.phi='Beta(2,3)', 
                      stan.tau = 'PC(1,0.01)', stan.phi = 'PC(0.5,2/3)')
   alg_res <- list(truth=list(a=a,b=b,tau_b=tau_b,phi=phi),results=summary_results,prior_info=prior_info)
   save(alg_res,file='bym2_res230711 fixed b N=5000x75.rda')

   }
## visualize results 
{
# 
 pdf(file='BYM2 Results 230711 fixed b N=5000x75.pdf')
{
  
  #compare r estimates (risk) for each area
    n_pages <- ceiling(n_areas/9)
    truth_dat <- data.frame(cbind(area=1:n_areas,r=exp(a+as.vector(b))))
    for(page in 1:n_pages){
      area_ids <- (9*(page-1)+1):(9*page)
      area_ids <- area_ids[area_ids<=n_areas]
      results_sub <- gather(summary_results$res_r[,c(area_ids,n_areas+1)],key='area',value='r',-model)
      results_sub$area <- as.numeric(str_remove(results_sub$area,'X'))
      g <- ggplot(results_sub) + geom_boxplot(aes(x=model,y=r)) + geom_hline(data=truth_dat[area_ids,], aes(yintercept=r),col='red') +
        geom_hline(yintercept = exp(a),lty=2,color='red') +
        ylab('') + facet_wrap(~area)
      print(g)
    }
  
#compare tau estimates
g <- ggplot(summary_results$res_tau) + geom_boxplot(aes(x=model,y=log(tau))) + geom_hline(yintercept=log(tau_b),col='red') +
    ggtitle('Spatial precision, tau')
  print(g)
  
#compare phi estimates
g <- ggplot(summary_results$res_phi) + geom_boxplot(aes(x=model,y=logitlink(phi))) + geom_hline(yintercept=logitlink(phi),col='red') +
    ggtitle('Spatial variance, phi')
print(g)

}
dev.off()

#Alg 3.1: P(Y,r,phi,tau|Z) and P(Y,r,phi,tau|Z,Y+)
}
#Alg 3.1 -----
## define ground truth 
{
  ## this section is calibrated to be relatively similar to Sierra Leone 2019 survey, in terms of sample sizes
  n_admin2 <- 20
  n_clusters_urban <- rnegbin(n_admin2,200,5) #number of urban clusters in each admin2 area
  n_clusters_rural <- rnegbin(n_admin2,300,5) #number of rural clusters in each admin2 area
  n_births_urban <- 75 #average number of births per urban cluster (0.3*number of households) 
  n_births_rural <-  100 #average number of births per rural cluster (0.45*number of households) 
  n_clusters_urban_samp <-  round(0.075*n_clusters_urban) #number of urban clusters sampled in each admin2 area
  n_clusters_rural_samp <-  round(0.05*n_clusters_rural) #number of rural clusters sampled in each admin2 area
  n_births_urban_samp <- 20 # average number of births sampled from each urban cluster 
  n_births_rural_samp <- 30 # average number of births sampled from each rural cluster 
  
  #spatial relationship
  Q <- matrix(NA,n_admin2,n_admin2)
  for(i in 1:(n_admin2-1)){
    Q[i,(i+1):n_admin2] <- -rbinom(n_admin2-i,1,0.25)
    if(sum(Q[i,(i+1):n_admin2])==0){
      Q[i,i+1] <- -1
    }
    Q[(i+1):n_admin2,i] <- Q[i,(i+1):n_admin2]
  }
  for(i in 1:n_admin2){
    Q[i,i] <- -sum(Q[i,],na.rm = T)
  }
  
  # intercept (fixed effect) 
  alphaU <- rnorm(1,-3.5,0.25)
  alphaR <- rnorm(1,-3.2,0.25)
  
  # draw hyperparameters
  tau_b <- rgamma(1,10,0.25)
  
  #draw spatial variation from beta hyperprior
  phi <- rbeta(1,2,3)
  
  # draw overdispersion parameter
  d <- rlnorm(1,0,0.25)
 
}

## arguments for BYM2 
{
## reparametrized matrix for BYM2 parameters
Q_scaled <- inla.scale.model(Q, constr=list(A=t(eigen(Q)$vectors[,eigen(Q)$values<1e-10]), e=rep(0,sum(eigen(Q)$values<1e-10))))
Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))
Var_b <- (1/tau_b)*(diag((1-phi),n_admin2) + phi*Q_scaled_inv)
b <- Rfast::rmvnorm(1,rep(0,n_admin2),Var_b)

## find lambda to properly calibrate PC prior on phi
alpha <- 2/3
U <- 0.5

kld_phi <- function(phi, Q_inv){
  out <- phi*tr(Q_inv) - phi*nrow(Q_inv) - as.numeric(determinant(diag(1-phi,nrow(Q_inv)) + phi*Q_inv,logarithm = T)$modulus)
  return(out/2)
}

gamma <- eigen(Q_scaled)$values
gamma_til <- 1/gamma
gamma_til[gamma_til>1e+10] <- 0

pc_pdf_phi<- function(phi,lambda,Q_inv,gamma_til){
  out <- lambda/sqrt(8*kld_phi(phi,Q_inv))*exp(-lambda*sqrt(2*kld_phi(phi,Q_inv)))*
    (tr(Q_inv)-nrow(Q_inv)-sum((gamma_til-1)/(1+phi*(gamma_til-1))))
  return(out)
}

#use prior from INLA to calibrate to proper lambda for given alpha and u
pc_pdf_phi_inla <- INLA:::inla.pc.bym.phi(Q=Q_scaled,rankdef = 1,alpha = alpha,u=U)
prior_diff <- function(lambda,Q_scaled_inv,gamma_til){
  phi_seq <- seq(0.01,0.99,0.01)
  sum((sapply(phi_seq,function(x){exp(pc_pdf_phi_inla(x))}) - sapply(phi_seq,function(x){pc_pdf_phi(x,lambda,Q_scaled_inv,gamma_til)}))^2)
}

lambda_opt <- optim(1,prior_diff,Q_scaled_inv=Q_scaled_inv,gamma_til=gamma_til,method = 'Brent',lower = 0,upper=10)$par

}

## generate data 
{
  admin_key <- data.frame(A2=1:n_admin2)
  all_dat <- data.frame(
    # (true) number of births in each cluster
    N = round(c(rnorm(sum(n_clusters_urban),n_births_urban,10),rnorm(sum(n_clusters_rural),n_births_rural,12))),
    # urban or rural strata (U=1 if urban, otherwise 0)
    U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
    # admin area of each cluster
    A2 = c(unlist(sapply(1:n_admin2,function(i){rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_admin2,function(i){rep(i,n_clusters_rural[i])}))))
  all_dat <- merge(admin_key,all_dat)
  all_dat$cluster <- 1:nrow(all_dat)
  all_dat$Y <- sapply(1:nrow(all_dat),function(i){
    rnbinom(1, size = 1/d*all_dat$N[i]*exp(alphaU*all_dat$U[i] + alphaR*(1-all_dat$U[i]) + b[all_dat$A2[i]] + rnorm(1,0,0.05)), prob=1/(1+d))})
  
  ### sample clusters
  obs_dat <- NULL
  for(area in 1:n_admin2){
    # randomly select clusters from urban strata based on size of cluster
    obs_dat <- rbind(obs_dat,all_dat[all_dat$A2==area & all_dat$U==1,][sample(1:n_clusters_urban[area],n_clusters_urban_samp[area],prob = all_dat[all_dat$U==1 & all_dat$A2==area,]$N),],
                     # randomly select cluster from rural strata  based on size of cluster
                     all_dat[all_dat$A2==area & all_dat$U==0,][sample(1:n_clusters_rural[area],n_clusters_rural_samp[area],prob = all_dat[all_dat$U==0& all_dat$A2==area,]$N),])
  }
  
  # sample births
  obs_dat$n <- obs_dat$Z <- NA
  obs_dat[obs_dat$U==1,]$n <- round(min(rnorm(sum(obs_dat$U),n_births_urban_samp,n_births_urban_samp/10),obs_dat[obs_dat$U==1,]$N))
  obs_dat[obs_dat$U==0,]$n <- round(min(rnorm(sum(1-obs_dat$U),n_births_rural_samp,n_births_rural_samp/10),obs_dat[obs_dat$U==0,]$N))
  for(i in 1:nrow(obs_dat)){
    obs_dat$Z[i] <- sum(sample(c(rep(1,obs_dat$Y[i]),rep(0,(obs_dat$N[i]-obs_dat$Y[i]))), obs_dat$n[i]))
  }
  
  # assign sampling weights
  #births in strata/(# cluster sampled * #births sampled in cluster)
  pop_dat <- all_dat %>% group_by(A2,U) %>% summarise(N=sum(N))
  pop_dat$num_clust_samp <- pop_dat$num_births_samp <- NA
  pop_dat$num_clust_samp[pop_dat$U==0] <- n_clusters_rural_samp
  pop_dat$num_clust_samp[pop_dat$U==1] <- n_clusters_urban_samp
  pop_dat$num_births_samp[pop_dat$U==0] <- n_births_rural_samp
  pop_dat$num_births_samp[pop_dat$U==1] <- n_births_urban_samp
  pop_dat$wt <- pop_dat$N/(pop_dat$num_births_samp*pop_dat$num_clust_samp)
  pop_dat$wt <- pop_dat$wt/100
  pop_dat <- pop_dat[order(pop_dat$U,pop_dat$A2),]
  obs_dat <- merge(obs_dat,pop_dat[c('A2','U','wt')])
  obs_dat <- obs_dat[order(obs_dat$A2),]
  
  data_list <- list(obs_dat = obs_dat, #observed data
                    Yplus = sum(all_dat$Y), #ALL deaths in each admin1 area
                    pop_strata_dat = pop_dat) # number of births per strata (admin2 x urban)
  
}
## run model
{

# run MCMC Alg 3.1
list_dat <- list(lenA=n_admin2,
                 lenC=nrow(obs_dat),
                 admin2_id = obs_dat$A2,
                 urban_id = 2-obs_dat$U,
                 Ntot = pop_dat$N, #MAKE SURE THIS IS IN CORRECT ORDER (ALL URBAN THEN ALL RURAL)
                 N=obs_dat$N,
                 Y=obs_dat$Y,
                 Yplus=sum(all_dat$Y), 
                 Q_scaled_inv=Q_scaled_inv,
                 eigen = gamma_til,
                 lambda = lambda_opt,
                 alpha_prior_mu = c(-3.5,-3.2))
  # init_fun <- function(){
  #   list(mu=4.5, sigma2=0.05)}
  # 
  fit.out <- stan(file = 'Alg4.1.stan', data = list_dat)
  fit <- rstan::extract(fit.out)
  
  ## save results
  alpha.res.stan <- fit$alpha
  b.res.stan <- fit$b
  r.res.stan <- fit$r
  logd.res.stan <- fit$log_d
  tau.res.stan <- fit$tau
  phi.res.stan <- fit$phi
  
}


# APPENDIX -----

## to optimize likelihood for BYM2 parameters (MLE for phi is undefined) ----
#define BYM2 likelihood function to estimate MLE of hyperparameters
loglik_BYM2 <- function(params,b,Q_scaled_inv,m){
  log_tau <- params[1]
  logit_phi <- params[2]
  phi <- expit(logit_phi)
  logdet <- as.numeric(determinant(diag(1-phi,m) + phi*Q_scaled_inv,logarithm = T)$modulus)
  out <- m/2*log_tau - 1/2*logdet - (exp(log_tau)/2)*(b%*%solve(diag(1-phi,m) + phi*Q_scaled_inv)%*%t(b))
  return(-out)
}

#define BYM2 likelihood function to estimate MLE of hyperparameters -- just optimized as function of phi because mle of tau is written and function of phi
loglik_BYM2_onevar <- function(params,b,Q_scaled_inv,m){
  logit_phi <- params
  phi <- expit(logit_phi)
  tau <- m/(b%*%solve(diag(1-phi,m) + phi*Q_scaled_inv)%*%t(b))
  logdet <- as.numeric(determinant(diag(1-phi,m) + phi*Q_scaled_inv,logarithm = T)$modulus)
  out <- m/2*log(tau) - 1/2*logdet - (tau/2)*(b%*%solve(diag(1-phi,m) + phi*Q_scaled_inv)%*%t(b))
  return(-out)
}

# plot(seq(0,1,0.01),sapply(seq(0,1,0.01),function(x){-loglik_BYM2_onevar(logitlink(x),b,Q_scaled_inv,n_areas)}),type='l')








## how to implement PC prior for phi in stan ------

## helper functions
kld_phi <- function(phi, Q_inv){
  out <- phi*tr(Q_inv) - phi*nrow(Q_inv) - as.numeric(determinant(diag(1-phi,nrow(Q_inv)) + phi*Q_inv,logarithm = T)$modulus)
  return(out/2)
}
tr_W <- function(phi, Q_inv){
  Id_mat <- diag(1,nrow(Q_inv))
  W <- (Id_mat - INLA:::inla.ginv((1-phi)*Id_mat + phi*Q_inv)) %*% (Q_inv - Id_mat)
  return(tr(W))  
}

## prior for phi
pc_pdf_phi <- function(phi, lambda, Q_inv){
  out <- lambda/sqrt(8*kld_phi(phi,Q_inv))*exp(-lambda*sqrt(2*kld_phi(phi,Q_inv)))*abs(tr_W(phi,Q_inv))
  return(out)
}

#use prior from INLA to calibrate to proper lambda for given alpha and u
pc_pdf_phi_inla <- INLA:::inla.pc.bym.phi(Q=Q_scaled,rankdef = 1,alpha = 2/3,u=1/2)
prior_diff <- function(lambda,Q_scaled_inv){
  phi_seq <- seq(0.01,0.99,0.01)
  sum((sapply(phi_seq,function(x){exp(pc_pdf_phi_inla(x))}) - sapply(phi_seq,function(x){pc_pdf_phi(x,lambda,Q_scaled_inv)}))^2)
}

lambda_opt <- optim(1,prior_diff,Q_scaled_inv=Q_scaled_inv,method = 'Brent',lower = 0,upper=10)$par


### serge's way (gives same result) ----
gamma <- eigen(Q_scaled)$values
gamma_til <- 1/gamma
gamma_til[gamma_til>1e+10] <- 0

pc_pdf_phi_serge <- function(phi,lambda,Q_inv,gamma_til){
  out <- lambda/sqrt(8*kld_phi(phi,Q_inv))*exp(-lambda*sqrt(2*kld_phi(phi,Q_inv)))*
    (tr(Q_inv)-nrow(Q_inv)-sum((gamma_til-1)/(1+phi*(gamma_til-1))))
  return(out)
}


