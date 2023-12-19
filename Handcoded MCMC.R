library(locfit)
library(LaplacesDemon)
library(Rfast)
library(MASS)
library(INLA)
library(SUMMER)
library(VGAM)
library(tidyverse)
library(graphicalExtremes)
library(igraph)
library(parallel)

code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]
setwd(paste(code.path.splitted[1: (length(code.path.splitted)-1)], collapse = "/"))
source('Data Generation.R')

## helper functions -----
W <- function(phi,Q_inv){
  return(diag(1-phi,nrow(Q_inv)) + phi*Q_inv)
} # returns (1-phi)I + phiQ^(-1) (matrix component of BYM2 variance)

prop_mean_bym2_eta <- function(eta_star,tau,phi,alpha,data,Q_inv,Y=NULL){
  if(is.null(Y)){
    Y <- data$Y}
  out <- (Y + alpha*tau*colSums(solve(W(phi,Q_inv))) - data$N*exp(eta_star)*(1-eta_star))
  return(out%*%prop_sigma_bym2_eta(eta_star,tau,phi,data,Q_inv))
} #mean of prop dist for eta in BYM2 model

prop_sigma_bym2_eta <- function(eta_star,tau,phi,data,Q_inv){
  out <- tau*solve(W(phi,Q_inv)) + diag(data$N*exp(eta_star))
  return(solve(out))} #var of prop dist for eta in BYM2 model

# return n_areas-length vector, Niri/sum(Njrj)
theta <- function(eta,data){
  data$N*exp(eta)/sum(data$N*exp(eta))
}

kld_phi <- function(phi, Q_inv){
  out <- phi*tr(Q_inv) - phi*nrow(Q_inv) - as.numeric(determinant(diag(1-phi,nrow(Q_inv)) + phi*Q_inv,logarithm = T)$modulus)
  return(out/2)
} #for PC prior on phi
pc_log_pdf_phi<- function(phi,lambda,Q_inv,gamma_til){
  out <- log(lambda) - 0.5*log(8*kld_phi(phi,Q_inv)) - lambda*sqrt(2*kld_phi(phi,Q_inv)) +
    log(tr(Q_inv)-nrow(Q_inv)-sum((gamma_til-1)/(1+phi*(gamma_til-1))))
  return(out)
} #log-density of PC prior for phi


#######################################################################
#######################################################################
# Alg 2.1: MCMC to get P(r|Y) (fixed spatial effect)  -------
#target distribution
log_target_b_Y <- function(b,data,area){
  #Y|b
  dpois(sum(data[data$A==area,]$Y),sum(data[data$A==area,]$N)*exp(b[area]), log=T) + 
    #prior on b + jacobian
    dnorm(b[area],-6,10, log=T) - b[area]
}

#MCMC -- accepts and rejects each b independently (likelihoods are totally independent)
postsamp_b_Y_MCMC <- function(b_init,prop_sd,data,n_iter,burnin=n_iter/2,chains=2){
  n_areas <- length(b_init)
  b_postsamp_t <- rbind(b_init, matrix(NA,n_iter-1,n_areas))
  b_postsamp <- NULL
  for(chain in 1:chains){
    for(i in 2:n_iter){
      b_current <- b_postsamp_t[(i-1),]
      b_proposals <- as.vector(Rfast::rmvnorm(1,b_current,diag(prop_sd^2,n_areas)))
      a_prob <- sapply(1:n_areas,function(area){
        b_proposed <- b_current
        b_proposed[area] <- b_proposals[area]
        exp(log_target_b_Y(b_proposed,data,area) - log_target_b_Y(b_current,data,area) + 
              dnorm(b_current[area],b_proposed[area],prop_sd,log=T) - dnorm(b_proposed[area],b_current[area],prop_sd,log=T))})
      u <- runif(1)
      #determine whether we accept or reject each b
      accept <- u<a_prob
      b_current <- accept*b_proposals + (1-accept)*b_current
      b_postsamp_t[i,] <- b_current
    } 
    #get rid of burn-in and add to other chains
    b_postsamp <- rbind(b_postsamp,b_postsamp_t[(burnin+1):n_iter,])
  }
  return(list(b=b_postsamp))
}

# #MCMC -- accepts and rejects b's together
postsamp_b_Y_MCMC2 <- function(b_init,prop_sd,data,n_iter,burnin=n_iter/2,chains=2){
  n_areas <- length(b_init)
  b_postsamp_t <- rbind(b_init, matrix(NA,n_iter-1,n_areas))
  b_postsamp <- NULL
  for(chain in 1:chains){
    for(i in 2:n_iter){
      b_current <- b_postsamp_t[(i-1),]
      b_proposals <- as.vector(Rfast::rmvnorm(1,b_current,diag(prop_sd^2,n_areas)))
      a_prob <- exp(sum(sapply(1:n_areas,function(area){
        log_target_b_Y(b_proposals,data,area) - log_target_b_Y(b_current,data,area) +
          dnorm(b_current[area],b_proposals[area],prop_sd,log=T) - dnorm(b_proposals[area],b_current[area],prop_sd,log=T)})))
      u <- runif(1)
      #determine whether we accept or reject each b
      if(u<a_prob){
        b_current <- b_proposals
      }
      b_postsamp_t[i,] <- b_current
    } 
    #get rid of burn-in and add to other chains
    b_postsamp <- rbind(b_postsamp,b_postsamp_t[(burnin+1):n_iter,])
  }
  return(list(b=b_postsamp))
}

postsamp_test <- postsamp_b_Y_MCMC2(rep(-6,n_areas),0.05,dat1,2000)
plot(postsamp_test$b[,7],type='l')
plot(postsamp_test$b[,14],type='l')
plot(postsamp_test$b[,24],type='l')

#######################################################################
#######################################################################
# Alg 2.2: MCMC to get P(r|Ys,Y+) (fixed spatial effect) -------  
#target distribution
log_target_b_YYplus <- function(b,data,area){
  #Yplus|b
  dpois(sum(data$Y), sum(sapply(1:n_areas, function(i){sum(data[data$A==i,]$N)*exp(b[i])})),log=T) +
    #Y|b (only terms that don't drop out)
    sum(data[data$A==area,]$Y)*log(sum(data[data$A==area,]$N)*exp(b[area])) -
    sum(data$Y)*log(sum(sapply(1:n_areas, function(i){sum(data[data$A==i,]$N)*exp(b[i])}))) +
    #prior on b w/ Jacobian
    dnorm(b[area],-6,10, log=T) - b[area]
}

#MCMC -- Gibbs sampler for b's with MH to draw from each conditional
postsamp_b_YYplus_MCMC <- function(b_init,prop_sd,data,n_iter=2000,burnin=n_iter/2,chains=2){
  n_areas <- length(b_init)
  b_postsamp_t <- rbind(b_init, matrix(NA,n_iter-1,n_areas))
  b_postsamp <- NULL
  b_current <- b_init
  for(chain in 1:chains){
    for(i in 2:n_iter){
      #generate all proposals for MH (can do because proposal dist is IID)
      b_proposals <- as.vector(Rfast::rmvnorm(1,b_current,diag(prop_sd^2,n_areas)))
      #gibbs with internal MH
      b_current <- sapply(1:n_areas,function(j){
        #only change the b for this area
        b_proposed <- b_current
        b_proposed[j] <- b_proposals[j]
        
        a_prob <- exp(log_target_b_YYplus(b_proposed,data,j) - log_target_b_YYplus(b_current,data,j) +
                        dnorm(b_current[j],b_proposed[j],prop_sd,log=T) - dnorm(b_proposed[j],b_current[j],prop_sd,log=T))
        u <- runif(1)
        #determine whether we accept or reject each b
        if(u<a_prob){
          return(b_proposed[j])
        }else{
          return(b_current[j])
        }
      })
      b_postsamp_t[i,] <- b_current
    } 
    #get rid of burn-in and add to other chains
    b_postsamp <- rbind(b_postsamp,b_postsamp_t[(burnin+1):n_iter,])
  }
  return(list(b=b_postsamp))
}


#######################################################################
#######################################################################
# Fixed effect simulations (Alg 2.1 and 2.2) ------------

#generate data
n_areas <- 15
n_clusters <- rep(10,n_areas)
n_births <- 100
b <- rnorm(n_areas,-6,0.25)
nsims <- 10

#run simulations
results2.1_b <- matrix(NA,nsims,n_areas)
results2.2_b <- matrix(NA,nsims,n_areas)
for(k in 1:nsims){
  #status check
  if(k %% (nsims/10)==0){cat('Starting simulation', k, '\n')}
  
  #generate data
  dat1 <- simdat_fsp(n_areas,n_clusters,n_births,b)
  
  #run MCMC
  postsamp_2.1 <- postsamp_b_Y_MCMC(rep(-6,n_areas),0.5,dat1,2000)
  #postsamp_2.2 <- postsamp_b_YYplus_MCMC(rep(-6,n_areas),0.5,dat1,2000)
  
  #record result
  results2.1_b[k,] <- colMeans(postsamp_2.1$b)
  #results2.2_b[k,] <- colMeans(postsamp_2.2$b)
}

#organize and add true value for comparison
results2.1_b_long <- gather(data.frame(results2.1_b),'area','b')
results2.1_b_long$model <- 1
results2.2_b_long <- gather(data.frame(results2.2_b),'area','b')
results2.2_b_long$model <- 2
results2_b <- rbind(results2.1_b_long,results2.2_b_long)
results2_b$area <- as.numeric(str_remove(results2_b$area,'X'))
true_b <- data.frame(area=1:n_areas,true_b=b)
results2_b <- left_join(results2_b,true_b,by='area')

#plot
pdf("Handcoded MCMC Simulations/Fixed spatial effect 230726.pdf")
{
  results2_b %>% ggplot() + geom_boxplot(aes(x=factor(model),y=exp(b),col=factor(model))) + 
    geom_hline(aes(yintercept = exp(true_b))) + scale_color_discrete('Algorithm',labels=c('2.1','2.2')) +
    facet_wrap(~area) + ggtitle('10 simulated datasets of N=50000 per area, fixed spatial effect')
}
dev.off()

#######################################################################
#######################################################################
# Alg 2.3a: MCMC to get P(r,a,tau|Y) (IID random effect on b) -----------

# mean and precision of proposal distribution for vector r (calculated from Rue 2002)
prop_mean_iid_eta <- function(eta_star,tau,alpha,data){
  num <- sapply(1:n_areas,function(area){tau*alpha + data[area,]$Y - data[area,]$N*(1-eta_star[area])*exp(eta_star[area])})
  denom <- tau + sapply(1:n_areas,function(area){data[area,]$N*exp(eta_star[area])})
  return(num/denom)
}
prop_prec_iid_eta <- function(eta_star,tau,data){
  out <-  tau + sapply(1:n_areas,function(area){data[area,]$N*exp(eta_star[area])})
  return(out)
}

#MCMC -- block updates (hyperparameters, then parameters) 
## consistently underestimating alpha by a few hundredths, but works for some reason when I get rid of the Jacobian (must be wrong..)?
postsamp_artau_Y_MCMC <- function(tau_init,alpha_init, #initial values
                                  data,n_iter,burnin=n_iter/2,chains=2){
  n_areas <- length(unique(data$A))
  #collapse data into areas only
  newdat <- data.frame(data %>% group_by(A) %>% summarise(Y=sum(Y),N=sum(N)))
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  alpha_postsamp_t <- c(alpha_init, rep(NA,n_iter-1))
  eta_init <- rnorm(n_areas,0,sqrt(1/tau_init)) + alpha_init
  eta_postsamp_t <- rbind(eta_init, matrix(NA,n_iter-1,n_areas))
  tau_postsamp <- alpha_postsamp <- b_postsamp <- eta_postsamp <- r_postsamp <- NULL
  for(chain in 1:chains){
    for(i in 2:n_iter){
      tau_current <- tau_postsamp_t[i-1]
      alpha_current <- alpha_postsamp_t[i-1]
      eta_current <- eta_postsamp_t[(i-1),]
      
      # Draw new tau from full conditional
      tau_current <- rgamma(1,0.01+n_areas/2,0.01+0.5*sum((eta_current - alpha_current)^2))
      
      # Draw new alpha from full conditional
      alpha_current <- rnorm(1,(-6+3^2*tau_current*sum(eta_current))/(1+3^2*n_areas*tau_current),
                             1/sqrt(n_areas*tau_current + 1/(3^2)))
      
      #Draw proposal for eta for Taylor approximation of full conditional
      eta_proposals <- as.vector(Rfast::rmvnorm(1, prop_mean_iid_eta(eta_current,tau_current,alpha_current,newdat),
                                                diag(1/prop_prec_iid_eta(eta_current,tau_current,newdat))))
      
      a_prob_eta <- exp(sum(sapply(1:n_areas, function(area){dpois(newdat[area,]$Y,newdat[area,]$N*exp(eta_proposals[area]),log=T)})) - #target likelihood
                          sum(sapply(1:n_areas, function(area){dpois(newdat[area,]$Y,newdat[area,]$N*exp(eta_current[area]),log=T)})) +
                          Rfast::dmvnorm(eta_proposals - alpha_current,rep(0,n_areas),diag(1/tau_current,n_areas),logged = T) - #target prior on eta
                          Rfast::dmvnorm(eta_current - alpha_current,rep(0,n_areas),diag(1/tau_current,n_areas),logged = T) + 
                          Rfast::dmvnorm(eta_current, prop_mean_iid_eta(eta_proposals,tau_current,alpha_current,newdat),
                                         diag(1/prop_prec_iid_eta(eta_proposals,tau_current,newdat)), logged=T) - 
                          Rfast::dmvnorm(eta_proposals, prop_mean_iid_eta(eta_current,tau_current,alpha_current,newdat),
                                         diag(1/prop_prec_iid_eta(eta_current,tau_current,newdat)), logged=T))
      
      if(runif(1) < a_prob_eta){
        eta_current <- eta_proposals}
      
      # record state
      tau_postsamp_t[i] <- tau_current
      alpha_postsamp_t[i] <- alpha_current
      eta_postsamp_t[i,] <- eta_current
    } 
    #get rid of burn-in and add to other chains
    eta_postsamp <- rbind(eta_postsamp,eta_postsamp_t[(burnin+1):n_iter,])
    alpha_postsamp <- c(alpha_postsamp,alpha_postsamp_t[(burnin+1):n_iter])
    tau_postsamp <- c(tau_postsamp,tau_postsamp_t[(burnin+1):n_iter])
  }
  
  b_postsamp <- eta_postsamp-alpha_postsamp
  r_postsamp <- exp(eta_postsamp)
  return(list(alpha = alpha_postsamp,b=b_postsamp,eta=eta_postsamp,r=r_postsamp,tau=tau_postsamp))
}

#######################################################################
#######################################################################
# IID random effect simulations (Alg 2.3a) ------------
n_areas <- 50
n_clusters <- rep(10,n_areas)
n_births <- 10
tau <- rgamma(1,10,0.25)
alpha <- rnorm(1,-6,0.25)
b <- Rfast::rmvnorm(1,rep(0,n_areas),diag(1/tau,n_areas))

### get posterior distributions for 1 realization
data <- simdat_iid(n_areas,n_clusters,n_births,alpha,tau,b)$dat
postsamp_2.3a <- postsamp_artau_Y_MCMC(tau_init = 40, alpha_init = -6, data,5000,chains=2)

pdf("Handcoded MCMC Simulations/IID spatial effect 230801, posterior distributions.pdf") 
{
  hist(log(postsamp_2.3a$tau))
  abline(v=log(tau),col='red')
  hist(postsamp_2.3a$alpha)
  abline(v=alpha,col='red')
  
  results2.3a_b <- gather(data.frame(postsamp_2.3a$b),'area','b')
  results2.3a_b$area <- as.numeric(str_remove(results2.3a_b$area,'X'))
  true_b <- data.frame(area=1:n_areas,true_b=as.vector(b))
  results2.3a_b <- left_join(results2.3a_b,true_b,by='area')
  g <- results2.3a_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~area) + ggtitle('10000 births per area')
  print(g)
  
  plot(log(postsamp_2.3a$tau),type='l')
  abline(h=log(tau),col='red')
  plot(postsamp_2.3a$alpha,type='l')
  abline(h=alpha,col='red')
  plot(postsamp_2.3a$b[,6],type='l')
  abline(h=b[6],col='red')
  plot(postsamp_2.3a$b[,15],type='l')
  abline(h=b[15],col='red')
  plot(postsamp_2.3a$b[,42],type='l')
  abline(h=b[42],col='red')
} 
dev.off()  


### get distribution of posterior mean for several realizations
nsims <- 25
results2.3a_tau <- rep(NA,nsims)
results2.3a_alpha <- rep(NA,nsims)
results2.3a_eta <- matrix(NA,nsims,n_areas)
for(k in 1:nsims){
  #status check
  if(k %% (nsims/5)==0){cat('Starting simulation', k, '\n')}
  
  #generate data
  data <- simdat_iid(n_areas,n_clusters,n_births,alpha,tau,b)$dat
  
  #run MCMC
  postsamp_2.3a <- postsamp_artau_Y_MCMC(tau_init = 40, alpha_init = -6, 
                                         data,n_iter=5000,chains=2)
  
  #record result
  results2.3a_tau[k] <- mean(postsamp_2.3a$tau)
  results2.3a_alpha[k] <- mean(postsamp_2.3a$alpha)
  results2.3a_eta[k,] <- colMeans(postsamp_2.3a$eta)
}

pdf("Handcoded MCMC Simulations/IID spatial effect 230731, same b.pdf")
{
  hist(log(results2.3a_tau))
  abline(v=log(tau),col='red')
  hist(results2.3a_alpha)
  abline(v=alpha,col='red')
  
  results_eta <- gather(data.frame(results2.3a_eta),'area','eta')
  results_eta$area <- as.numeric(str_remove(results_eta$area,'X'))
  true_eta <- data.frame(area=1:n_areas,true_eta=alpha + as.vector(b))
  results_eta <- left_join(results_eta,true_eta,by='area')
  results_eta %>% ggplot() + geom_boxplot(aes(y=eta)) +
    geom_hline(aes(yintercept = true_eta),col='red') + facet_grid(~area) + ggtitle('Posterior means from 25 simulated data sets with same b')
}

dev.off()

#######################################################################
#######################################################################
# Alg 2.3: MCMC to get P(r,a,tau,phi|Y) (BYM2 spatial effect)  -------

# updates all parameters except alpha in 1 block
postsamp_artauphi_Y_MCMC <- function(tau_init,phi_init,alpha_init, #initial values
                                     prop_sd_phi,prop_sd_tau,
                                     data,Q_inv,n_iter,burnin=n_iter/2,chains=2){
  n_areas <- length(unique(data$A))
  #collapse data into areas only
  newdat <- data.frame(data %>% group_by(A) %>% summarise(Y=sum(Y),N=sum(N)))
  #hyperpriors
  a_tau <- 0.01
  b_tau <-  0.01
  a_phi <- 2
  b_phi <- 3
  mu_a <- -6
  sd_a <- 3
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  phi_postsamp_t <- c(phi_init, rep(NA,n_iter-1))
  alpha_postsamp_t <- c(alpha_init, rep(NA,n_iter-1))
  acc <- 0
  
  eta_init <- Rfast::rmvnorm(1,rep(alpha_init,n_areas),1/tau_init*(W(phi_init,Q_inv))) 
  eta_postsamp_t <- rbind(eta_init, matrix(NA,n_iter-1,n_areas))
  tau_postsamp <- phi_postsamp <- alpha_postsamp <- b_postsamp <- eta_postsamp <- r_postsamp <- NULL
  for(chain in 1:chains){
    for(i in 2:n_iter){
      tau_current <- tau_postsamp_t[i-1]
      phi_current <- phi_postsamp_t[i-1]
      alpha_current <- alpha_postsamp_t[i-1]
      eta_current <- eta_postsamp_t[(i-1),]
      
      # Draw proposal for phi
      phi_proposal <- expit(rnorm(1,logitlink(phi_current),prop_sd_phi)) 
      
      # Draw proposal for tau
      tau_proposal <- rlnorm(1,log(tau_current),prop_sd_tau)
      
      # Draw proposal for eta for Taylor approximation of full conditional
      eta_proposal <- as.vector(Rfast::rmvnorm(1,prop_mean_bym2_eta(eta_current,tau_proposal,phi_proposal,alpha_current,newdat,Q_inv),
                                               prop_sigma_bym2_eta(eta_current,tau_proposal,phi_proposal,newdat,Q_inv)))
      
      aprob <- exp( #Y|eta
        sum(sapply(1:n_areas,function(area){
          dpois(newdat$Y[area],newdat$N[area]*exp(eta_proposal[area]),log=T) -
            dpois(newdat$Y[area],newdat$N[area]*exp(eta_current[area]),log=T)})) +
          #eta-alpha|phi,tau
          Rfast::dmvnorm(eta_proposal,rep(alpha_current,n_areas),sigma = 1/tau_proposal*(W(phi_proposal,Q_inv)), logged=T) -
          Rfast::dmvnorm(eta_current,rep(alpha_current,n_areas),sigma = 1/tau_current*(W(phi_current,Q_inv)), logged=T) +
          #prior on tau 
          dgamma(tau_proposal,a_tau,b_tau,log = T) - dgamma(tau_current,a_tau,b_tau,log = T) + 
          #prior on phi
          dbeta(phi_proposal,a_phi,b_phi,log=T) - dbeta(phi_current,a_phi,b_phi,log=T) +
          #Jacobian to go from logitphi to phi
          (-log(phi_current*(1-phi_current)) + log(phi_proposal*(1-phi_proposal))) + 
          #proposal distribution for logitphi is symmetric
          #proposal dist for tau
          dlnorm(tau_current,log(tau_proposal),prop_sd_tau,log=T) - 
          dlnorm(tau_proposal,log(tau_current),prop_sd_tau,log=T) +
          #proposal dist for alpha is symmetric
          #proposal dist for eta
          Rfast::dmvnorm(eta_current,prop_mean_bym2_eta(eta_proposal,tau_current,phi_current,alpha_current,newdat,Q_inv),
                         prop_sigma_bym2_eta(eta_proposal,tau_current,phi_current,newdat,Q_inv),logged=T) - 
          Rfast::dmvnorm(eta_proposal,prop_mean_bym2_eta(eta_current,tau_proposal,phi_proposal,alpha_current,newdat,Q_inv),
                         prop_sigma_bym2_eta(eta_current,tau_proposal,phi_proposal,newdat,Q_inv),logged=T)
      )
      if(runif(1) < aprob){
        phi_current <- phi_proposal
        tau_current <- tau_proposal
        eta_current <- eta_proposal
        if(i>burnin){
          acc <- acc+1 
        }
      }
      
      #draw alpha from full conditional
      alpha_current <- rnorm(1,(mu_a+sd_a^2*tau_current*sum(t(eta_current)%*%solve(W(phi_current,Q_inv))))/(1+sd_a^2*tau_current*sum(solve(W(phi_current,Q_inv)))),
                             1/sqrt((1/sd_a^2) + tau_current*sum(solve(W(phi_current,Q_inv)))))
      
      # record state
      tau_postsamp_t[i] <- tau_current
      phi_postsamp_t[i] <- phi_current
      alpha_postsamp_t[i] <- alpha_current
      eta_postsamp_t[i,] <- eta_current
      
    } 
    #get rid of burn-in and add to other chains
    eta_postsamp <- rbind(eta_postsamp,eta_postsamp_t[(burnin+1):n_iter,])
    alpha_postsamp <- c(alpha_postsamp,alpha_postsamp_t[(burnin+1):n_iter])
    tau_postsamp <- c(tau_postsamp,tau_postsamp_t[(burnin+1):n_iter])
    phi_postsamp <- c(phi_postsamp,phi_postsamp_t[(burnin+1):n_iter])
  }
  
  b_postsamp <- eta_postsamp-alpha_postsamp
  r_postsamp <- exp(eta_postsamp)
  acc_rate <- acc/((n_iter-burnin)*chains)
  
  return(list(alpha = alpha_postsamp,b=b_postsamp, eta=eta_postsamp,r=r_postsamp, tau=tau_postsamp,phi=phi_postsamp,acc_rate=acc_rate))
}

#######################################################################
#######################################################################
# BYM2 random effect simulations (Alg 2.3) ------------

## define ground truth
{
  n_areas <- 50
  n_clusters <- rep(10,n_areas)
  n_births <- 500
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
  Q_scaled <- inla.scale.model(Q, constr=list(A=t(eigen(Q)$vectors[,eigen(Q)$values<1e-10]), e=rep(0,sum(eigen(Q)$values<1e-10))))
  Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))
  #get approximation of eigenvalues for Q_scaled_inv
  gamma <- eigen(Q_scaled)$values
  gamma_til <- 1/gamma
  gamma_til[gamma_til>1e+10] <- 0
  
  # intercept (fixed effect) 
  alpha <- rnorm(1,-6,0.25)
  
  # draw hyperparameters
  tau <- rgamma(1,10,.25)
  phi <- rbeta(1,2,3)
}
## generate data 
{
  Var_b <- (1/tau)*W(phi,Q_scaled_inv)
  b <- Rfast::rmvnorm(1,rep(0,n_areas),Var_b)
  dat2 <- simdat_bym2sp(n_areas=n_areas,n_clusters=n_clusters,n_births=n_births,a=alpha,Var_b=Var_b,b=b)$dat
}
## parameters for PC prior on phi
{
  alpha_phi <- 2/3
  U_phi <- 0.5
  
  ## find lambda to properly calibrate PC prior on phi
  #use prior from INLA to calibrate to proper lambda for given alpha and u
  pc_pdf_phi_inla <- INLA:::inla.pc.bym.phi(Q=Q_scaled,rankdef = 1,alpha = alpha_phi,u=U_phi)
  prior_diff <- function(lambda,Q_scaled_inv,gamma_til){
    phi_seq <- seq(0.01,0.99,0.01)
    sum((sapply(phi_seq,function(x){pc_pdf_phi_inla(x)}) - sapply(phi_seq,function(x){pc_log_pdf_phi(x,lambda,Q_scaled_inv,gamma_til)}))^2)
  }
  
  lambda_opt <- optim(1,prior_diff,Q_scaled_inv=Q_scaled_inv,gamma_til=gamma_til,method = 'Brent',lower = 0,upper=10)$par
}

## run MCMC for one dataset -- good mixing and estimation after 50k iterations
postsamp_2.3 <- postsamp_artauphi_Y_MCMC(50,0.15,-6,
                                         prop_sd_phi=.25,prop_sd_tau = .25,
                                         dat2,Q_inv = Q_scaled_inv,50000,chains=2)

pdf("Handcoded MCMC Simulations/BYM2 spatial effect 230825, block updates except alpha, 50k iterations posterior distributions.pdf") 
{
  hist(log(postsamp_2.3$tau))
  abline(v=log(tau),col='red')
  hist(logitlink(postsamp_2.3$phi))
  abline(v=logitlink(phi),col='red')
  hist(postsamp_2.3$alpha)
  abline(v=alpha,col='red')
  
  results2.3_b <- gather(data.frame(postsamp_2.3$b),'area','b')
  results2.3_b$area <- as.numeric(str_remove(results2.3_b$area,'X'))
  true_b <- data.frame(area=1:n_areas,true_b=as.vector(b))
  true_b <- true_b[order(true_b$true_b),]
  true_b$b_rank <- 1:nrow(true_b)
  results2.3_b <- left_join(results2.3_b,true_b,by='area')
  g <- results2.3_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~b_rank) + ggtitle('5000 births per area') +
    geom_hline(yintercept=0,col='grey50')
  print(g)
  
  plot(log(postsamp_2.3$tau[1:50000%%50==0]),type='l')
  abline(h=log(tau),col='red')
  plot(logitlink(postsamp_2.3$phi[1:50000%%50==0]),type='l')
  abline(h=logitlink(phi),col='red')
  plot(postsamp_2.3$alpha[1:50000%%50==0],type='l')
  abline(h=alpha,col='red')
  plot(postsamp_2.3$b[1:50000%%50==0,7],type='l')
  abline(h=b[7],col='red')
  plot(postsamp_2.3$b[1:50000%%50==0,15],type='l')
  abline(h=b[15],col='red')
  plot(postsamp_2.3$b[1:50000%%50==0,28],type='l')
  abline(h=b[28],col='red')
  plot(postsamp_2.3$b[1:50000%%50==0,42],type='l')
  abline(h=b[42],col='red')
} 
dev.off()  

# check against INLA with all PC priors
hyper.inla <- list(prec = list(prior = "loggamma", param = c(0.01,0.01)),
                   phi = list(prior= "logitbeta", param=c(2,3)))

inlafit <- INLA::inla(Y ~ f(A, graph = abs(Q),
                            model = "bym2",
                            hyper = hyper.inla,
                            scale.model = TRUE,
                            adjust.for.con.comp = TRUE),
                      data=dat2, family='poisson', E=N,
                      control.predictor = list(compute = FALSE, link = 1))

### get distribution of posterior mean for several data realizations
{
  nsims <- 50
  results2.3_tau <- rep(NA,nsims)
  results2.3_phi <- rep(NA,nsims)
  results2.3_alpha <- rep(NA,nsims)
  results2.3_eta <- matrix(NA,nsims,n_areas)
  for(k in 1:nsims){
    #status check
    if(k %% (nsims/10)==0){cat('Starting simulation', k, '\n')}
    
    #generate data
    dat2 <- simdat_bym2sp(n_areas=n_areas,n_clusters=n_clusters,n_births=n_births,
                          a=alpha,Var_b=Var_b,b=b)$dat
    
    #run MCMC
    postsamp_2.3 <- postsamp_artauphi_Y_MCMC(50,0.25,-6,
                                             prop_sd_phi=0.15, prop_sd_tau = 0.1, prop_sd_alpha = 0.15,
                                             dat2,Q_inv = Q_scaled_inv,50000,chains=2)
    
    #record result
    results2.3_tau[k] <- mean(postsamp_2.3$tau)
    results2.3_phi[k] <- mean(postsamp_2.3$phi)
    results2.3_alpha[k] <- mean(postsamp_2.3$alpha)
    results2.3_eta[k,] <- colMeans(postsamp_2.3$eta)
  }
  
}

pdf("Handcoded MCMC Simulations/BYM2 spatial effect 230824 same b, block updates 50k iterations.pdf")
{
  hist(log(results2.3_tau))
  abline(v=log(tau),col='red')
  hist(logitlink(results2.3_phi))
  abline(v=logitlink(phi),col='red')
  hist(results2.3_alpha)
  abline(v=alpha,col='red')
  
  results_eta <- gather(data.frame(results2.3_eta),'area','eta')
  results_eta$area <- as.numeric(str_remove(results_eta$area,'X'))
  results_eta$b <- results_eta$eta - results2.3_alpha
  true_b <- data.frame(area=1:n_areas,true_b=as.vector(b))
  results_b <- left_join(results_eta,true_b,by='area')
  results_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~area) + ggtitle('Posterior means from 50 simulated data sets with same b')
}
dev.off()



#######################################################################
#######################################################################
# Alg 2.5: P(r,alphaU,alphaR,b,tau,d|Y) --------

#gradient of log target wrt b (for MALA)-- not sure this is right
grad_log_target_resp_b <- function(b,alphaU,alphaR,tau,d,data){
  Nr <- sapply(1:nrow(data),function(i){data$N[i]*exp(alphaU*data$U[i]+alphaR*(1-data$U[i])+b[data$A[i]])})
  grad_log_target_tt <- sapply(1:nrow(data), function(i){
    1/d*Nr[i]*(digamma(data$Y[i]+1/d*Nr[i]) - digamma(1/d*Nr[i]) - 1/(1+d) - log(1+d))})
  return(sapply(1:n_areas, function(area){
    sum(grad_log_target_tt[which(data$A==area)]) - tau*b[area]
  }))
}

tau_init <- 25
d_init <- 0.1
alphaU_init <- -6
alphaR_init <- -5
prop_sd_d <- 0.15
prop_k_b <- 0.01
prop_sd_alpha <- 0.1
data <- all_dat
Q_inv <- Q_scaled_inv
n_iter <- 1000
burnin <- 500
postsamp_Alg2.5_MCMC <- function(tau_init, d_init, alphaU_init, alphaR_init, #initial values
                                 prop_sd_d, prop_sd_alpha, # for random walk proposals
                                 prop_sd_b,
                                 data,Q_inv,n_iter,burnin=n_iter/2,chains=2){
  n_areas <- length(unique(data$A))
  #hyperpriors
  a_tau <- a_sig2 <- 0.01
  b_tau <- b_sig2 <-  0.01
  mu_U <- -6
  sd_U <- 3
  mu_R <- -5
  sd_R <- 3
  scale_d <- 5
  
  #set initial values
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  d_postsamp_t <- c(d_init, rep(NA,n_iter-1))
  alphaU_postsamp_t <- c(alphaU_init, rep(NA,n_iter-1))
  alphaR_postsamp_t <- c(alphaR_init, rep(NA,n_iter-1))
  
  b_init <- rnorm(n_areas,0,1/sqrt(tau_init))
  b_postsamp_t <- rbind(b_init, matrix(NA,n_iter-1,n_areas))
  
  tau_postsamp <- d_postsamp <- alphaU_postsamp <- alphaR_postsamp <- b_postsamp <- NULL
  acc <- 0
  
  # set current states
  tau_current <- tau_init
  d_current <- d_init
  alphaU_current <- alphaU_init
  alphaR_current <- alphaR_init
  b_current <- b_init
  
  for(chain in 1:chains){
    for(i in 2:n_iter){
      
      # draw new tau from full conditional
      tau_current <- rgamma(1,n_areas/2+a_tau,0.5*sum(b_current^2)+b_tau)
      
      # draw proposal for b 
      # random walk
      b_proposal <- sapply(b_current, function(x){(rnorm(1,x,prop_sd_b))})
      # MALA
      ## b_proposal <- b_current + prop_k_b*grad_log_target_resp_b(b_current,alphaU_current,alphaR_current,tau_current,d_current,data) + sqrt(2*prop_k_b)*rnorm(n_areas)
      
      # draw proposal for d
      d_proposal <- rlnorm(1,log(d_current),prop_sd_d)
      
      # draw proposal for alphas
      alphaU_proposal <- rnorm(1,alphaU_current,prop_sd_alpha)
      alphaR_proposal <- rnorm(1,alphaR_current,prop_sd_alpha)
      
      # accept or reject proposal for (alpha, b, d)
      a_prob <- exp({
        # Y|alpha,b,d
        sum(vapply(1:nrow(data),function(i){
          dnbinom(data$Y[i], size=1/d_proposal*data$N[i]*exp(alphaU_proposal*data$U[i]+alphaR_proposal*(1-data$U[i])+b_proposal[data$A[i]]),
                  prob=1/(1+d_proposal),log=T) - 
            dnbinom(data$Y[i], size=1/d_current*data$N[i]*exp(alphaU_current*data$U[i]+alphaR_current*(1-data$U[i])+b_current[data$A[i]]),
                    prob=1/(1+d_current),log=T)
        },numeric(1))) +
          # b|tau
          sum(vapply(1:n_areas,function(i){
            dnorm(b_proposal[i],0,sqrt(1/tau_current),log=T) - dnorm(b_current[i],0,sqrt(1/tau_current),log=T)
          },numeric(1))) +
          # prior on alphas
          dnorm(alphaU_proposal,mu_U,sd_U,log=T) - dnorm(alphaU_current,mu_U,sd_U,log=T) +
          dnorm(alphaR_proposal,mu_R,sd_R,log=T) - dnorm(alphaR_current,mu_R,sd_R,log=T) +
          # prior on d
          dhalfnorm(d_proposal,scale_d,log=T) - dhalfnorm(d_current,scale_d,log=T) +
          # proposal dists for alpha are symmetric
          # prop dist for b (MALA)
          #Rfast::dmvnorm(b_current,b_proposal + prop_k_b*grad_log_target_resp_b(b_proposal,alphaU_current,alphaR_current,tau_current,d_current,data), diag(2*prop_k_b,n_areas), logged=T) - 
          #Rfast::dmvnorm(b_proposal,b_current + prop_k_b*grad_log_target_resp_b(b_current,alphaU_current,alphaR_current,tau_current,d_current,data), diag(2*prop_k_b,n_areas), logged=T) +
          # proposal dist for d
          dlnorm(d_current,log(d_proposal),prop_sd_d,log=T) - dlnorm(d_proposal,log(d_current),prop_sd_d,log=T)
      })
      
      if(a_prob>runif(1)){
        d_current <- d_proposal
        b_current <- b_proposal
        alphaU_current <- alphaU_proposal
        alphaR_current <- alphaR_proposal
        acc <- acc + 1
      }
      
      # record current state
      tau_postsamp_t[i] <- tau_current
      d_postsamp_t[i] <- d_current
      alphaU_postsamp_t[i] <- alphaU_current
      alphaR_postsamp_t[i] <- alphaR_current
      b_postsamp_t[i,] <- b_current
    } 
    
    #get rid of burn-in and add to other chains
    b_postsamp <- rbind(b_postsamp,b_postsamp_t[(burnin+1):n_iter,])
    alphaU_postsamp <- c(alphaU_postsamp,alphaU_postsamp_t[(burnin+1):n_iter])
    alphaR_postsamp <- c(alphaR_postsamp,alphaR_postsamp_t[(burnin+1):n_iter])
    tau_postsamp <- c(tau_postsamp,tau_postsamp_t[(burnin+1):n_iter])
    d_postsamp <- c(d_postsamp,d_postsamp_t[(burnin+1):n_iter])
  }
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp, alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  acc_rate <- acc/((n_iter-burnin)*chains)
  
  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              tau=tau_postsamp, d=d_postsamp, 
              acc_rate=acc_rate))
}
#######################################################################
#######################################################################
# Alg 2.6: P(alphaU,alphaR,b,tau,d|Y) --------

#gradients for MALA
grad_log_target_2.6_alphaU <- function(data,alphaU,alphaR,b,d,mu_U,sd_U){
  Nr_sum_urb <- sum(data$N*data$U*exp(alphaU+b[data$A]))
  Nr_sum_rur <- sum(data$N*(1-data$U)*exp(alphaR+b[data$A]))
  Nr_sum <- Nr_sum_urb + Nr_sum_rur
  term1 <- sum(data$Y*data$U)
  term2 <- (1/d*(digamma(sum(data$Y)+1/d*Nr_sum)-digamma(1/d*Nr_sum)-log(1+d)) - sum(data$Y)/Nr_sum)*Nr_sum_urb #negbin
  term3 <- -1/sd_U^2*(alphaU-mu_U) #prior
  return(term1+term2+term3)
}

grad_log_target_2.6_alphaR <- function(data,alphaU,alphaR,b,d,mu_R,sd_R){
  Nr_sum_urb <- sum(data$N*data$U*exp(alphaU+b[data$A]))
  Nr_sum_rur <- sum(data$N*(1-data$U)*exp(alphaR+b[data$A]))
  Nr_sum <- Nr_sum_urb + Nr_sum_rur
  term1 <- sum(data$Y*(1-data$U))
  term2 <- (1/d*(digamma(sum(data$Y)+1/d*Nr_sum)-digamma(1/d*Nr_sum)-log(1+d)) - sum(data$Y)/Nr_sum)*Nr_sum_rur #negbin
  term3 <- -1/sd_R^2*(alphaR-mu_R) #prior
  return(term1+term2+term3)
}

grad_log_target_2.6_logd <- function(data,alphaU,alphaR,b,d,mu_d,sd_d){
  Nr_sum <- sum(data$N*exp(alphaU*data$U + alphaR*(1-data$U) + b[data$A]))
  term1 <- (sum(data$Y) - Nr_sum)/(1+d)
  term2 <- (1/d)*Nr_sum*(log(1+d) - digamma(sum(data$Y) + 1/d*Nr_sum) + digamma(1/d*Nr_sum))
  term3 <- -(log(d)-mu_d)/sd_d^2
  return(term1+term2+term3)
}

grad_log_target_2.6_b <- function(data,area,alphaU,alphaR,b,d,tau){
  subdat <- data[data$A==area,]
  Nr_sum_area <- sum(subdat$N*exp(alphaU*subdat$U + alphaR*(1-subdat$U) + b[area]))
  Nr_sum <- sum(data$N*exp(alphaU*data$U + alphaR*(1-data$U) + b[data$A]))
  term1 <- sum(subdat$Y)
  term2 <- Nr_sum_area*(1/d*(digamma(sum(data$Y)+1/d*Nr_sum) - digamma(1/d*Nr_sum) - log(1+d)) - sum(data$Y)/Nr_sum)
  term3 <- -tau*b[area]
  return(term1+term2+term3)
}

tau_init <- 25
d_init <- 0.1
alphaU_init <- -3
alphaR_init <- -3
prop_sd_d <- 0.15
prop_k_b <- 1e-2
prop_sd_b <- 0.05
prop_sd_alpha <- 0.1
prop_k_alpha <- 2e-3
prop_k_logd <- 0.05
data <- all_dat
n_iter <- 1000
burnin <- 500

postsamp_Alg2.6_MCMC <- function(tau_init, d_init, alphaU_init, alphaR_init, #initial values
                                 prop_sd_d=NULL, prop_sd_alpha=NULL, prop_sd_b=NULL, # for random walk proposals
                                 prop_k_alpha=NULL, prop_k_logd=NULL, prop_k_b=NULL, # for MALA proposals
                                 data,n_iter,burnin=n_iter/2,chains=2){
  n_areas <- length(unique(data$A))
  Yplus <- sum(data$Y)
  #hyperpriors
  a_tau <- a_sig2 <- 0.01
  b_tau <- b_sig2 <-  0.01
  mu_U <- -3
  sd_U <- 3
  mu_R <- -3
  sd_R <- 3
  #scale_d <- 5
  mu_d <- 0
  sd_d <- 1
  
  #set initial values
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  d_postsamp_t <- c(d_init, rep(NA,n_iter-1))
  alphaU_postsamp_t <- c(alphaU_init, rep(NA,n_iter-1))
  alphaR_postsamp_t <- c(alphaR_init, rep(NA,n_iter-1))
  
  b_init <- rnorm(n_areas,0,1/sqrt(tau_init))
  b_postsamp_t <- rbind(b_init, matrix(NA,n_iter-1,n_areas))
  
  tau_postsamp <- d_postsamp <- alphaU_postsamp <- alphaR_postsamp <- b_postsamp <- NULL
  acc <- 0
  
  # set current states
  tau_current <- tau_init
  d_current <- d_init
  alphaU_current <- alphaU_init
  alphaR_current <- alphaR_init
  b_current <- b_init
  
  for(chain in 1:chains){
    for(i in 2:n_iter){
      
      # draw new tau from full conditional
      tau_current <- rgamma(1,n_areas/2+a_tau,0.5*sum(b_current^2)+b_tau)
      
      # draw proposal for b 
      # random walk
       b_proposal <- sapply(b_current, function(x){(rnorm(1,x,prop_sd_b))})

      # b_proposal <- b_current +
      #   sapply(1:n_areas,function(A){prop_k_b*grad_log_target_2.6_b(data,A,alphaU_current,alphaR_current,b_current,d_current,tau_current)}) +
      #   sqrt(2*prop_k_b)*rnorm(n_areas)
      
      # draw proposal for d
      d_proposal <- rlnorm(1,log(d_current),prop_sd_d)
      # logd_proposal <- log(d_current) + prop_k_logd*grad_log_target_2.6_logd(data,alphaU_current,alphaR_current,b_current,d_current,mu_d,sd_d) +
      #   sqrt(2*prop_k_logd)*rnorm(1)
      # d_proposal <- exp(logd_proposal)
      
      # draw proposal for alphas
      alphaU_proposal <- rnorm(1,alphaU_current,prop_sd_alpha)
      alphaR_proposal <- rnorm(1,alphaR_current,prop_sd_alpha)
      # alphaU_proposal <- alphaU_current + prop_k_alpha*grad_log_target_2.6_alphaU(data,alphaU_current,alphaR_current,b_current,d_current,mu_U,sd_U) +
      #   sqrt(2*prop_k_alpha)*rnorm(1)
      # alphaR_proposal <- alphaR_current + prop_k_alpha*grad_log_target_2.6_alphaR(data,alphaU_current,alphaR_current,b_current,d_current,mu_R,sd_R) +
      #   sqrt(2*prop_k_alpha)*rnorm(1)
      
      #calculate the alphas parameter for the Dir-Multinom distributions
      alpha_dmc_proposal <- sapply(1:nrow(data),function(i){data$N[i]*exp(alphaU_proposal*data$U[i]+alphaR_proposal*(1-data$U[i])+b_proposal[data$A[i]])/d_proposal})
      alpha_dmc_current <- sapply(1:nrow(data),function(i){data$N[i]*exp(alphaU_current*data$U[i]+alphaR_current*(1-data$U[i])+b_current[data$A[i]])/d_current})
      
      # accept or reject proposal for (alpha, b, d)
      a_prob <- exp({
        # Y|Y+,alpha,b
        lgamma(sum(alpha_dmc_proposal)) - lgamma(Yplus + sum(alpha_dmc_proposal)) + sum(lgamma(data$Y+alpha_dmc_proposal)) -  sum(lgamma(alpha_dmc_proposal)) - 
          (lgamma(sum(alpha_dmc_current)) - lgamma(Yplus + sum(alpha_dmc_current)) + sum(lgamma(data$Y+alpha_dmc_current)) -  sum(lgamma(alpha_dmc_current))) +
          # Yplus|alpha,b,d
          dnbinom(Yplus, size = sum(alpha_dmc_proposal), prob=1/(1+d_proposal),log=T) - 
          dnbinom(Yplus, size= sum(alpha_dmc_current), prob=1/(1+d_current),log=T) +
          # b|tau
          sum(sapply(1:n_areas,function(i){
            dnorm(b_proposal[i],0,sqrt(1/tau_current),log=T) - dnorm(b_current[i],0,sqrt(1/tau_current),log=T)})) +
          # prior on alphas
          dnorm(alphaU_proposal,mu_U,sd_U,log=T) - dnorm(alphaU_current,mu_U,sd_U,log=T) +
          dnorm(alphaR_proposal,mu_R,sd_R,log=T) - dnorm(alphaR_current,mu_R,sd_R,log=T) +
          # prior on d
          # dhalfnorm(d_proposal,scale_d,log=T) - dhalfnorm(d_current,scale_d,log=T) +
          dlnorm(d_proposal,mu_d,sd_d,log=T) - dlnorm(d_current,mu_d,sd_d,log=T) +
          # # proposal dist for alpha (MALA)
          # dnorm(alphaU_current,alphaU_proposal + prop_k_alpha*grad_log_target_2.6_alphaU(data,alphaU_proposal,alphaR_current,b_current,d_current,mu_U,sd_U),sqrt(2*prop_k_alpha),log=T) -
          # dnorm(alphaU_proposal,alphaU_current + prop_k_alpha*grad_log_target_2.6_alphaU(data,alphaU_current,alphaR_current,b_current,d_current,mu_U,sd_U),sqrt(2*prop_k_alpha),log=T) +
          # dnorm(alphaR_current,alphaR_proposal + prop_k_alpha*grad_log_target_2.6_alphaR(data,alphaU_current,alphaR_proposal,b_current,d_current,mu_R,sd_R),sqrt(2*prop_k_alpha),log=T) -
          # dnorm(alphaR_proposal,alphaR_current + prop_k_alpha*grad_log_target_2.6_alphaR(data,alphaU_current,alphaR_current,b_current,d_current,mu_R,sd_R),sqrt(2*prop_k_alpha),log=T) +
          # # MALA proposal dist for b
          # sum(vapply(1:n_areas,function(i){
          #   dnorm(b_current[i],b_proposal[i] + prop_k_b*grad_log_target_2.6_b(data,i,alphaU_current,alphaR_current,b_proposal,d_current,tau_current),sqrt(2*prop_k_b),log=T) -
          #     dnorm(b_proposal[i],b_current[i] + prop_k_b*grad_log_target_2.6_b(data,i,alphaU_current,alphaR_current,b_current,d_current,tau_current),sqrt(2*prop_k_b),log=T)
          # },numeric(1))) +
          # random walk proposal dist for d
          dlnorm(d_current,log(d_proposal),prop_sd_d,log=T) - dlnorm(d_proposal,log(d_current),prop_sd_d,log=T)
          # MALA proposal dist for d
          # dlnorm(d_current,log(d_proposal) + prop_k_logd*grad_log_target_2.6_logd(data,alphaU_current,alphaR_current,b_current,d_proposal,mu_d,sd_d),sqrt(2*prop_k_logd),log=T) -
          # dlnorm(d_proposal,log(d_current) + prop_k_logd*grad_log_target_2.6_logd(data,alphaU_current,alphaR_current,b_current,d_current,mu_d,sd_d),sqrt(2*prop_k_logd),log=T)

      })
      
      if(a_prob>runif(1)){
        d_current <- d_proposal
        b_current <- b_proposal
        alphaU_current <- alphaU_proposal
        alphaR_current <- alphaR_proposal
        if(i>burnin){acc <- acc + 1}
      }
      
      # record current state
      tau_postsamp_t[i] <- tau_current
      d_postsamp_t[i] <- d_current
      alphaU_postsamp_t[i] <- alphaU_current
      alphaR_postsamp_t[i] <- alphaR_current
      b_postsamp_t[i,] <- b_current
      
    } 
    
    #get rid of burn-in and add to other chains
    b_postsamp <- rbind(b_postsamp,b_postsamp_t[(burnin+1):n_iter,])
    alphaU_postsamp <- c(alphaU_postsamp,alphaU_postsamp_t[(burnin+1):n_iter])
    alphaR_postsamp <- c(alphaR_postsamp,alphaR_postsamp_t[(burnin+1):n_iter])
    tau_postsamp <- c(tau_postsamp,tau_postsamp_t[(burnin+1):n_iter])
    d_postsamp <- c(d_postsamp,d_postsamp_t[(burnin+1):n_iter])
  }
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp, alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  acc_rate <- acc/((n_iter-burnin)*chains)
  
  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              tau=tau_postsamp, d=d_postsamp, 
              acc_rate=acc_rate))
}
#######################################################################
#######################################################################
# Alg 2.5 and 2.6 simulations ----

## settings
{
  n_areas <- 50 #numbers of areas
  n_clusters_urban <- rnegbin(n_areas,15,5) #number of urban clusters in each area (more accurate is around 250)
  n_clusters_rural <- rnegbin(n_areas,15,5) #number of rural clusters in each area (more accurate is around 250)
  n_births_urban <- 100 #average number of births per urban cluster
  n_births_rural <- 75 #average number of births per rural cluster
}

## generate data 
{
  # intercept (fixed effect) 
  alphaU <- rnorm(1,-3.5,0.25)
  alphaR <- rnorm(1,-3,0.25)
  
  # draw hyperparameter
  tau <- rgamma(1,10,.25)
  
  # draw overdispersion parameter
  d <- rhalfnorm(1,5)
  
  # draw random effects
  b <- as.vector(Rfast::rmvnorm(1,rep(0,n_areas),diag(1/tau,n_areas)))
  
  all_dat <- data.frame(
    # (true) number of births in each cluster
    N = c(rpois(sum(n_clusters_urban),n_births_urban),rpois(sum(n_clusters_rural),n_births_rural)),
    # urban or rural strata (U=1 if urban, otherwise 0)
    U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
    # admin area of each cluster
    A = c(unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_rural[i])}))))
  all_dat$cluster <- 1:nrow(all_dat)
  all_dat$Y <- sapply(1:nrow(all_dat),function(i){rnbinom(1,
                                                          size = 1/d*all_dat$N[i]*exp(alphaU*all_dat$U[i] + alphaR*(1-all_dat$U[i]) + b[all_dat$A[i]] + rnorm(1,0,0.05)),
                                                          prob=1/(1+d))})
  
}

## run MCMC
postsamp_2.5 <- postsamp_Alg2.5_MCMC(tau_init = 40, d_init = 0.1, alphaU_init = -6, alphaR_init = -5,
                                     prop_sd_d = 0.15, prop_sd_b = 0.05, prop_sd_alpha = 0.05,
                                     data = all_dat, Q_inv = Q_scaled_inv, n_iter = 50000, chains = 2)

pdf("Handcoded MCMC Simulations/Alg 2.5 231005, block updates except tau, 50k iterations posterior distributions.pdf") 
{
  hist(log(postsamp_2.5$tau))
  abline(v=log(tau),col='red')
  hist(postsamp_2.5$d)
  abline(v=d,col='red')
  
  hist(postsamp_2.5$alphaU)
  abline(v=alphaU,col='red')
  hist(postsamp_2.5$alphaR)
  abline(v=alphaR,col='red')
  
  results2.5_b <- gather(data.frame(postsamp_2.5$b),'area','b')
  results2.5_b$area <- as.numeric(str_remove(results2.5_b$area,'X'))
  true_b <- data.frame(area=1:n_areas,true_b=as.vector(b))
  true_b <- true_b[order(true_b$true_b),]
  true_b$b_rank <- 1:nrow(true_b)
  results2.5_b <- left_join(results2.5_b,true_b,by='area')
  g <- results2.5_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~b_rank) + ggtitle('XX births per area') +
    geom_hline(yintercept=0,col='grey50')
  print(g)
  
  plot(log(postsamp_2.5$tau[1:50000%%50==0]),type='l')
  abline(h=log(tau),col='red')
  plot(postsamp_2.5$d[1:50000%%50==0],type='l')
  abline(h=d,col='red')
  plot(postsamp_2.5$alphaU[1:50000%%50==0],type='l')
  abline(h=alphaU,col='red')
  plot(postsamp_2.5$alphaR[1:50000%%50==0],type='l')
  abline(h=alphaR,col='red')
  plot(postsamp_2.5$b[1:50000%%50==0,7],type='l')
  abline(h=b[7],col='red')
  plot(postsamp_2.5$b[1:50000%%50==0,15],type='l')
  abline(h=b[15],col='red')
  plot(postsamp_2.5$b[1:50000%%50==0,28],type='l')
  abline(h=b[28],col='red')
  plot(postsamp_2.5$b[1:50000%%50==0,42],type='l')
  abline(h=b[42],col='red')
} 
dev.off()  

postsamp_2.6 <- postsamp_Alg2.6_MCMC(tau_init = 25, d_init = 0.1, alphaU_init = -3, alphaR_init = -3,
                                     prop_sd_d = 0.05, prop_sd_b = 0.01, prop_sd_alpha = 0.05,
                                     #prop_k_logd = 0.1, prop_k_b = 2e-3, prop_k_alpha = 1e-3,
                                     data = all_dat, n_iter = 1000, chains = 1, burnin=0)

pdf("Handcoded MCMC Simulations/Alg 2.6 231010, MALA block updates except tau, 50k iterations posterior distributions.pdf") 
{
  hist(log(postsamp_2.6$tau))
  abline(v=log(tau),col='red')
  hist(postsamp_2.6$d)
  abline(v=d,col='red')
  
  hist(postsamp_2.6$alphaU)
  abline(v=alphaU,col='red')
  hist(postsamp_2.6$alphaR)
  abline(v=alphaR,col='red')
  
  results2.6_b <- gather(data.frame(postsamp_2.6$b),'area','b')
  results2.6_b$area <- as.numeric(str_remove(results2.6_b$area,'X'))
  true_b <- data.frame(area=1:n_areas,true_b=as.vector(b))
  true_b <- true_b[order(true_b$true_b),]
  true_b$b_rank <- 1:nrow(true_b)
  results2.6_b <- left_join(results2.6_b,true_b,by='area')
  g <- results2.6_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~b_rank) + ggtitle('Approx 2500 births per area') +
    geom_hline(yintercept=0,col='grey50')
  print(g)
  
  plot(log(postsamp_2.6$tau[1:50000%%50==0]),type='l')
  abline(h=log(tau),col='red')
  plot(postsamp_2.6$d[1:50000%%50==0],type='l')
  abline(h=d,col='red')
  plot(postsamp_2.6$alphaU[1:50000%%50==0],type='l')
  abline(h=alphaU,col='red')
  plot(postsamp_2.6$alphaR[1:50000%%50==0],type='l')
  abline(h=alphaR,col='red')
  plot(postsamp_2.6$b[1:50000%%50==0,7],type='l')
  abline(h=b[7],col='red')
  plot(postsamp_2.6$b[1:50000%%50==0,15],type='l')
  abline(h=b[15],col='red')
  plot(postsamp_2.6$b[1:50000%%50==0,28],type='l')
  abline(h=b[28],col='red')
  plot(postsamp_2.6$b[1:50000%%50==0,42],type='l')
  abline(h=b[42],col='red')
} 
dev.off()  

# get distribution of posterior mean over several iterations
{
  nsims <- 10
  results2.6_tau <- rep(NA,nsims)
  results2.6_d <- rep(NA,nsims)
  results2.6_alphaU <- rep(NA,nsims)
  results2.6_alphaR <- rep(NA,nsims)
  results2.6_b <- matrix(NA,nsims,n_areas)
  
  for(k in 1:nsims){
    #status check
    if(k %% (nsims/10)==0){cat('Starting simulation', k, '\n')}
    
    #generate data
    all_dat <- data.frame(
      # (true) number of births in each cluster
      N = c(rpois(sum(n_clusters_urban),n_births_urban),rpois(sum(n_clusters_rural),n_births_rural)),
      # urban or rural strata (U=1 if urban, otherwise 0)
      U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
      # admin area of each cluster
      A = c(unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_rural[i])}))))
    all_dat$cluster <- 1:nrow(all_dat)
    all_dat$Y <- sapply(1:nrow(all_dat),function(i){rnbinom(1,
                                                            size = 1/d*all_dat$N[i]*exp(alphaU*all_dat$U[i] + alphaR*(1-all_dat$U[i]) + b[all_dat$A[i]] + rnorm(1,0,0.05)),
                                                            prob=1/(1+d))})
    
    #run MCMC
    postsamp_2.6 <- postsamp_Alg2.6_MCMC(tau_init = 25, d_init = 0.1, alphaU_init = -6, alphaR_init = -5,
                                         prop_k_logd = 0.1, prop_k_b = 1e-2, prop_k_alpha = 2e-3,
                                         data = all_dat, n_iter = 50000, chains = 2)
    
    #record result
    results2.6_tau[k] <- mean(postsamp_2.6$tau)
    results2.6_d[k] <- mean(postsamp_2.6$d)
    results2.6_alphaU[k] <- mean(postsamp_2.6$alphaU)
    results2.6_alphaR[k] <- mean(postsamp_2.6$alphaR)
    results2.6_b[k,] <- colMeans(postsamp_2.6$b)
    
  }
  
  results_2.6 <- list(tau_res = results2.6_tau, d_res = results2.6_d,
                      alphaU_res = results2.6_alphaU, alphaR_res = results2.6_alphaR, b_res = results2.6_b, 
                      tau_truth = tau, d_truth=d, alphaU_truth=alphaU, alphaR_truth=alphaR, b_truth = b)
  
  save(results_2.6,file = 'Handcoded MCMC Simulations/Alg 2.6 231010 50k iterations for 10 datasets.rda')
}

pdf("Handcoded MCMC Simulations/Alg 2.6 231010 same b, 50k iterations.pdf")
{
  hist(log(results2.6_tau), main="Posterior means from 10 data realizations")
  abline(v=log(tau),col='red')
  hist(log(results2.6_d), main="Posterior means from 10 data realizations")
  abline(v=log(d),col='red')
  hist(results2.6_alphaU, main="Posterior means from 10 data realizations")
  abline(v=alphaU,col='red')
  hist(results2.6_alphaR, main="Posterior means from 10 data realizations")
  abline(v=alphaR,col='red')
  
  results_b <- gather(data.frame(results_2.6$b_res),'area','b')
  results_b$area <- as.numeric(str_remove(results_b$area,'X'))
  true_b <- data.frame(area=1:n_areas,true_b=results_2.6$b_truth)
  results_b <- left_join(results_b,true_b,by='area')
  g <- results_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~area) + ggtitle('Posterior means from 10 data realizations')
  print(g)
  
}
dev.off()


#######################################################################
#######################################################################
# Alg 3.1: MCMC to get P(Y,r,a,tau,phi|Z) (BYM2 spatial effect, SRS random sampling) ------
postsamp_Yartauphi_Z_MCMC <- function(tau_init,phi_init,alpha_init, #initial values
                                      prop_sd_phi, prop_sd_tau, prop_Y_J,
                                      data,Q_inv,
                                      n_iter,burnin=n_iter/2,chains=2){
  if(length(unique(data$A))!=nrow(data)){
    stop('Data should be at area level')
  }
  n_areas <- length(unique(data$A))
  #hyperpriors
  a_tau <- b_tau <- 0.01
  a_phi <- 2
  b_phi <- 3
  mu_a <- -6
  sd_a <- 3
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  phi_postsamp_t <- c(phi_init, rep(NA,n_iter-1))
  alpha_postsamp_t <- c(alpha_init, rep(NA,n_iter-1))
  acc_bym2 <- acc_Y <- 0
  
  eta_init <- Rfast::rmvnorm(1,rep(alpha_init,n_areas),1/tau_init*(W(phi_init,Q_inv))) 
  eta_postsamp_t <- rbind(eta_init, matrix(NA,n_iter-1,n_areas))
  
  Y_init <- sapply(1:n_areas,function(area){rpois(1,data$N[area]*exp(eta_init[area]))})
  Y_postsamp_t <- rbind(Y_init, matrix(NA,n_iter-1,n_areas))
  
  tau_postsamp <- phi_postsamp <- alpha_postsamp <- eta_postsamp <- Y_postsamp <- NULL
  for(chain in 1:chains){
    for(i in 2:n_iter){
      tau_current <- tau_postsamp_t[i-1]
      phi_current <- phi_postsamp_t[i-1]
      alpha_current <- alpha_postsamp_t[i-1]
      eta_current <- eta_postsamp_t[(i-1),]
      Y_current <- Y_postsamp_t[(i-1),]
      
      # Draw proposal for phi
      phi_proposal <- expit(rnorm(1,logitlink(phi_current),prop_sd_phi)) 
      
      # Draw proposal for tau
      tau_proposal <- rlnorm(1,log(tau_current), prop_sd_tau)
      
      # Draw proposal for eta for Taylor approximation of full conditional
      eta_proposal <- as.vector(Rfast::rmvnorm(1,prop_mean_bym2_eta(eta_current,tau_proposal,phi_proposal,alpha_current,data,Q_inv,Y=Y_current),
                                               prop_sigma_bym2_eta(eta_current,tau_proposal,phi_proposal,data,Q_inv)))
      
      aprob_bym2 <- exp( #Y|eta
        sum(sapply(1:n_areas,function(area){
          dpois(Y_current[area],data$N[area]*exp(eta_proposal[area]),log=T) - #Y|eta
            dpois(Y_current[area],data$N[area]*exp(eta_current[area]),log=T)})) +
          #eta-alpha|phi,tau
          Rfast::dmvnorm(eta_proposal,rep(alpha_current,n_areas),sigma= 1/tau_proposal*W(phi_proposal,Q_inv),logged = T) - 
          Rfast::dmvnorm(eta_current,rep(alpha_current,n_areas),sigma= 1/tau_current*W(phi_current,Q_inv),logged = T) +
          #prior on tau 
          dgamma(tau_proposal,a_tau,b_tau,log = T) - dgamma(tau_current,a_tau,b_tau,log = T) + 
          #prior on phi
          dbeta(phi_proposal,a_phi,b_phi,log=T) - dbeta(phi_current,a_phi,b_phi,log=T) +
          #Jacobian to go from logitphi to phi
          (-log(phi_current*(1-phi_current)) + log(phi_proposal*(1-phi_proposal))) + 
          #proposal distribution for logitphi is symmetric
          #proposal dist for tau
          dlnorm(tau_current,log(tau_proposal),prop_sd_tau,log=T) - 
          dlnorm(tau_proposal,log(tau_current),prop_sd_tau,log=T) +
          #proposal dist for alpha is symmetric
          #proposal dist for eta
          Rfast::dmvnorm(eta_current,prop_mean_bym2_eta(eta_proposal,tau_current,phi_current,alpha_current,data,Q_inv,Y=Y_current),
                         prop_sigma_bym2_eta(eta_proposal,tau_current,phi_current,data,Q_inv),logged=T) - 
          Rfast::dmvnorm(eta_proposal,prop_mean_bym2_eta(eta_current,tau_proposal,phi_proposal,alpha_current,data,Q_inv,Y=Y_current),
                         prop_sigma_bym2_eta(eta_current,tau_proposal,phi_proposal,data,Q_inv),logged=T)
      )
      if(runif(1) < aprob_bym2){
        phi_current <- phi_proposal
        tau_current <- tau_proposal
        eta_current <- eta_proposal
        if(i>burnin){
          acc_bym2 <- acc_bym2+1 
        }
      }
      
      #draw alpha from full conditional
      alpha_current <- rnorm(1,(mu_a+sd_a^2*tau_current*sum(t(eta_current)%*%solve(W(phi_current,Q_inv))))/(1+sd_a^2*tau_current*sum(solve(W(phi_current,Q_inv)))),
                             1/sqrt((1/sd_a^2) + tau_current*sum(solve(W(phi_current,Q_inv)))))
      
      
      ## Draw proposal for Y from discerete uniform centered at current state
      Y_proposal <- Y_current + sample((-prop_Y_J):prop_Y_J,n_areas,replace=T)
      
      if(sum(Y_proposal<0)>0){
        aprob_Y  <- 0
      }else{
        aprob_Y <- exp(sum(sapply(1:n_areas,function(area){
          dhyper(data$Z[area],Y_proposal[area],data$N[area] - Y_proposal[area],data$n[area],log = T) - #Z|Y
            dhyper(data$Z[area],Y_current[area],data$N[area] - Y_current[area],data$n[area],log = T)})) + 
            sum(sapply(1:n_areas,function(area){
              dpois(Y_proposal[area],data$N[area]*exp(eta_current[area]),log=T) - #Y|eta
                dpois(Y_current[area],data$N[area]*exp(eta_current[area]),log=T)}))
          #proposal distribution is symmetric
        )}
      
      if(runif(1) < aprob_Y){
        Y_current <- Y_proposal
        if(i>burnin){
          acc_Y <- acc_Y+1 
        }}
      
      # record state
      tau_postsamp_t[i] <- tau_current
      phi_postsamp_t[i] <- phi_current
      alpha_postsamp_t[i] <- alpha_current
      eta_postsamp_t[i,] <- eta_current
      Y_postsamp_t[i,] <- Y_current
      
    } 
    #get rid of burn-in and add to other chains
    Y_postsamp <- rbind(Y_postsamp,Y_postsamp_t[(burnin+1):n_iter,])
    eta_postsamp <- rbind(eta_postsamp,eta_postsamp_t[(burnin+1):n_iter,])
    alpha_postsamp <- c(alpha_postsamp,alpha_postsamp_t[(burnin+1):n_iter])
    tau_postsamp <- c(tau_postsamp,tau_postsamp_t[(burnin+1):n_iter])
    phi_postsamp <- c(phi_postsamp,phi_postsamp_t[(burnin+1):n_iter])
  }
  
  b_postsamp <- eta_postsamp-alpha_postsamp
  r_postsamp <- exp(eta_postsamp)
  acc_rates <- matrix(c(acc_bym2/((n_iter-burnin)*chains),
                        acc_Y/((n_iter-burnin)*chains)), nrow=1)
  colnames(acc_rates) <- c('bym2','Y') 
  
  return(list(alpha = alpha_postsamp,b=b_postsamp, eta=eta_postsamp,r=r_postsamp, Y=Y_postsamp, tau=tau_postsamp,phi=phi_postsamp,acc_rates=acc_rates))
}

#######################################################################
#######################################################################
# Alg 3.2: MCMC to get P(Y,r,a,tau,phi|Z,Y+) (BYM2 spatial effect, SRS random sampling) ----

tau_init = 25
phi_init = 0.25
alpha_init = -6
prop_sd_tau = 0.25
prop_sd_phi = 0.25
prop_Y_f_max = 5
data = obs_dat
Q_inv = Q_scaled_inv
n_iter = 2000
chains=1


postsamp_Yartauphi_ZYplus_MCMC <- function(tau_init,phi_init,alpha_init, #initial values
                                           prop_sd_phi, prop_sd_tau, prop_Y_f_max,
                                           data,Q_inv,
                                           n_iter,burnin=n_iter/2,chains=2){
  if(length(unique(data$A))!=nrow(data)){
    stop('Data should be at area level')
  }
  n_areas <- length(unique(data$A))
  Yplus <- sum(data$Y)
  #hyperpriors
  a_tau <- b_tau <- 0.01
  a_phi <- 2
  b_phi <- 3
  mu_a <- -6
  sd_a <- 3
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  phi_postsamp_t <- c(phi_init, rep(NA,n_iter-1))
  alpha_postsamp_t <- c(alpha_init, rep(NA,n_iter-1))
  acc_bym2 <- acc_Y <- 0
  
  eta_init <- Rfast::rmvnorm(1,rep(alpha_init,n_areas),1/tau_init*(W(phi_init,Q_inv))) 
  eta_postsamp_t <- rbind(eta_init, matrix(NA,n_iter-1,n_areas))
  
  Y_init <- t(rmultinom(1,Yplus,theta(eta_init,data)))
  Y_postsamp_t <- rbind(Y_init, matrix(NA,n_iter-1,n_areas))
  
  tau_postsamp <- phi_postsamp <- alpha_postsamp <- eta_postsamp <- Y_postsamp <- NULL
  for(chain in 1:chains){
    for(i in 2:n_iter){
      tau_current <- tau_postsamp_t[i-1]
      phi_current <- phi_postsamp_t[i-1]
      alpha_current <- alpha_postsamp_t[i-1]
      eta_current <- eta_postsamp_t[(i-1),]
      Y_current <- Y_postsamp_t[(i-1),]
      
      # Draw proposal for phi
      phi_proposal <- expit(rnorm(1,logitlink(phi_current),prop_sd_phi)) 
      
      # Draw proposal for tau
      tau_proposal <- rlnorm(1,log(tau_current), prop_sd_tau)
      
      # Draw proposal for eta for Taylor approximation of full conditional
      eta_proposal <- as.vector(Rfast::rmvnorm(1,prop_mean_bym2_eta(eta_current,tau_proposal,phi_proposal,alpha_current,data,Q_inv,Y=Y_current),
                                               prop_sigma_bym2_eta(eta_current,tau_proposal,phi_proposal,data,Q_inv)))
      
      log_aprob_bym2 <- #Y|eta
        sum(sapply(1:n_areas,function(area){
          dpois(Y_current[area],data$N[area]*exp(eta_proposal[area]),log=T) - #Y|eta
            dpois(Y_current[area],data$N[area]*exp(eta_current[area]),log=T)})) +
        #eta-alpha|phi,tau
        Rfast::dmvnorm(eta_proposal,rep(alpha_current,n_areas),sigma= 1/tau_proposal*W(phi_proposal,Q_inv),logged = T) - 
        Rfast::dmvnorm(eta_current,rep(alpha_current,n_areas),sigma= 1/tau_current*W(phi_current,Q_inv),logged = T) +
        #prior on tau 
        dgamma(tau_proposal,a_tau,b_tau,log = T) - dgamma(tau_current,a_tau,b_tau,log = T) + 
        #prior on phi
        dbeta(phi_proposal,a_phi,b_phi,log=T) - dbeta(phi_current,a_phi,b_phi,log=T) +
        #Jacobian to go from logitphi to phi
        (-log(phi_current*(1-phi_current)) + log(phi_proposal*(1-phi_proposal))) + 
        #proposal distribution for logitphi is symmetric
        #proposal dist for tau
        dlnorm(tau_current,log(tau_proposal),prop_sd_tau,log=T) - 
        dlnorm(tau_proposal,log(tau_current),prop_sd_tau,log=T) +
        #proposal dist for alpha is symmetric
        #proposal dist for eta
        Rfast::dmvnorm(eta_current,prop_mean_bym2_eta(eta_proposal,tau_current,phi_current,alpha_current,data,Q_inv,Y=Y_current),
                       prop_sigma_bym2_eta(eta_proposal,tau_current,phi_current,data,Q_inv),logged=T) - 
        Rfast::dmvnorm(eta_proposal,prop_mean_bym2_eta(eta_current,tau_proposal,phi_proposal,alpha_current,data,Q_inv,Y=Y_current),
                       prop_sigma_bym2_eta(eta_current,tau_proposal,phi_proposal,data,Q_inv),logged=T)
      
      if(runif(1) < exp(log_aprob_bym2)){
        phi_current <- phi_proposal
        tau_current <- tau_proposal
        eta_current <- eta_proposal
        if(i>burnin){
          acc_bym2 <- acc_bym2+1 
        }
      }
      
      #draw alpha from full conditional
      alpha_current <- rnorm(1,(mu_a+sd_a^2*tau_current*sum(t(eta_current)%*%solve(W(phi_current,Q_inv))))/(1+sd_a^2*tau_current*sum(solve(W(phi_current,Q_inv)))),
                             1/sqrt((1/sd_a^2) + tau_current*sum(solve(W(phi_current,Q_inv)))))
      
      ## Draw proposal for Y while maintaining Yplus
      inc_size <- sample(0:prop_Y_f_max,floor(n_areas/2),replace = T) #choose change in increment for each pair of Ys
      Y_change_id <- sample(1:n_areas,n_areas,replace = F) #choose areas which will be increased and which will be decreased
      if(n_areas%%2==0){
        Y_changes <- data.frame(area=Y_change_id,inc_size=c(inc_size,-1*inc_size))
      }else{
        Y_changes <- data.frame(area=Y_change_id,inc_size=c(inc_size,-1*inc_size,0))
      }
      
      Y_changes <- Y_changes[order(Y_changes$area),]
      Y_proposal <- Y_current + Y_changes$inc_size
      
      if(sum(Y_proposal<0)>0){
        aprob_Y  <- 0
      }else{
        aprob_Y <-  exp(sum(sapply(1:n_areas,function(area){
          dhyper(data$Z[area],Y_proposal[area],data$N[area] - Y_proposal[area],data$n[area],log = T) - #Z|Y
            dhyper(data$Z[area],Y_current[area],data$N[area] - Y_current[area],data$n[area],log = T)})) + 
            dmultinom(Y_proposal,prob=theta(eta_current,data),log=T) - #Y|Y+,eta
            dmultinom(Y_current,prob=theta(eta_current,data),log=T)
          #proposal distribution is symmetric
        )
      }
      
      if(runif(1) < aprob_Y){
        Y_current <- Y_proposal
        if(i>burnin){
          acc_Y <- acc_Y+1 
        }
      }
      
      # record state
      tau_postsamp_t[i] <- tau_current
      phi_postsamp_t[i] <- phi_current
      alpha_postsamp_t[i] <- alpha_current
      eta_postsamp_t[i,] <- eta_current
      Y_postsamp_t[i,] <- Y_current
      
    } 
    #get rid of burn-in and add to other chains
    Y_postsamp <- rbind(Y_postsamp,Y_postsamp_t[(burnin+1):n_iter,])
    eta_postsamp <- rbind(eta_postsamp,eta_postsamp_t[(burnin+1):n_iter,])
    alpha_postsamp <- c(alpha_postsamp,alpha_postsamp_t[(burnin+1):n_iter])
    tau_postsamp <- c(tau_postsamp,tau_postsamp_t[(burnin+1):n_iter])
    phi_postsamp <- c(phi_postsamp,phi_postsamp_t[(burnin+1):n_iter])
    
  }
  
  b_postsamp <- eta_postsamp-alpha_postsamp
  r_postsamp <- exp(eta_postsamp)
  acc_rates <- matrix(c(acc_bym2/((n_iter-burnin)*chains),
                        acc_Y/((n_iter-burnin)*chains)), nrow=1)
  colnames(acc_rates) <- c('bym2','Y') 
  
  return(list(alpha = alpha_postsamp,b=b_postsamp, eta=eta_postsamp,r=r_postsamp, 
              Y=Y_postsamp, tau=tau_postsamp,phi=phi_postsamp,acc_rates=acc_rates))
}
#######################################################################
#######################################################################
# Alg 3.3: MCMC to estimate P(Y,Y_+,r,a,tau,phi|Z,Y+_hat) (BYM2 spatial effect, SRS random sampling) ------

postsamp_Alg3.3_MCMC <- function(tau_init,phi_init,alpha_init, #initial values
                                 prop_sd_phi, prop_sd_tau, prop_Y_f_max,
                                 data,Q_inv,
                                 n_iter,burnin=n_iter/2,chains=2){
  if(length(unique(data$A))!=nrow(data)){
    stop('Data should be at area level')
  }
  n_areas <- length(unique(data$A))
  ## calculate direct estimate of Y_+
  Yplus_hat <- sum(data$Z*data$N/data$n)
  Yplus_var_hat <- sum(data$N*(data$N-data$n)/(data$n^2*(data$n-1))*data$Z*(data$n-data$Z))
  
  #hyperpriors
  a_tau <- b_tau <- 0.01
  a_phi <- 2
  b_phi <- 3
  mu_a <- -6
  sd_a <- 3
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  phi_postsamp_t <- c(phi_init, rep(NA,n_iter-1))
  alpha_postsamp_t <- c(alpha_init, rep(NA,n_iter-1))
  acc_bym2 <- acc_Y <- 0
  
  eta_init <- Rfast::rmvnorm(1,rep(alpha_init,n_areas),1/tau_init*(W(phi_init,Q_inv))) 
  eta_postsamp_t <- rbind(eta_init, matrix(NA,n_iter-1,n_areas))
  
  Yplus_init <- round(rnorm(1,Yplus_hat,sqrt(Yplus_var_hat)))
  Yplus_postsamp_t <- c(Yplus_init, rep(NA,n_iter-1))
  
  Y_init <- t(rmultinom(1,Yplus_init,theta(eta_init,data)))
  Y_postsamp_t <- rbind(Y_init, matrix(NA,n_iter-1,n_areas))
  
  tau_postsamp <- phi_postsamp <- alpha_postsamp <- eta_postsamp <- Y_postsamp <- Yplus_postsamp <- NULL
  for(chain in 1:chains){
    for(i in 2:n_iter){
      tau_current <- tau_postsamp_t[i-1]
      phi_current <- phi_postsamp_t[i-1]
      alpha_current <- alpha_postsamp_t[i-1]
      eta_current <- eta_postsamp_t[(i-1),]
      Y_current <- Y_postsamp_t[(i-1),]
      Yplus_current <- Yplus_postsamp_t[i-1]
      
      # Draw proposal for phi
      phi_proposal <- expit(rnorm(1,logitlink(phi_current),prop_sd_phi)) 
      
      # Draw proposal for tau
      tau_proposal <- rlnorm(1,log(tau_current), prop_sd_tau)
      
      # Draw proposal for eta for Taylor approximation of full conditional
      eta_proposal <- as.vector(Rfast::rmvnorm(1,prop_mean_bym2_eta(eta_current,tau_proposal,phi_proposal,alpha_current,data,Q_inv,Y=Y_current),
                                               prop_sigma_bym2_eta(eta_current,tau_proposal,phi_proposal,data,Q_inv)))
      
      log_aprob_bym2 <- #Y|eta
        sum(sapply(1:n_areas,function(area){
          dpois(Y_current[area],data$N[area]*exp(eta_proposal[area]),log=T) - #Y|eta
            dpois(Y_current[area],data$N[area]*exp(eta_current[area]),log=T)})) +
        #eta-alpha|phi,tau
        Rfast::dmvnorm(eta_proposal,rep(alpha_current,n_areas),sigma= 1/tau_proposal*W(phi_proposal,Q_inv),logged = T) - 
        Rfast::dmvnorm(eta_current,rep(alpha_current,n_areas),sigma= 1/tau_current*W(phi_current,Q_inv),logged = T) +
        #prior on tau 
        dgamma(tau_proposal,a_tau,b_tau,log = T) - dgamma(tau_current,a_tau,b_tau,log = T) + 
        #prior on phi
        dbeta(phi_proposal,a_phi,b_phi,log=T) - dbeta(phi_current,a_phi,b_phi,log=T) +
        #Jacobian to go from logitphi to phi
        (-log(phi_current*(1-phi_current)) + log(phi_proposal*(1-phi_proposal))) + 
        #proposal distribution for logitphi is symmetric
        #proposal dist for tau
        dlnorm(tau_current,log(tau_proposal),prop_sd_tau,log=T) - 
        dlnorm(tau_proposal,log(tau_current),prop_sd_tau,log=T) +
        #proposal dist for alpha is symmetric
        #proposal dist for eta
        Rfast::dmvnorm(eta_current,prop_mean_bym2_eta(eta_proposal,tau_current,phi_current,alpha_current,data,Q_inv,Y=Y_current),
                       prop_sigma_bym2_eta(eta_proposal,tau_current,phi_current,data,Q_inv),logged=T) - 
        Rfast::dmvnorm(eta_proposal,prop_mean_bym2_eta(eta_current,tau_proposal,phi_proposal,alpha_current,data,Q_inv,Y=Y_current),
                       prop_sigma_bym2_eta(eta_current,tau_proposal,phi_proposal,data,Q_inv),logged=T)
      
      if(runif(1) < exp(log_aprob_bym2)){
        phi_current <- phi_proposal
        tau_current <- tau_proposal
        eta_current <- eta_proposal
        if(i>burnin){
          acc_bym2 <- acc_bym2+1 
        }
      }
      
      #draw alpha from full conditional
      alpha_current <- rnorm(1,(mu_a+sd_a^2*tau_current*sum(t(eta_current)%*%solve(W(phi_current,Q_inv))))/(1+sd_a^2*tau_current*sum(solve(W(phi_current,Q_inv)))),
                             1/sqrt((1/sd_a^2) + tau_current*sum(solve(W(phi_current,Q_inv)))))
      
      #draw Yplus from full conditional
      Yplus_current <- round(rnorm(1,Yplus_hat,sqrt(Yplus_var_hat)))
      
      ## Randomly adjust Ys to get new Yplus
      Y_adj_inc <- as.vector(rmultinom(1,abs(sum(Y_current)-Yplus_current),rep(1/n_areas,n_areas)))
      Y_proposal_t <- Y_current + sign(Yplus_current-sum(Y_current))*Y_adj_inc
      
      ## Draw proposal for Y while maintaining new Yplus
      inc_size <- sample(0:prop_Y_f_max,floor(n_areas/2),replace = T) #choose change in increment for each pair of Ys
      Y_change_id <- sample(1:n_areas,n_areas,replace = F) #choose areas which will be increased and which will be decreased
      if(n_areas%%2==0){
        Y_changes <- data.frame(area=Y_change_id,inc_size=c(inc_size,-1*inc_size))
      }else{
        Y_changes <- data.frame(area=Y_change_id,inc_size=c(inc_size,-1*inc_size,0))
      }
      
      Y_changes <- Y_changes[order(Y_changes$area),]
      Y_proposal <- Y_proposal_t + Y_changes$inc_size
      
      if(sum(Y_proposal<0)>0){
        aprob_Y  <- 0
      }else{
        aprob_Y <-  exp(sum(sapply(1:n_areas,function(area){
          dhyper(data$Z[area],Y_proposal[area],data$N[area] - Y_proposal[area],data$n[area],log = T) - #Z|Y
            dhyper(data$Z[area],Y_current[area],data$N[area] - Y_current[area],data$n[area],log = T)})) + 
            dmultinom(Y_proposal,prob=theta(eta_current,data),log=T) - #Y|Y+,eta
            dmultinom(Y_current,prob=theta(eta_current,data),log=T)
          #proposal distribution is symmetric
        )
      }
      
      if(runif(1) < aprob_Y){
        Y_current <- Y_proposal
        if(i>burnin){
          acc_Y <- acc_Y+1 
        }
      }
      
      # record state
      tau_postsamp_t[i] <- tau_current
      phi_postsamp_t[i] <- phi_current
      alpha_postsamp_t[i] <- alpha_current
      eta_postsamp_t[i,] <- eta_current
      Y_postsamp_t[i,] <- Y_current
      Yplus_postsamp_t[i] <- Yplus_current
      
    } 
    #get rid of burn-in and add to other chains
    Y_postsamp <- rbind(Y_postsamp,Y_postsamp_t[(burnin+1):n_iter,])
    Yplus_postsamp <- c(Yplus_postsamp,Yplus_postsamp_t[(burnin+1):n_iter])
    eta_postsamp <- rbind(eta_postsamp,eta_postsamp_t[(burnin+1):n_iter,])
    alpha_postsamp <- c(alpha_postsamp,alpha_postsamp_t[(burnin+1):n_iter])
    tau_postsamp <- c(tau_postsamp,tau_postsamp_t[(burnin+1):n_iter])
    phi_postsamp <- c(phi_postsamp,phi_postsamp_t[(burnin+1):n_iter])
    
  }
  
  b_postsamp <- eta_postsamp-alpha_postsamp
  r_postsamp <- exp(eta_postsamp)
  acc_rates <- matrix(c(acc_bym2/((n_iter-burnin)*chains),
                        acc_Y/((n_iter-burnin)*chains)), nrow=1)
  colnames(acc_rates) <- c('bym2','Y') 
  
  return(list(alpha = alpha_postsamp,b=b_postsamp, eta=eta_postsamp,r=r_postsamp, 
              Y=Y_postsamp, Yplus=Yplus_postsamp, tau=tau_postsamp,phi=phi_postsamp,acc_rates=acc_rates))
}


#######################################################################
#######################################################################
# BYM2 random effect, SRS simulations (Alg 3.1, 3.2, 3.3) ------------

## define ground truth
{
  n_areas <- 50 #numbers of areas
  n_clusters <- rep(500,n_areas) #number of clusters in each area
  n_births <- 100 #average number of births per cluster
  n_births_samp <- rep(5000,n_areas) #number of births sampled from each area
  
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
  Q_scaled <- inla.scale.model(Q, constr=list(A=t(eigen(Q)$vectors[,eigen(Q)$values<1e-10]), e=rep(0,sum(eigen(Q)$values<1e-10))))
  Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))
  #get approximation of eigenvalues for Q_scaled_inv
  gamma <- eigen(Q_scaled)$values
  gamma_til <- 1/gamma
  gamma_til[gamma_til>1e+10] <- 0
  
  # intercept (fixed effect) 
  alpha <- rnorm(1,-6,0.25)
  
  # draw hyperparameters
  tau <- rgamma(1,10,.25)
  phi <- rbeta(1,2,3)
}

## generate data 
{
  Var_b <- (1/tau)*(diag((1-phi),n_areas) + phi*Q_scaled_inv)
  b <- Rfast::rmvnorm(1,rep(0,n_areas),Var_b)
  all_dat <- simdat_bym2sp(n_areas=n_areas,n_clusters=n_clusters,n_births=n_births,a=alpha,Var_b=Var_b,b=b)$dat
  obs_dat <- get_srs(all_dat,n_births_samp,'A')
}

## Run MCMC for 1 dataset -- doesn't mix well but curious to see if this will change when Y is sampled differently in 3.2
postsamp_3.1 <- postsamp_Yartauphi_Z_MCMC(tau_init = 25,phi_init = 0.25,alpha_init = -6,
                                          prop_sd_tau = 0.25, prop_sd_phi = 0.25,
                                          prop_Y_J = 2, data= obs_dat, Q_inv = Q_scaled_inv, n_iter = 5000)

pdf("Handcoded MCMC Simulations/BYM2 spatial effect SRS Alg3.1 230814, posterior distributions.pdf") 
{
  hist(log(postsamp_3.1$tau))
  abline(v=log(tau),col='red')
  hist(logitlink(postsamp_3.1$phi))
  abline(v=logitlink(phi),col='red')
  hist(postsamp_3.1$alpha)
  abline(v=alpha,col='red')
  
  results3.1_b <- gather(data.frame(postsamp_3.1$b),'area','b')
  results3.1_b$area <- as.numeric(str_remove(results3.1_b$area,'X'))
  true_b <- data.frame(area=1:n_areas,true_b=as.vector(b))
  true_b <- true_b[order(true_b$true_b),]
  true_b$b_rank <- 1:nrow(true_b)
  results3.1_b <- left_join(results3.1_b,true_b,by='area')
  g <- results3.1_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~b_rank) + ggtitle('5000 births sampled per area') +
    geom_hline(yintercept=0,col='grey50')
  print(g)
  
  results3.1_Y <- gather(data.frame(postsamp_3.1$Y),'area','Y')
  results3.1_Y$area <- as.numeric(str_remove(results3.1_Y$area,'X'))
  true_Y <- data.frame(area=1:n_areas,true_Y=obs_dat$Y)
  results3.1_Y <- left_join(results3.1_Y,true_Y,by='area')
  g <- results3.1_Y %>% ggplot() + geom_boxplot(aes(y=Y)) +
    geom_hline(aes(yintercept = true_Y),col='red') + facet_grid(~area) + ggtitle('5000 births sampled per area')
  print(g)
  
  plot(log(postsamp_3.1$tau)[1:50000%%50==0],type='l')
  abline(h=log(tau),col='red')
  plot(logitlink(postsamp_3.1$phi[1:50000%%50==0]),type='l')
  abline(h=logitlink(phi),col='red')
  plot(postsamp_3.1$alpha[1:50000%%50==0],type='l')
  abline(h=alpha,col='red')
  
  plot(postsamp_3.1$b[,6],type='l')
  abline(h=b[6],col='red')
  plot(postsamp_3.1$b[,15],type='l')
  abline(h=b[15],col='red')
  plot(postsamp_3.1$b[,42],type='l')
  abline(h=b[42],col='red')
  
  plot(postsamp_3.1$Y[,6],type='l')
  abline(h=obs_dat$Y[6],col='red')
  plot(postsamp_3.1$Y[,15],type='l')
  abline(h=obs_dat$Y[15],col='red')
  plot(postsamp_3.1$Y[,42],type='l')
  abline(h=obs_dat$Y[42],col='red')
  
} 
dev.off()  

## Mixes fairly well after 50k iterations
postsamp_3.2 <- postsamp_Yartauphi_ZYplus_MCMC(tau_init = 25,phi_init = 0.25,alpha_init = -6,
                                               prop_sd_tau = 0.25, prop_sd_phi = 0.25,
                                               prop_Y_f_max = 5, data = obs_dat, Q_inv = Q_scaled_inv, n_iter = 50000)

pdf("Handcoded MCMC Simulations/BYM2 spatial effect SRS Alg3.2 230830 dataset 2, block updates for BYM2, 50k iterations posterior distributions.pdf") 
{
  hist(log(postsamp_3.2$tau))
  abline(v=log(tau),col='red')
  hist(logitlink(postsamp_3.2$phi))
  abline(v=logitlink(phi),col='red')
  hist(postsamp_3.2$alpha)
  abline(v=alpha,col='red')
  
  results3.2_b <- gather(data.frame(postsamp_3.2$b),'area','b')
  results3.2_b$area <- as.numeric(str_remove(results3.2_b$area,'X'))
  true_b <- data.frame(area=1:n_areas,true_b=as.vector(b))
  true_b <- true_b[order(true_b$true_b),]
  true_b$b_rank <- 1:nrow(true_b)
  results3.2_b <- left_join(results3.2_b,true_b,by='area')
  g <- results3.2_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~b_rank) + ggtitle('5000 births sampled per area') +
    geom_hline(yintercept=0,col='grey50')
  print(g)
  
  results3.2_Y <- gather(data.frame(postsamp_3.2$Y),'area','Y')
  results3.2_Y$area <- as.numeric(str_remove(results3.2_Y$area,'X'))
  true_Y <- data.frame(area=1:n_areas,true_Y=obs_dat$Y)
  results3.2_Y <- left_join(results3.2_Y,true_Y,by='area')
  g <- results3.2_Y %>% ggplot() + geom_boxplot(aes(y=Y)) +
    geom_hline(aes(yintercept = true_Y),col='red') + facet_grid(~area) + ggtitle('5000 births sampled per area')
  print(g)
  
  plot(log(postsamp_3.2$tau[1:50000%%50==0]),type='l')
  abline(h=log(tau),col='red')
  plot(logitlink(postsamp_3.2$phi[1:50000%%50==0]),type='l')
  abline(h=logitlink(phi),col='red')
  plot(postsamp_3.2$alpha[1:50000%%50==0],type='l')
  abline(h=alpha,col='red')
  
  plot(postsamp_3.2$b[1:50000%%50==0,6],type='l')
  abline(h=b[6],col='red')
  plot(postsamp_3.2$b[1:50000%%50==0,15],type='l')
  abline(h=b[15],col='red')
  plot(postsamp_3.2$b[1:50000%%50==0,27],type='l')
  abline(h=b[27],col='red')
  
  plot(postsamp_3.2$Y[1:50000%%50==0,6],type='l')
  abline(h=obs_dat$Y[6],col='red')
  plot(postsamp_3.2$Y[1:50000%%50==0,15],type='l')
  abline(h=obs_dat$Y[15],col='red')
  plot(postsamp_3.2$Y[1:50000%%50==0,27],type='l')
  abline(h=obs_dat$Y[27],col='red')
} 
dev.off()  

# not great mixing after 50k iterations, better after 100k -- probably good enough because estimation looks good
postsamp_3.3 <- postsamp_Alg3.3_MCMC(tau_init = 25,phi_init = 0.25,alpha_init = -6,
                                     prop_sd_tau = 0.5, prop_sd_phi = 0.25,
                                     prop_Y_f_max = 2, data = obs_dat, Q_inv = Q_scaled_inv, n_iter = 100000)

pdf("Handcoded MCMC Simulations/BYM2 spatial effect SRS Alg3.3 230901 dataset 1, block updates for BYM2, 100k iterations posterior distributions.pdf") 
{
  hist(log(postsamp_3.3$tau))
  abline(v=log(tau),col='red')
  hist(logitlink(postsamp_3.3$phi))
  abline(v=logitlink(phi),col='red')
  hist(postsamp_3.3$alpha)
  abline(v=alpha,col='red')
  
  results3.3_b <- gather(data.frame(postsamp_3.3$b),'area','b')
  results3.3_b$area <- as.numeric(str_remove(results3.3_b$area,'X'))
  true_b <- data.frame(area=1:n_areas,true_b=as.vector(b))
  true_b <- true_b[order(true_b$true_b),]
  true_b$b_rank <- 1:nrow(true_b)
  results3.3_b <- left_join(results3.3_b,true_b,by='area')
  g <- results3.3_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~b_rank) + ggtitle('5000 births sampled per area') +
    geom_hline(yintercept=0,col='grey50')
  print(g)
  
  results3.3_Y <- gather(data.frame(postsamp_3.3$Y),'area','Y')
  results3.3_Y$area <- as.numeric(str_remove(results3.3_Y$area,'X'))
  true_Y <- data.frame(area=1:n_areas,true_Y=obs_dat$Y)
  results3.3_Y <- left_join(results3.3_Y,true_Y,by='area')
  g <- results3.3_Y %>% ggplot() + geom_boxplot(aes(y=Y)) +
    geom_hline(aes(yintercept = true_Y),col='red') + facet_grid(~area) + ggtitle('5000 births sampled per area')
  print(g)
  
  plot(log(postsamp_3.3$tau[1:100000%%100==0]),type='l')
  abline(h=log(tau),col='red')
  plot(logitlink(postsamp_3.3$phi[1:100000%%100==0]),type='l')
  abline(h=logitlink(phi),col='red')
  plot(postsamp_3.3$alpha[1:100000%%100==0],type='l')
  abline(h=alpha,col='red')
  
  plot(postsamp_3.3$b[1:100000%%100==0,6],type='l')
  abline(h=b[6],col='red')
  plot(postsamp_3.3$b[1:100000%%100==0,15],type='l')
  abline(h=b[15],col='red')
  plot(postsamp_3.3$b[1:100000%%100==0,27],type='l')
  abline(h=b[27],col='red')
  
  plot(postsamp_3.3$Y[1:100000%%100==0,6],type='l')
  abline(h=obs_dat$Y[6],col='red')
  plot(postsamp_3.3$Y[1:100000%%100==0,15],type='l')
  abline(h=obs_dat$Y[15],col='red')
  plot(postsamp_3.3$Y[1:100000%%100==0,27],type='l')
  abline(h=obs_dat$Y[27],col='red')
} 
dev.off() 


### get distribution of posterior mean for several data realizations
{
  nsims <- 100
  results3.3_tau <- rep(NA,nsims)
  results3.3_phi <- rep(NA,nsims)
  results3.3_alpha <- rep(NA,nsims)
  results3.3_eta <- matrix(NA,nsims,n_areas)
  results3.3_Y <- matrix(NA,nsims,n_areas)
  for(k in 1:nsims){
    #status check
    if(k %% (nsims/10)==0){cat('Starting simulation', k, '\n')}
    
    #generate data
    all_dat <- simdat_bym2sp(n_areas=n_areas,n_clusters=n_clusters,n_births=n_births,a=alpha,Var_b=Var_b,b=b)$dat
    obs_dat <- get_srs(all_dat,n_births_samp,'A')
    
    #run MCMC
    postsamp_3.3 <- postsamp_Alg3.3_MCMC(tau_init = 25,phi_init = 0.25,alpha_init = -6,
                                         prop_sd_tau = 0.5, prop_sd_phi = 0.25,
                                         prop_Y_f_max = 2, data = obs_dat, Q_inv = Q_scaled_inv, 
                                         n_iter = 100000)
    
    #record result
    results3.3_tau[k] <- mean(postsamp_3.3$tau)
    results3.3_phi[k] <- mean(postsamp_3.3$phi)
    results3.3_alpha[k] <- mean(postsamp_3.3$alpha)
    results3.3_eta[k,] <- colMeans(postsamp_3.3$eta)
    results3.3_Y[k,] <- colMeans(postsamp_3.3$Y)
  }
  
  results_3.3 <- list(tau_res = results3.3_tau, phi_res = results3.3_phi,
                      alpha_res = results3.3_alpha,eta_res = results3.3_eta, Y_res  = results3.3_Y,
                      tau_truth = tau, phi_truth=phi, alpha_truth=alpha, b_truth = b)
  
  save(results_3.3,file = 'Handcoded MCMC Simulations/Alg 3.3 230901 100k iterations for 100 datasets.rda')
}

pdf("Handcoded MCMC Simulations/BYM2 spatial effect SRS Alg 3.3 230901 same b, block updates for BYM2, 100k iterations.pdf")
{
  hist(log(results3.3_tau))
  abline(v=log(tau),col='red')
  hist(logitlink(results3.3_phi))
  abline(v=logitlink(phi),col='red')
  hist(results3.3_alpha)
  abline(v=alpha,col='red')
  
  results_eta <- gather(data.frame(results3.3_eta),'area','eta')
  results_eta$area <- as.numeric(str_remove(results_eta$area,'X'))
  results_eta$b <- results_eta$eta - results3.3_alpha
  true_b <- data.frame(area=1:n_areas,true_b=as.vector(b))
  results_b <- left_join(results_eta,true_b,by='area')
  g <- results_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~area) + ggtitle('Posterior means from 100 simulated data sets with same b')
  print(g)
  
  results_Y <- gather(data.frame(results3.3_Y),'area','Y')
  results_Y$area <- as.numeric(str_remove(results_Y$area,'X'))
  true_Y <- data.frame(area=1:n_areas,true_Y_mean=n_births*n_clusters*exp(alpha + as.vector(b)))
  results_Y <- left_join(results_Y,true_Y,by='area')
  g <- results_Y %>% ggplot() + geom_boxplot(aes(y=Y)) +
    geom_hline(aes(yintercept = true_Y_mean),col='red') + facet_grid(~area) + ggtitle('Posterior means from 100 simulated data sets with same b')
  print(g)
}
dev.off()

#######################################################################
#######################################################################
# Alg 4.1: P(alphaU,alphaR,b,tau,d|Y+,Y^s,delta) ------

grad_log_target_4.1_logd <- function(data,Yplus,pop_dat,alphaU,alphaR,b,d,mu_d,sd_d){
  Nr_data <- data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A])
  Nr_sum_obs <- sum(Nr_data)
  Nr_sum <- sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  term1 <- (Yplus - Nr_sum)/(d+1)
  term2 <- -(1/d)*(Nr_sum-Nr_sum_obs)*(digamma(Yplus - sum(data$Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - digamma((1/d)*(Nr_sum-Nr_sum_obs)) - log(1+d))
  term3 <- -sum(1/d*Nr_data*(digamma(data$Y+1/d*Nr_data) - digamma(1/d*Nr_data) - log(1+d)))
  term4 <- -(log(d)-mu_d)/sd_d^2
  return(term1 + term2 + term3 + term4)
}

grad_log_target_4.1_b <- function(data,Yplus,pop_dat,area,alphaU,alphaR,b,d,tau){
  Nr_sum_obs <- sum(data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A]))
  Nr_sum <- sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  Nr_sum_area_obs <- sum(I(data$A==area)*data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[area]))
  Nr_sum_area <- sum(I(pop_dat$A==area)*pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[area]))
  
  term1 <- (1/d)*(Nr_sum_area-Nr_sum_area_obs)*(digamma(Yplus - sum(data$Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - digamma((1/d)*(Nr_sum-Nr_sum_obs)) - log(d+1))
  Nr_tmp <- data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[area])
  term2 <- sum(I(data$A==area)*1/d*Nr_tmp*(digamma(data$Y+1/d*Nr_tmp) - digamma(1/d*Nr_tmp) - log(1+d)),na.rm = T)
  term3 <- -tau*b[area]
  return(term1 + term2 + term3)
}

grad_log_target_4.1_alphaU <- function(data,Yplus,pop_dat,alphaU,alphaR,b,d,mu_U,sd_U){
  Nr_sum_obs <- sum(data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A]))
  Nr_sum <- sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  Nr_sum_urb_obs <- sum(data$U*data$N*exp(data$U*alphaU + b[data$A]))
  Nr_sum_urb <- sum(pop_dat$U*pop_dat$N*exp(pop_dat$U*alphaU + b[pop_dat$A]))
  
  term1 <- (1/d)*(Nr_sum_urb-Nr_sum_urb_obs)*(digamma(Yplus - sum(data$Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - digamma((1/d)*(Nr_sum-Nr_sum_obs)) - log(d+1))
  Nr_tmp <- data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A])
  term2 <- sum(data$U*1/d*Nr_tmp*(digamma(data$Y+1/d*Nr_tmp) - digamma(1/d*Nr_tmp) - log(1+d)),na.rm = T)
  term3 <- -(alphaU - mu_U)/sd_U^2
  return(term1 + term2 + term3)
  
}

grad_log_target_4.1_alphaR <- function(data,Yplus,pop_dat,alphaU,alphaR,b,d,mu_R,sd_R){
  Nr_sum_obs <- sum(data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A]))
  Nr_sum <- sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  Nr_sum_rur_obs <- sum((1-data$U)*data$N*exp((1-data$U)*alphaR + b[data$A]))
  Nr_sum_rur <- sum((1-pop_dat$U)*pop_dat$N*exp((1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  term1 <- (1/d)*(Nr_sum_rur-Nr_sum_rur_obs)*(digamma(Yplus - sum(data$Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - digamma((1/d)*(Nr_sum-Nr_sum_obs)) - log(d+1))
  Nr_tmp <- data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A])
  term2 <- sum((1-data$U)*1/d*Nr_tmp*(digamma(data$Y+1/d*Nr_tmp) - digamma(1/d*Nr_tmp) - log(1+d)),na.rm = T)
  term3 <- -(alphaR - mu_R)/sd_R^2
  return(term1 + term2 + term3)
  
}

postsamp_Alg4.1_MCMC <- function(tau_init, d_init, alphaU_init, alphaR_init, #initial values
                                 prop_sd_d=NULL, prop_sd_alpha=NULL, prop_sd_b=NULL, # for random walk proposals
                                 prop_k_alpha=NULL, prop_k_logd=NULL, prop_k_b=NULL, # for MALA proposals
                                 data_list, n_iter,burnin=n_iter/2,chains=2){
  
  data <- data_list$obs_dat
  Yplus <- data_list$Yplus
  pop_dat <- data_list$pop_strata_dat
  n_areas <- length(unique(data_list$pop_strata_dat$A))
  #hyperpriors
  a_tau <- 0.01
  b_tau <- 0.01
  mu_U <- -3.5
  sd_U <- 3
  mu_R <- -3
  sd_R <- 3
  mu_d <- 0
  sd_d <- 1
  
  #set initial values
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  d_postsamp_t <- c(d_init, rep(NA,n_iter-1))
  alphaU_postsamp_t <- c(alphaU_init, rep(NA,n_iter-1))
  alphaR_postsamp_t <- c(alphaR_init, rep(NA,n_iter-1))
  
  b_init <- rnorm(n_areas,0,1/sqrt(10*tau_init))
  b_postsamp_t <- rbind(b_init, matrix(NA,n_iter-1,n_areas))
  
  tau_postsamp <- d_postsamp <- alphaU_postsamp <- alphaR_postsamp <- b_postsamp <- NULL
  acc <- acc_b <- acc_alpha <- acc_d <- 0
  
  # set current states
  tau_current <- tau_init
  d_current <- d_init
  alphaU_current <- alphaU_init
  alphaR_current <- alphaR_init
  b_current <- b_init
  
  for(chain in 1:chains){
    for(i in 2:n_iter){
      
      # draw new tau from full conditional
      tau_current <- rgamma(1,n_areas/2+a_tau,0.5*sum(b_current^2)+b_tau)
      
      # draw proposal for b 
      # random walk
      # b_proposal <- sapply(b_current, function(x){(rnorm(1,x,prop_sd_b))})
      
      b_proposal <- b_current +
        sapply(1:n_areas,function(A){prop_k_b*grad_log_target_4.1_b(data,Yplus,pop_dat,A,alphaU_current,alphaR_current,b_current,d_current,tau_current)}) +
        sqrt(2*prop_k_b)*rnorm(n_areas)

      # draw proposal for d
      #d_proposal <- rlnorm(1,log(d_current),prop_sd_d)
      logd_proposal <- log(d_current) + prop_k_logd*grad_log_target_4.1_logd(data,Yplus,pop_dat,alphaU_current,alphaR_current,b_current,d_current,mu_d,sd_d) +
        sqrt(2*prop_k_logd)*rnorm(1)
      d_proposal <- exp(logd_proposal)
      
      # draw proposal for alphas
      # alphaU_proposal <- rnorm(1,alphaU_current,prop_sd_alpha)
      # alphaR_proposal <- rnorm(1,alphaR_current,prop_sd_alpha)
      alphaU_proposal <- alphaU_current + prop_k_alpha*grad_log_target_4.1_alphaU(data,Yplus,pop_dat,alphaU_current,alphaR_current,b_current,d_current,mu_U,sd_U) +
        sqrt(2*prop_k_alpha)*rnorm(1)
      alphaR_proposal <- alphaR_current + prop_k_alpha*grad_log_target_4.1_alphaR(data,Yplus,pop_dat,alphaU_current,alphaR_current,b_current,d_current,mu_R,sd_R) +
        sqrt(2*prop_k_alpha)*rnorm(1)

      #calculate the alphas parameter for the Dir-Multinom distributions
      Nr_obs_proposal <- sapply(1:nrow(data),function(i){data$N[i]*exp(alphaU_proposal*data$U[i]+alphaR_proposal*(1-data$U[i])+b_proposal[data$A[i]])})
      Nr_obs_current <- sapply(1:nrow(data),function(i){data$N[i]*exp(alphaU_current*data$U[i]+alphaR_current*(1-data$U[i])+b_current[data$A[i]])})
      
      Nr_sum_proposal <- sum(pop_dat$N*exp(alphaU_proposal*pop_dat$U+alphaR_proposal*(1-pop_dat$U)+b_proposal[pop_dat$A]))
      Nr_sum_current <- sum(pop_dat$N*exp(alphaU_current*pop_dat$U+alphaR_current*(1-pop_dat$U)+b_current[pop_dat$A]))
      # accept or reject proposal for (alpha, b, d)
      a_prob <- exp({
        # Y|Y+,alpha,b
        lgamma(Nr_sum_proposal/d_proposal) - lgamma(Yplus + Nr_sum_proposal/d_proposal) + 
          lgamma(Yplus - sum(data$Y) + (Nr_sum_proposal - sum(Nr_obs_proposal))/d_proposal) -  lgamma((Nr_sum_proposal - sum(Nr_obs_proposal))/d_proposal) +
          sum(lgamma(data$Y+Nr_obs_proposal/d_proposal)) -  sum(lgamma(Nr_obs_proposal/d_proposal)) - 
        (lgamma(Nr_sum_current/d_current) - lgamma(Yplus + Nr_sum_current/d_current) + 
           lgamma(Yplus - sum(data$Y) + (Nr_sum_current - sum(Nr_obs_current))/d_current) -  lgamma((Nr_sum_current - sum(Nr_obs_current))/d_current) +
           sum(lgamma(data$Y+Nr_obs_current/d_current)) -  sum(lgamma(Nr_obs_current/d_current))) +
          # Yplus|alpha,b,d
          dnbinom(Yplus, size=1/d_proposal*Nr_sum_proposal, prob=1/(1+d_proposal),log=T) - 
          dnbinom(Yplus, size=1/d_current*Nr_sum_current, prob=1/(1+d_current),log=T) +
          # b|tau
          sum(sapply(1:n_areas,function(i){
            dnorm(b_proposal[i],0,sqrt(1/tau_current),log=T) - dnorm(b_current[i],0,sqrt(1/tau_current),log=T)})) +
          # prior on alphas
          dnorm(alphaU_proposal,mu_U,sd_U,log=T) - dnorm(alphaU_current,mu_U,sd_U,log=T) +
          dnorm(alphaR_proposal,mu_R,sd_R,log=T) - dnorm(alphaR_current,mu_R,sd_R,log=T) +
          # prior on d
          dlnorm(d_proposal,mu_d,sd_d,log=T) - dlnorm(d_current,mu_d,sd_d,log=T) +
          # proposal dist for alpha (MALA)
          dnorm(alphaU_current,alphaU_proposal + prop_k_alpha*grad_log_target_4.1_alphaU(data,Yplus,pop_dat,alphaU_proposal,alphaR_current,b_current,d_current,mu_U,sd_U),sqrt(2*prop_k_alpha),log=T) -
          dnorm(alphaU_proposal,alphaU_current + prop_k_alpha*grad_log_target_4.1_alphaU(data,Yplus,pop_dat,alphaU_current,alphaR_current,b_current,d_current,mu_U,sd_U),sqrt(2*prop_k_alpha),log=T) +
          dnorm(alphaR_current,alphaR_proposal + prop_k_alpha*grad_log_target_4.1_alphaR(data,Yplus,pop_dat,alphaU_current,alphaR_proposal,b_current,d_current,mu_R,sd_R),sqrt(2*prop_k_alpha),log=T) -
          dnorm(alphaR_proposal,alphaR_current + prop_k_alpha*grad_log_target_4.1_alphaR(data,Yplus,pop_dat,alphaU_current,alphaR_current,b_current,d_current,mu_R,sd_R),sqrt(2*prop_k_alpha),log=T) +
          # MALA proposal dist for b
          sum(vapply(1:n_areas,function(i){
            dnorm(b_current[i],b_proposal[i] + prop_k_b*grad_log_target_4.1_b(data,Yplus,pop_dat,i,alphaU_current,alphaR_current,b_proposal,d_current,tau_current),sqrt(2*prop_k_b),log=T) -
              dnorm(b_proposal[i],b_current[i] + prop_k_b*grad_log_target_4.1_b(data,Yplus,pop_dat,i,alphaU_current,alphaR_current,b_current,d_current,tau_current),sqrt(2*prop_k_b),log=T)
          },numeric(1))) +
          # random walk proposal dist for d
          #dlnorm(d_current,log(d_proposal),prop_sd_d,log=T) - dlnorm(d_proposal,log(d_current),prop_sd_d,log=T)
          # MALA proposal dist for d
          dlnorm(d_current,log(d_proposal) + prop_k_logd*grad_log_target_4.1_logd(data,Yplus,pop_dat,alphaU_current,alphaR_current,b_current,d_proposal,mu_d,sd_d),sqrt(2*prop_k_logd),log=T) -
          dlnorm(d_proposal,log(d_current) + prop_k_logd*grad_log_target_4.1_logd(data,Yplus,pop_dat,alphaU_current,alphaR_current,b_current,d_current,mu_d,sd_d),sqrt(2*prop_k_logd),log=T)

      })
      
      if(a_prob>runif(1)){
        d_current <- d_proposal
        b_current <- b_proposal
        alphaU_current <- alphaU_proposal
        alphaR_current <- alphaR_proposal
        if(i>burnin){acc <- acc + 1}
      }
      
      # record current state
      tau_postsamp_t[i] <- tau_current
      d_postsamp_t[i] <- d_current
      alphaU_postsamp_t[i] <- alphaU_current
      alphaR_postsamp_t[i] <- alphaR_current
      b_postsamp_t[i,] <- b_current
      
    } 
    
    #get rid of burn-in and add to other chains
    b_postsamp <- rbind(b_postsamp,b_postsamp_t[(burnin+1):n_iter,])
    alphaU_postsamp <- c(alphaU_postsamp,alphaU_postsamp_t[(burnin+1):n_iter])
    alphaR_postsamp <- c(alphaR_postsamp,alphaR_postsamp_t[(burnin+1):n_iter])
    tau_postsamp <- c(tau_postsamp,tau_postsamp_t[(burnin+1):n_iter])
    d_postsamp <- c(d_postsamp,d_postsamp_t[(burnin+1):n_iter])
  }
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp, alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  acc_rate <- acc/((n_iter-burnin)*chains)
  #acc_b_rate <- acc_b/((n_iter-burnin)*chains)
  #acc_alpha_rate <- acc_alpha/((n_iter-burnin)*chains)
  #acc_d_rate <- acc_d/((n_iter-burnin)*chains)
  
  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              tau=tau_postsamp, d=d_postsamp, 
              acc_rate = acc_rate))
  #acc_b_rate=acc_b_rate,acc_d_rate=acc_d_rate,acc_alpha_rate=acc_alpha_rate))
}


#######################################################################
#######################################################################
# Simulations for Alg 4.1 ------

## population/parameter settings
{
  n_areas <- 20 #numbers of areas
  #n_clusters_urban <- rnegbin(n_areas,250,5) #number of urban clusters in each area
  #n_clusters_rural <- rnegbin(n_areas,250,5) #number of rural clusters in each area
  n_clusters_urban <- rnegbin(n_areas,250,5) #number of urban clusters in each area
  n_clusters_rural <- rnegbin(n_areas,250,5) #number of rural clusters in each area
  n_births_urban <- 30 #average number of births per urban cluster
  n_births_rural <- 20 #average number of births per rural cluster
  n_clusters_urban_samp <-  round(0.1*n_clusters_urban)
  n_clusters_rural_samp <-  round(0.1*n_clusters_rural)
  # n_births_urban_samp <- 25 # number of births sampled from each urban cluster
  # n_births_rural_samp <- 20 # number of births sampled from each rural cluster
  
  # intercept (fixed effect) 
  alphaU <- rnorm(1,-3.5,0.25)
  alphaR <- rnorm(1,-3,0.25)
  
  # draw hyperparameters
  tau <- rgamma(1,10,.25)
  
  # draw overdispersion parameter
  d <- rlnorm(1,0,0.25)
  
  # draw random effects
  b <- as.vector(Rfast::rmvnorm(1,rep(0,n_areas),diag(1/tau,n_areas)))
}

## generate data 
{
  all_dat <- data.frame(
    # (true) number of births in each cluster
    N = c(rpois(sum(n_clusters_urban),n_births_urban),rpois(sum(n_clusters_rural),n_births_rural)),
    # urban or rural strata (U=1 if urban, otherwise 0)
    U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
    # admin area of each cluster
    A = c(unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_rural[i])}))))
  all_dat$cluster <- 1:nrow(all_dat)
  all_dat$Y <- sapply(1:nrow(all_dat),function(i){rnbinom(1,
                                                          size = 1/d*all_dat$N[i]*exp(alphaU*all_dat$U[i] + alphaR*(1-all_dat$U[i]) + b[all_dat$A[i]] + rnorm(1,0,0.05)),
                                                          prob=1/(1+d))})
  
  ### sample clusters
  obs_dat <- NULL
    for(area in 1:n_areas){
      # randomly select clusters from urban strata based on size of cluster
      obs_dat <- rbind(obs_dat,all_dat[all_dat$A==area & all_dat$U==1,][sample(1:n_clusters_urban[area],n_clusters_urban_samp[area],prob = all_dat[all_dat$U==1 & all_dat$A==area,]$N),],
                       # randomly select cluster from rural strata  based on size of cluster
                       all_dat[all_dat$A==area & all_dat$U==0,][sample(1:n_clusters_rural[area],n_clusters_rural_samp[area],prob = all_dat[all_dat$U==0& all_dat$A==area,]$N),])
    }
  
  
  data_list <- list(obs_dat = obs_dat, #observed data
                    Yplus = sum(all_dat$Y), #ALL deaths
                    pop_strata_dat = all_dat %>% group_by(A,U) %>% summarise(N=sum(N))) # number of births per strata (area x urban)
}

## check data generation
{
  # Yplus <- NULL
  # for(i in 1:500){
  #   all_dat$Y <- sapply(1:nrow(all_dat),function(i){rnbinom(1,
  #                                                           size = 1/d*all_dat$N[i]*exp(alphaU*all_dat$U[i] + alphaR*(1-all_dat$U[i]) + b[all_dat$A[i]] + rnorm(1,0,0.05)),
  #                                                           prob=1/(1+d))})
  #   Yplus[i] <- sum(all_dat$Y)
  # }
  # 
  # # # theoretical mean:
  # sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  # # sample mean:
  # mean(Yplus)
  # 
  # # theoretical variance:
  # (d+1)*sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  # # sample variance:
  # var(Yplus)
}

## run model
# very sensitive to choice of k (particularly for alpha)
# if too large, for certain data realizations, acceptance rate drops straight to 0
# always run for 1000 iterations first to properly tune ks
postsamp_4.1 <- postsamp_Alg4.1_MCMC(tau_init = 40, d_init = 1, alphaU_init = -3.5, alphaR_init = -3,
                                   #  prop_sd_alpha = 0.05, 
                                    # prop_sd_b = 0.01,
                                    # prop_sd_d = 0.1,
                                     prop_k_alpha = 2e-4, 
                                     prop_k_b = 2e-4,
                                     prop_k_logd = 0.01, 
                                     data_list = data_list, n_iter=1000)

pdf("/Users/alanamcgovern/Desktop/Research/New Benchmarking/Handcoded MCMC Simulations/Alg 4.1 231128, MALA block updates except tau, 50k iterations posterior distributions.pdf") 
{
  hist(log(postsamp_4.1$tau))
  abline(v=log(tau),col='red')
  hist(log(postsamp_4.1$d))
  abline(v=log(d),col='red')
  
  hist(postsamp_4.1$alphaU)
  abline(v=alphaU,col='red')
  hist(postsamp_4.1$alphaR)
  abline(v=alphaR,col='red')
  
  results4.1_b <- gather(data.frame(postsamp_4.1$b),'area','b')
  results4.1_b$area <- as.numeric(str_remove(results4.1_b$area,'X'))
  true_b <- data.frame(area=1:n_areas,true_b=as.vector(b))
  true_b <- true_b[order(true_b$true_b),]
  true_b$b_rank <- 1:nrow(true_b)
  results4.1_b <- left_join(results4.1_b,true_b,by='area')
  g <- results4.1_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    geom_hline(aes(yintercept = true_b),col='red') + facet_grid(~b_rank) + ggtitle('Approx 2500 births per area') +
    geom_hline(yintercept=0,col='grey50')
  print(g)
  
  plot(log(postsamp_4.1$tau[1:50000%%50==0]),type='l')
  abline(h=log(tau),col='red')
  plot(log(postsamp_4.1$d[1:50000%%50==0]),type='l')
  abline(h=log(d),col='red')
  plot(postsamp_4.1$alphaU[1:50000%%50==0],type='l')
  abline(h=alphaU,col='red')
  plot(postsamp_4.1$alphaR[1:50000%%50==0],type='l')
  abline(h=alphaR,col='red')
  plot(postsamp_4.1$b[1:50000%%50==0,2],type='l')
  abline(h=b[2],col='red')
  plot(postsamp_4.1$b[1:50000%%50==0,7],type='l')
  abline(h=b[7],col='red')
  plot(postsamp_4.1$b[1:50000%%50==0,15],type='l')
  abline(h=b[15],col='red')
  plot(postsamp_4.1$b[1:50000%%50==0,20],type='l')
  abline(h=b[20],col='red')
} 
dev.off()  

#######################################################################
#######################################################################
# Alg 4.1 with real data! --------------
## Because we are using 1 stage cluster sampling, we are assuming the only households in a sampled cluster are the ones we observed
# and each unsampled cluster has same number of households that were observed in sampled clusters

#load in data
{
  load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Sierra_Leone/Sierra_Leone_cluster_dat.rda")
  admin2_key <- mod.dat %>% dplyr::select(admin2,admin2.char,admin2.name) %>% unique()
  
  obs_dat <- mod.dat %>% filter(age==0,years==2010,survey==2019) %>% group_by(cluster) %>% 
    reframe(N=sum(total),
            U=as.numeric(I(urban=='urban')),
            urban=urban,
            A=unique(admin2),
            A.name = unique(admin2.name),
            wt = v005,
            Y = sum(Y)) %>% unique()
  
  
  #calculate (obs) births per household in each strata
  pop_dat <- obs_dat %>% group_by(A,A.name,U) %>% summarise(sum(N)) %>% arrange(A,U) %>% rename(admin2.name = A.name)
  # from DHS report -- number of observed HH
  pop_dat$obsHH <- c(628, 364, 498, 390,
                     572, 130, 576, 312, 
                     472, 260, 628, 52,
                     524, 468, 498, 52,
                     468, 208, 646, 286, 
                     732, 78, 628, 312,
                     654, 78, 680, 208,
                     104, 936, 1430)
  pop_dat$bHHratio <- pop_dat$`sum(N)`/pop_dat$obsHH
  
  # from DHS report -- total number of clusters
  pop_dat$nclusters <- c(708, 325,464,266,
                         390,71, 615,276,
                         376, 200, 409, 26,
                         678, 441, 330, 24,
                         271, 123, 586, 201,
                         579, 37, 706, 294,
                         549, 33, 834, 207,
                         65, 635, 2139)
  # from DHS report -- choose 24 from each cluster
  pop_dat$N <- round(pop_dat$nclusters*24*pop_dat$bHHratio)
  pop_dat <- pop_dat %>% dplyr::select(A,admin2.name,U,N)
  
  igme.nmr <- read_csv("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/IGME/igme2022_nmr.csv")
  igme.est <- (igme.nmr %>% filter(iso=='SLE',Quantile=='Median'))$`2010.5`/1000
  Yplus_est <- round(igme.est*sum(pop_dat$N))
  
  data_list <- list(obs_dat = obs_dat, #observed data
                    Yplus = Yplus_est, #ALL deaths
                    pop_strata_dat = pop_dat) # number of births per strata (area x urban)
}

postsamp_4.1 <- postsamp_Alg4.1_MCMC(tau_init = 25, d_init = 0.1, alphaU_init = -4, alphaR_init = -3,
                                     #prop_sd_alpha = 0.05, 
                                     #prop_sd_b = 0.05,
                                     prop_k_alpha = 2e-3, 
                                     prop_k_b = 5e-3,
                                     prop_k_logd = 0.05,
                                     data_list = data_list, n_iter=50000)

pdf("/Users/alanamcgovern/Desktop/Research/New Benchmarking/Handcoded MCMC Simulations/Alg 4.1 231031, SLE 2010, 50k iterations posterior distributions.pdf") 
{
  hist(log(postsamp_4.1$tau))
  hist(log(postsamp_4.1$d))
  
  hist(postsamp_4.1$alphaU)
  hist(postsamp_4.1$alphaR)
  
  results4.1_b <- gather(data.frame(postsamp_4.1$b),'area','b')
  results4.1_b$area <- as.numeric(str_remove(results4.1_b$area,'X'))
  g <- results4.1_b %>% ggplot() + geom_boxplot(aes(y=b)) +
    facet_grid(~area) + ggtitle('Approx 100 births per area') +
    geom_hline(yintercept=0,col='grey50')
  print(g)
  
  plot(log(postsamp_4.1$tau[1:50000%%50==0]),type='l')
  plot(log(postsamp_4.1$d[1:50000%%50==0]),type='l')
  plot(postsamp_4.1$alphaU[1:50000%%50==0],type='l')
  plot(postsamp_4.1$alphaR[1:50000%%50==0],type='l')
  plot(postsamp_4.1$b[1:50000%%50==0,2],type='l')
  plot(postsamp_4.1$b[1:50000%%50==0,7],type='l')
  plot(postsamp_4.1$b[1:50000%%50==0,15],type='l')
  
} 
dev.off()  


#compare to similar INLA models
{
  alg4.1.res <- expand_grid(admin2=1:nrow(admin2_key),U=c(0,1))
  alg4.1.res$eta <- median(postsamp_4.1$alphaU)*alg4.1.res$U + median(postsamp_4.1$alphaR)*(1-alg4.1.res$U) + Rfast::colMedians(postsamp_4.1$b)[alg4.1.res$admin2]
  alg4.1.res$model <- 'Alg 4.1'
  
  mod.dat.inla <- mod.dat %>% filter(age==0,years==2010,survey==2019)
  
  hyperpc1 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  inla.fit1 <- INLA::inla(Y ~ urban -1 +
                            f(admin2, model = "iid",hyper=hyperpc1),
                          data=mod.dat.inla, family='nbinomial', E=total,
                          control.predictor = list(compute = F, link = 1),
                          control.family = list(link='log'))
  inla.res1 <- expand_grid(admin2=1:nrow(admin2_key),U=c(0,1))
  inla.res1$eta <- inla.fit1$summary.fixed[1,1]*inla.res1$U + inla.fit1$summary.fixed[2,1]*(1-inla.res1$U) + inla.fit1$summary.random$admin2[inla.res1$admin2,2]
  inla.res1$model <- 'IID Adm2'
  
  inla.fit2 <- INLA::inla(Y ~ urban + factor(admin1) -1 +
                            f(admin2, model = "iid",hyper=hyperpc1),
                          data=mod.dat.inla, family='nbinomial', E=total,
                          control.predictor = list(compute = F, link = 1),
                          control.family = list(link='log'))
  inla.res2 <- expand_grid(admin2=1:nrow(admin2_key),U=c(0,1))
  inla.res2$eta <- inla.fit2$summary.fixed[1,1]*inla.res2$U + inla.fit2$summary.fixed[2,1]*(1-inla.res2$U) + inla.fit2$summary.random$admin2[inla.res2$admin2,2]
  inla.res2$model <- 'FE Adm1 + IID Adm2'
  
  all.res <- rbind(alg4.1.res,inla.res1,inla.res2)
  
  pdf("/Users/alanamcgovern/Desktop/Research/New Benchmarking/Handcoded MCMC Simulations/Compare Alg 4.1 to INLA 231031, SLE 2010.pdf") 
  
  all.res %>% ggplot(aes(factor(U),exp(eta),color=model)) + geom_point() + facet_grid(~admin2) + 
    ggtitle('Comparing results from Alg 4.1 MCMC to NB INLA models with urban intercept and \n various spatial effects using Sierra Leone 2010 data from 2019 survey')
  
  dev.off()
}

#######################################################################
#######################################################################
# Alg 5.1: P(alphaU,alphaR,b,tau,d|Y+,Y^s,delta) ------

grad_log_target_5.1_logd <- function(data,Yplus,pop_dat,Y,alphaU,alphaR,b,d,mu_d,sd_d){
  Nr_data <- data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A])
  Nr_sum_obs <- sum(Nr_data)
  Nr_sum <- sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  term1 <- (Yplus - Nr_sum)/(d+1)
  term2 <- -(1/d)*(Nr_sum-Nr_sum_obs)*(digamma(Yplus - sum(Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - digamma((1/d)*(Nr_sum-Nr_sum_obs)) - log(1+d))
  term3 <- -sum(1/d*Nr_data*(digamma(Y+1/d*Nr_data) - digamma(1/d*Nr_data) - log(1+d)))
  term4 <- -(log(d)-mu_d)/sd_d^2
  return(term1 + term2 + term3 + term4)
}

grad_log_target_5.1_b <- function(data,Yplus,pop_dat,Y,area,alphaU,alphaR,b,d,tau){
  Nr_sum_obs <- sum(data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A]))
  Nr_sum <- sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  Nr_sum_area_obs <- sum(I(data$A==area)*data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[area]))
  Nr_sum_area <- sum(I(pop_dat$A==area)*pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[area]))
  
  term1 <- (1/d)*(Nr_sum_area-Nr_sum_area_obs)*(digamma(Yplus - sum(Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - digamma((1/d)*(Nr_sum-Nr_sum_obs)) - log(1+d))
  Nr_tmp <- data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[area])
  term2 <- sum(I(data$A==area)*1/d*Nr_tmp*(digamma(Y+1/d*Nr_tmp) - digamma(1/d*Nr_tmp) - log(1+d)),na.rm = T)
  term3 <- -tau*b[area]
  return(term1 + term2 + term3)
  
}

grad_log_target_5.1_alphaU <- function(data,Yplus,pop_dat,Y,alphaU,alphaR,b,d,mu_U,sd_U){
  Nr_sum_obs <- sum(data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A]))
  Nr_sum <- sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  Nr_sum_urb_obs <- sum((data$U)*data$N*exp((data$U)*alphaU + b[data$A]))
  Nr_sum_urb <- sum(pop_dat$U*pop_dat$N*exp(pop_dat$U*alphaU + b[pop_dat$A]))
  
  term1 <- (1/d)*(Nr_sum_urb-Nr_sum_urb_obs)*(digamma(Yplus - sum(Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - digamma((1/d)*(Nr_sum-Nr_sum_obs)) - log(d+1))
  Nr_tmp <- data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A])
  term2 <- sum((data$U)*1/d*Nr_tmp*(digamma(Y+1/d*Nr_tmp) - digamma(1/d*Nr_tmp) - log(1+d)),na.rm = T)
  term3 <- -(alphaU - mu_U)/sd_U^2
  return(term1 + term2 + term3)
  
}

grad_log_target_5.1_alphaR <- function(data,Yplus,pop_dat,Y,alphaU,alphaR,b,d,mu_R,sd_R){
  Nr_sum_obs <- sum(data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A]))
  Nr_sum <- sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  Nr_sum_rur_obs <- sum((1-data$U)*data$N*exp((1-data$U)*alphaR + b[data$A]))
  Nr_sum_rur <- sum((1-pop_dat$U)*pop_dat$N*exp((1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  term1 <- (1/d)*(Nr_sum_rur-Nr_sum_rur_obs)*(digamma(Yplus - sum(Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - digamma((1/d)*(Nr_sum-Nr_sum_obs)) - log(d+1))
  Nr_tmp <- data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A])
  term2 <- sum((1-data$U)*1/d*Nr_tmp*(digamma(Y+1/d*Nr_tmp) - digamma(1/d*Nr_tmp) - log(1+d)),na.rm = T)
  term3 <- -(alphaR - mu_R)/sd_R^2
  return(term1 + term2 + term3)
  
}


# plot(data$Z[668]:round(0.25*data$N[668]),sapply(data$Z[668]:round(0.25*data$N[668]),function(x){
#   log_target_cond_5.1_y(x,data,Yplus,pop_dat,Y_current,668,Nr_sum_current,Nr_obs_current,d_current)}))

postsamp_Alg5.1_MCMC <- function(tau_init, d_init, alphaU_init, alphaR_init, #initial values
                                 b_init=NULL, Y_init=NULL,
                                 prop_sd_d=NULL, prop_sd_alpha=NULL, prop_sd_b=NULL, # for random walk proposals
                                 prop_k_alpha=NULL, prop_k_logd=NULL, prop_k_b=NULL, # for MALA proposals
                                 data_list, n_iter, burnin=n_iter/2,chains=2){
  
  data <- data_list$obs_dat
  Yplus <- data_list$Yplus
  pop_dat <- data_list$pop_strata_dat
  n_areas <- length(unique(data_list$pop_strata_dat$A))

  #hyperpriors
  a_tau <- 0.01
  b_tau <- 0.01
  mu_U <- -3.5
  sd_U <- 3
  mu_R <- -3
  sd_R <- 3
  mu_d <- 0
  sd_d <- 1
  
  #set initial values
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  d_postsamp_t <- c(d_init, rep(NA,n_iter-1))
  alphaU_postsamp_t <- c(alphaU_init, rep(NA,n_iter-1))
  alphaR_postsamp_t <- c(alphaR_init, rep(NA,n_iter-1))
  
  if(is.null(b_init))
    b_init <- rnorm(n_areas,0,1/sqrt(10*tau_init))
  b_postsamp_t <- rbind(b_init, matrix(NA,n_iter-1,n_areas))
  
  Nr_sum_init <- sum(pop_dat$N*exp(pop_dat$U*alphaU_init + (1-pop_dat$U)*alphaR_init + b_init[pop_dat$A]))
  Nr_obs_init <- data$N*exp(data$U*alphaU_init + (1-data$U)*alphaR_init + b_init[data$A])
  if(is.null(Y_init)){
    Y_init <- rmultinom(1,Yplus,prob=c(Nr_obs_init,Nr_sum_init-sum(Nr_obs_init)))
    Y_init <- Y_init[-length(Y_init)]
    Y_init[Y_init<data$Z] <- data$Z[Y_init<data$Z]
  }
  
  Y_postsamp_t <- rbind(Y_init, matrix(NA,n_iter-1,nrow(data)))
  
  tau_postsamp <- d_postsamp <- alphaU_postsamp <- alphaR_postsamp <- b_postsamp <- Y_postsamp <- NULL
  acc_reg <- 0
  acc_Y <- rep(0,nrow(data))
  
  # set current states
  tau_current <- tau_init
  d_current <- d_init
  alphaU_current <- alphaU_init
  alphaR_current <- alphaR_init
  b_current <- b_init
  Y_current <- Y_init
  Nr_sum_current <- Nr_sum_init
  Nr_obs_current <- Nr_obs_init
  
  parallel_results <- mclapply(1:chains, function(j){
    time.start <- Sys.time()
    for(i in 2:n_iter){
    
      # draw new tau from full conditional
      tau_current <- rgamma(1,n_areas/2+a_tau,0.5*sum(b_current^2)+b_tau)
      
      # draw proposal for b 
      # random walk
      # b_proposal <- sapply(b_current, function(x){(rnorm(1,x,prop_sd_b))})
      
      b_proposal <- b_current +
        prop_k_b*sapply(1:n_areas,function(A){grad_log_target_5.1_b(data,Yplus,pop_dat,Y_current,A,alphaU_current,alphaR_current,b_current,d_current,tau_current)}) +
        sqrt(2*prop_k_b)*rnorm(n_areas)
      
      # draw proposal for d
      #d_proposal <- rlnorm(1,log(d_current),prop_sd_d)
      logd_proposal <- log(d_current) + prop_k_logd*grad_log_target_5.1_logd(data,Yplus,pop_dat,Y_current,alphaU_current,alphaR_current,b_current,d_current,mu_d,sd_d) +
        sqrt(2*prop_k_logd)*rnorm(1)
      d_proposal <- exp(logd_proposal)
      
      # draw proposal for alphas
      # alphaU_proposal <- rnorm(1,alphaU_current,prop_sd_alpha)
      # alphaR_proposal <- rnorm(1,alphaR_current,prop_sd_alpha)
      alphaU_proposal <- alphaU_current + prop_k_alpha*grad_log_target_5.1_alphaU(data,Yplus,pop_dat,Y_current,alphaU_current,alphaR_current,b_current,d_current,mu_U,sd_U) +
        sqrt(2*prop_k_alpha)*rnorm(1)
      alphaR_proposal <- alphaR_current + prop_k_alpha*grad_log_target_5.1_alphaR(data,Yplus,pop_dat,Y_current,alphaU_current,alphaR_current,b_current,d_current,mu_R,sd_R) +
        sqrt(2*prop_k_alpha)*rnorm(1)
      
      #calculate the alphas parameter for the Dir-Multinom distributions
      Nr_obs_proposal <- sapply(1:nrow(data),function(i){data$N[i]*exp(alphaU_proposal*data$U[i]+alphaR_proposal*(1-data$U[i])+b_proposal[data$A[i]])})
      Nr_sum_proposal <- sum(pop_dat$N*exp(alphaU_proposal*pop_dat$U+alphaR_proposal*(1-pop_dat$U)+b_proposal[pop_dat$A]))
      
      # accept or reject proposal for (alpha, b, d)
      a_prob_reg <- exp({
        # Y|Y+,alpha,b
        lgamma(Nr_sum_proposal/d_proposal) - lgamma(Yplus + Nr_sum_proposal/d_proposal) + 
          lgamma(Yplus - sum(Y_current) + (Nr_sum_proposal - sum(Nr_obs_proposal))/d_proposal) -  lgamma((Nr_sum_proposal - sum(Nr_obs_proposal))/d_proposal) +
          sum(lgamma(Y_current+Nr_obs_proposal/d_proposal)) -  sum(lgamma(Nr_obs_proposal/d_proposal)) - 
          (lgamma(Nr_sum_current/d_current) - lgamma(Yplus + Nr_sum_current/d_current) + 
             lgamma(Yplus - sum(Y_current) + (Nr_sum_current - sum(Nr_obs_current))/d_current) -  lgamma((Nr_sum_current - sum(Nr_obs_current))/d_current) +
             sum(lgamma(Y_current+Nr_obs_current/d_current)) -  sum(lgamma(Nr_obs_current/d_current))) +
          # Yplus|alpha,b,d
          dnbinom(Yplus, size=1/d_proposal*Nr_sum_proposal, prob=1/(1+d_proposal),log=T) - 
          dnbinom(Yplus, size=1/d_current*Nr_sum_current, prob=1/(1+d_current),log=T) +
          # b|tau
          sum(dnorm(b_proposal,0,sqrt(1/tau_current),log=T) - dnorm(b_current,0,sqrt(1/tau_current),log=T)) +
          # prior on alphas
          dnorm(alphaU_proposal,mu_U,sd_U,log=T) - dnorm(alphaU_current,mu_U,sd_U,log=T) +
          dnorm(alphaR_proposal,mu_R,sd_R,log=T) - dnorm(alphaR_current,mu_R,sd_R,log=T) +
          # prior on d
          dlnorm(d_proposal,mu_d,sd_d,log=T) - dlnorm(d_current,mu_d,sd_d,log=T) +
          # proposal dist for alpha (MALA)
          dnorm(alphaU_current,alphaU_proposal + prop_k_alpha*grad_log_target_5.1_alphaU(data,Yplus,pop_dat,Y_current,alphaU_proposal,alphaR_current,b_current,d_current,mu_U,sd_U),sqrt(2*prop_k_alpha),log=T) -
          dnorm(alphaU_proposal,alphaU_current + prop_k_alpha*grad_log_target_5.1_alphaU(data,Yplus,pop_dat,Y_current,alphaU_current,alphaR_current,b_current,d_current,mu_U,sd_U),sqrt(2*prop_k_alpha),log=T) +
          dnorm(alphaR_current,alphaR_proposal + prop_k_alpha*grad_log_target_5.1_alphaR(data,Yplus,pop_dat,Y_current,alphaU_current,alphaR_proposal,b_current,d_current,mu_R,sd_R),sqrt(2*prop_k_alpha),log=T) -
          dnorm(alphaR_proposal,alphaR_current + prop_k_alpha*grad_log_target_5.1_alphaR(data,Yplus,pop_dat,Y_current,alphaU_current,alphaR_current,b_current,d_current,mu_R,sd_R),sqrt(2*prop_k_alpha),log=T) +
          # MALA proposal dist for b
          sum(vapply(1:n_areas,function(i){
            dnorm(b_current[i],b_proposal[i] + prop_k_b*grad_log_target_5.1_b(data,Yplus,pop_dat,Y_current,i,alphaU_current,alphaR_current,b_proposal,d_current,tau_current),sqrt(2*prop_k_b),log=T) -
              dnorm(b_proposal[i],b_current[i] + prop_k_b*grad_log_target_5.1_b(data,Yplus,pop_dat,Y_current,i,alphaU_current,alphaR_current,b_current,d_current,tau_current),sqrt(2*prop_k_b),log=T)
          },numeric(1))) +
          # random walk proposal dist for d
          #dlnorm(d_current,log(d_proposal),prop_sd_d,log=T) - dlnorm(d_proposal,log(d_current),prop_sd_d,log=T)
          # MALA proposal dist for d
          dlnorm(d_current,log(d_proposal) + prop_k_logd*grad_log_target_5.1_logd(data,Yplus,pop_dat,Y_current,alphaU_current,alphaR_current,b_current,d_proposal,mu_d,sd_d),sqrt(2*prop_k_logd),log=T) -
          dlnorm(d_proposal,log(d_current) + prop_k_logd*grad_log_target_5.1_logd(data,Yplus,pop_dat,Y_current,alphaU_current,alphaR_current,b_current,d_current,mu_d,sd_d),sqrt(2*prop_k_logd),log=T)
        
      })
      
      if(a_prob_reg>runif(1)){
        d_current <- d_proposal
        b_current <- b_proposal
        alphaU_current <- alphaU_proposal
        alphaR_current <- alphaR_proposal
        Nr_sum_current <- Nr_sum_proposal
        Nr_obs_current <- Nr_obs_proposal
        if(i>burnin){acc_reg <- acc_reg + 1}
      }
      
      # accept or reject Ys
      for(k in 1:nrow(data)){
        # determine poisson approximation
        domain <- data$Z[k]:round(0.25*data$N[k])
        f <- log_target_cond_5.1_y(domain,data,Yplus,pop_dat,Y_current,k,Nr_sum_current,Nr_obs_current,d_current)
        del_f <- f-mean(f)
        target_density <- (exp(del_f)/sum(exp(del_f)))
        par_tmp <- optimize(f=function(x) sum((target_density - dpois(domain,x))^2), 
                            lower=max(which.max(target_density)-2,0),upper=which.max(target_density)+2)$minimum

        # draw from poisson
        y_proposal <- rpois(1,par_tmp)

        # accept of reject draw
        if(y_proposal %in% domain)
          a_prob_Y_tmp <- exp(log_target_cond_5.1_y(y_proposal,data,Yplus,pop_dat,Y_current,k,Nr_sum_current,Nr_obs_current,d_current) -
                                log_target_cond_5.1_y(Y_current[k],data,Yplus,pop_dat,Y_current,k,Nr_sum_current,Nr_obs_current,d_current) +
                                dpois(Y_current[k],par_tmp,log=T) - dpois(y_proposal,par_tmp,log=T))
        else
          a_prob_Y_tmp <- 0

        if(a_prob_Y_tmp>runif(1)){
          Y_current[k] <- y_proposal
          if(i>burnin){acc_Y[k] <- acc_Y[k] + 1}
        }
      }
      
      
      # record current state
      tau_postsamp_t[i] <- tau_current
      d_postsamp_t[i] <- d_current
      alphaU_postsamp_t[i] <- alphaU_current
      alphaR_postsamp_t[i] <- alphaR_current
      b_postsamp_t[i,] <- b_current
      Y_postsamp_t[i,] <- Y_current
      
    } 
    Sys.time() - time.start
    
   return(list(alphaU_t = alphaU_postsamp_t[(burnin+1):n_iter], alphaR_t = alphaR_postsamp_t[(burnin+1):n_iter],
                b_t=b_postsamp_t[(burnin+1):n_iter,], 
                tau_t=tau_postsamp_t[(burnin+1):n_iter], d_t=d_postsamp_t[(burnin+1):n_iter], 
                Y_t=Y_postsamp_t[(burnin+1):n_iter,],
                acc_reg_t = acc_reg, acc_Y_t = acc_Y))
  })

  #combine chains
  b_postsamp <- Y_postsamp <- NULL
  for(j in 1:chains){
    b_postsamp <- rbind(b_postsamp,parallel_results[[j]]$b_t)
    Y_postsamp <- rbind(Y_postsamp,parallel_results[[j]]$Y_t)
  }
  alphaU_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$alphaU_t}))
  alphaR_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$alphaR_t}))
  tau_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$tau_t}))
  d_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$d_t}))
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp, alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  acc_reg_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_reg_t})))/((n_iter-burnin)*chains)
  acc_Y_rate <- Reduce('+', lapply(1:chains,function(x){parallel_results[[x]]$acc_Y_t}))/((n_iter-burnin)*chains)

  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              tau=tau_postsamp, d=d_postsamp, Y=Y_postsamp,
              acc_reg_rate = acc_reg_rate, acc_Y_rate = acc_Y_rate))
}

#kernel of loglikelihood with respect to all regression parameters (all except Y)
log_target_5.1_reg_params <- function(params,data_list,Y,tau,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R){
  d <- exp(params[1])
  alphaU <- params[2]
  alphaR <- params[3]
  b <- params[4:length(params)]
  
  Nr_obs <- data_list$obs_dat$N*exp(data_list$obs_dat$U*alphaU + (1-data_list$obs_dat$U)*alphaR + b[data_list$obs_dat$A])
  Nr_sum_obs <- sum(Nr_obs)
  Nr_sum <- sum(data_list$pop_strata_dat$N*exp(data_list$pop_strata_dat$U*alphaU + (1-data_list$pop_strata_dat$U)*alphaR + b[data_list$pop_strata_dat$A]))
  
  term1 <- data_list$Yplus*log(d/(1+d)) - log(1+d)*Nr_sum/d
  term2 <- lgamma(data_list$Yplus - sum(Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - lgamma((1/d)*(Nr_sum-Nr_sum_obs)) 
  term3 <- sum(lgamma(Y+1/d*Nr_obs) - lgamma(1/d*Nr_obs))
  term4 <- log(d) - tau*sum(b^2)/2 - (alphaU-mu_U)^2/(2*sd_U^2) - (alphaR-mu_R)^2/(2*sd_R^2) - (log(d)-mu_d)^2/(2*sd_d^2)
  
  return(term1 + term2 + term3 + term4)
}

log_target_cond_5.1_y <- function(y,data,Yplus,pop_dat,Y,i,Nr_sum,Nr_obs,d){
  term1 <- lgamma(Yplus - sum(Y[-i]) - y + 1/d*(Nr_sum - sum(Nr_obs))) - lgamma(Yplus - sum(Y[-i]) -y + 1)
  term2 <- lgamma(data$N[i]-y+1) - lgamma(y-data$Z[i]+1) - lgamma(data$N[i]-y-data$n[i]+data$Z[i] +1)
  term3 <- lgamma(y + 1/d*Nr_obs[i])
  
  return(term1 + term2 + term3)
}

eps_delta <- function(eps,params,w,L_inv,grad,data_list,Y,tau,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R){
  params_star <- as.vector(params + eps^2/2*t(L_inv)%*%L_inv%*%grad + eps*t(L_inv)%*%w)
  grad_star <- pracma::grad(log_target_5.1_reg_params,x=params_star,data_list=data_list,Y=Y,tau=tau,
                            mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
  r <- L_inv%*%(grad + grad_star)
  out <- -log_target_5.1_reg_params(params,data_list,Y,tau,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R) +
    log_target_5.1_reg_params(params_star,data_list,Y,tau,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R) - 
    eps/2*t(w)%*%r - eps^2/8*t(r)%*%r
  return(as.numeric(out))
}

est_eps0 <- function(params,L_inv,grad,w,eps_max,rho,beta,delta,data_list,Y,tau,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R){
  eps <- eps_max
  del_fn <- eps_delta(eps,params,w,L_inv,grad,data_list,Y,tau,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R)
  while(abs(del_fn)>delta){
    if(abs(del_fn)>beta){
      eps <- rho*eps
    }else{
      eps <- 0.95*(delta/abs(del_fn))^(1/3)*eps
    }
    del_fn <- eps_delta(eps,params,w,L_inv,grad,data_list,Y,tau,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R)
  }
  return(eps)
}

get.GMW.chol <- function(A,u){
  d <- nrow(A)
  L <- diag(1,d)
  D <- diag(diag(A))
  nu <- max(abs(diag(A)))
  zeta <- max(abs(A[lower.tri(A)]))
  phi_sq <- max(nu,zeta/sqrt(d^2-1),u)
  delta <- u*max(nu,zeta,1)
  
  for(j in 1:d){
    if(j>1){
      L[j,1:(j-1)] <- L[j,1:(j-1)]/diag(D)[1:(j-1)]}
    if(j<d){
      L[(j+1):d,j] <- A[(j+1):d,j]}
    if(1<j & j<d){
      if(j==2){
        L[(j+1):d,j] <- L[(j+1):d,j] - L[(j+1):d,1:(j-1)]*L[j,1:(j-1)]
      }else{
        L[(j+1):d,j] <- L[(j+1):d,j] - L[(j+1):d,1:(j-1)]%*%L[j,1:(j-1)]}
    }
    if(j<d){
      theta <- max(abs(L[(j+1):d,j]))
    }else{
      theta <- 0}
    D[j,j] <- max(delta,abs(D[j,j]),theta^2/phi_sq)
    if(j<d){
      diag(D)[(j+1):d] <- diag(D)[(j+1):d] - L[(j+1):d,j]^2/D[j,j]}
  }
  
  return(L = L%*%sqrt(D))
}


tau_init=40
d_init = seq(0.8,1.2,by=0.1) 
alphaU_init = seq(-4,-3,by=0.1)
alphaR_init= seq(-3.5,-2.5,by=0.1)  #initial values
eps_max=1
rho=0.75
beta=10 #params for adaptive step algorithm
a_tau = 0.01
b_tau = 0.01
mu_U = -3.5
sd_U = 9
mu_R = -3
sd_R = 9
mu_d = 0
sd_d = 1#hyperpriors

delta = 0.1
eps=0.7
n_iter = 2000
burnin=n_iter/2
chains=4

#smMALA as described in "Adaptive step size selection for Hessian-based manifold Langevian samplers" Tore Selland Kleppe
postsamp_Alg5.1_smMALA_MCMC <- function(tau_init=40, d_init = seq(0.8,1.2,by=0.1), alphaU_init = seq(-4,-3,by=0.1), alphaR_init= seq(-3.5,-2.5,by=0.1), #initial values
                                        eps=NULL,
                                        eps_max=1, rho=0.75, beta=10, delta=NULL, #params for adaptive step algorithm
                                 a_tau = 0.01,b_tau = 0.01,mu_U = -3.5,sd_U = 9,mu_R = -3,sd_R = 9,mu_d = 0,sd_d = 1,#hyperpriors
                                 data_list, n_iter, burnin=n_iter/2,chains=2){
  
  data <- data_list$obs_dat
  Yplus <- data_list$Yplus
  pop_dat <- data_list$pop_strata_dat
  n_areas <- length(unique(data_list$pop_strata_dat$A))
  
  b_init <- rnorm(n_areas,0,1/sqrt(tau_init))
  
  ## ensure we are starting with a PD negative hessian
  reg_params_init_grid <- as.matrix(expand.grid(log(d_init),alphaU_init,alphaR_init))
  search_order <- sample(1:nrow(reg_params_init_grid),nrow(reg_params_init_grid))
  
  for(j in 1:nrow(reg_params_init_grid)){
    reg_params_init <-  c(reg_params_init_grid[search_order[j],],b_init)
    Nr_sum_init <- sum(pop_dat$N*exp(pop_dat$U*reg_params_init[2] + (1-pop_dat$U)*reg_params_init[3] + reg_params_init[3+pop_dat$A]))
    Nr_obs_init <- data$N*exp(data$U*reg_params_init[2] + (1-data$U)*reg_params_init[3] + reg_params_init[3+data$A])
    Y_init <- rmultinom(1,Yplus,prob=c(Nr_obs_init,Nr_sum_init-sum(Nr_obs_init)))
    Y_init <- Y_init[-length(Y_init)]
    Y_init[Y_init<data$Z] <- data_list$obs_dat$Z[Y_init<data$Z]
    
    hess_init <- pracma::hessian(log_target_5.1_reg_params,x=reg_params_init,data_list=data_list,
                                 Y=Y_init,tau=tau_init,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
    if(is.positive.definite(-hess_init)){
      break
    }
  }

  #initialize data frames
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  reg_params_postsamp_t <- rbind(reg_params_init, matrix(NA,n_iter-1,length(reg_params_init)))
  Y_postsamp_t <- rbind(Y_init, matrix(NA,n_iter-1,nrow(data)))
  acc_reg <- 0
  acc_Y <- rep(0,nrow(data))
  
  # set current states
  tau_current <- tau_init
  d_current <- exp(reg_params_init[1])
  alphaU_current <- reg_params_init[2]
  alphaR_current <- reg_params_init[3]
  b_current <- b_init
  reg_params_current <- reg_params_init
  Y_current <- Y_init
  
   parallel_results <- mclapply(1:chains, function(j){

     #start.time <- Sys.time()
    for(i in 2:n_iter){
      
      # draw new tau from full conditional
      tau_current <- rgamma(1,n_areas/2+a_tau,0.5*sum(b_current^2)+b_tau)
      
      # draw proposal for regression params
      grad_current <- pracma::grad(log_target_5.1_reg_params,x=reg_params_current,data_list=data_list,
                                    Y=Y_current,tau=tau_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)

      hess_current <- pracma::hessian(log_target_5.1_reg_params,x=reg_params_current,data_list=data_list,
                                     Y=Y_current,tau=tau_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      if(is.positive.definite(-hess_current)){
        L_current <- t(chol(-hess_current))
      }else{
        message(paste0('current Hessian is not PD at iteration ',i))
        #L_current_inv <- diag(1,length(reg_params_current))
        L_current <- get.GMW.chol(-hess_current,0.01)
      }
      L_current_inv <- matlib::inv(L_current)
      
      # auxiliary variable for adaptive step size
      #w <- rnorm(length(reg_params_current))
      
      #calculate forward eps
      # eps0_f <- est_eps0(reg_params_current,L_current_inv,grad_current,w,
      #   eps_max,rho,beta,delta,data_list,Y_current,tau_current,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R)
      # eps_f <- min(eps_max,eps0_f)
      eps_f <- eps
      
      reg_params_proposal <- as.vector(reg_params_current + eps_f^2/2*t(L_current_inv)%*%L_current_inv%*%grad_current +
             t(Rfast::rmvnorm(1,rep(0,length(reg_params_current)),eps_f^2*t(L_current_inv)%*%L_current_inv))) #optimize later

      grad_proposal <- pracma::grad(log_target_5.1_reg_params,x=reg_params_proposal,data_list=data_list,
                                     Y=Y_current,tau=tau_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      hess_proposal <- pracma::hessian(log_target_5.1_reg_params,x=reg_params_proposal,data_list=data_list,
                                        Y=Y_current,tau=tau_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)

      if(is.symmetric(hess_proposal)){
        if(is.positive.definite(-hess_proposal)){
          L_proposal <- t(chol(-hess_proposal))
          L_proposal_inv <- matlib::inv(L_proposal)
          
      #calculate backward epsilon
       # eps0_b <- est_eps0(reg_params_proposal,L_proposal_inv,grad_proposal,w,
       #                        eps_max,rho,beta,delta,data_list,Y_current,tau_current,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R)
       #  eps_b <- min(eps_max,eps0_b)
          eps_b <- eps
        
        # accept or reject proposal for (alpha, b, d)
        a_prob_reg <- min(1,exp({
          log_target_5.1_reg_params(reg_params_proposal,data_list=data_list,Y=Y_current,tau=tau_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R) +
           Rfast::dmvnorm(reg_params_current,reg_params_proposal + eps_b^2/2*t(L_proposal_inv)%*%L_proposal_inv%*%grad_proposal,eps_b^2*t(L_proposal_inv)%*%L_proposal_inv, log=T) - #optimize
           log_target_5.1_reg_params(reg_params_current,data_list=data_list,Y=Y_current,tau=tau_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R) -
           Rfast::dmvnorm(reg_params_proposal,reg_params_current + eps_f^2/2*t(L_current_inv)%*%L_current_inv%*%grad_current,eps_f^2*t(L_current_inv)%*%L_current_inv, log=T) #optimize
         }))

        if(a_prob_reg>runif(1)){
          reg_params_current <- reg_params_proposal
          d_current <- exp(reg_params_proposal[1])
          alphaU_current <- reg_params_proposal[2]
          alphaR_current <- reg_params_proposal[3]
          b_current <- reg_params_proposal[4:length(reg_params_proposal)]
          if(i>burnin){acc_reg <- acc_reg + 1}
        }
      
     }else{
        message(paste0('Proposed Hessian not PD at iteration ',i,'\n'))}
      }else{
        message(paste0('Proposed Hessian not symmetric at iteration ',i,'\n'))}
      
      # accept or reject Ys
      Nr_sum_current <- sum(pop_dat$N*exp(pop_dat$U*alphaU_current + (1-pop_dat$U)*alphaR_current + b_current[pop_dat$A]))
      Nr_obs_current <- data$N*exp(data$U*alphaU_current + (1-data$U)*alphaR_current + b_current[data$A])
      
      for(k in 1:nrow(data)){
        # determine poisson approximation
        domain <- data$Z[k]:round(0.5*data$N[k])
        f <- log_target_cond_5.1_y(domain,data,Yplus,pop_dat,Y_current,k,Nr_sum_current,Nr_obs_current,d_current)
        # update support if too wide
        domain <- domain[f!=-Inf]
        f <- f[f!=-Inf]
        del_f <- f-mean(f)
        target_density <- (exp(del_f)/sum(exp(del_f)))
        par_tmp <- optimize(f=function(x) sum((target_density - dpois(domain,x))^2),
                              lower=max(domain[which.max(f)]-2,0),upper=domain[which.max(f)]+2)$minimum
        
        # draw from poisson
        y_proposal <- rpois(1,par_tmp)

        # accept of reject draw
        if(y_proposal %in% domain)
          a_prob_Y_tmp <- exp(log_target_cond_5.1_y(y_proposal,data,Yplus,pop_dat,Y_current,k,Nr_sum_current,Nr_obs_current,d_current) -
                                log_target_cond_5.1_y(Y_current[k],data,Yplus,pop_dat,Y_current,k,Nr_sum_current,Nr_obs_current,d_current) +
                                dpois(Y_current[k],par_tmp,log=T) - dpois(y_proposal,par_tmp,log=T))
        else
          a_prob_Y_tmp <- 0

        if(a_prob_Y_tmp>runif(1)){
          Y_current[k] <- y_proposal
          if(i>burnin){acc_Y[k] <- acc_Y[k] + 1}
        }
      }
      
      # record current state
      tau_postsamp_t[i] <- tau_current
      if(is.symmetric(hess_proposal)){
        if(is.positive.definite(-hess_proposal)){
          reg_params_postsamp_t[i,] <- reg_params_current
        }
      }
      Y_postsamp_t[i,] <- Y_current
    } 
     #Sys.time() - start.time
    
  results <- list(reg_params_t=reg_params_postsamp_t[(burnin+1):n_iter,], 
                  tau_t=tau_postsamp_t[(burnin+1):n_iter], 
                  Y_t=Y_postsamp_t[(burnin+1):n_iter,],
                  acc_reg_t = acc_reg, acc_Y_t = acc_Y)
  
    return(results)
  },mc.cores = chains)
  
  #combine chains
  b_postsamp <- reg_params_postsamp <- Y_postsamp <- NULL
  for(j in 1:chains){
    b_postsamp <- rbind(b_postsamp,parallel_results[[j]]$reg_params_t[,4:length(reg_params_init)])
    reg_params_postsamp <- rbind(reg_params_postsamp,parallel_results[[j]]$reg_params_t)
    Y_postsamp <- rbind(Y_postsamp,parallel_results[[j]]$Y_t)
  }
  tau_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$tau_t}))
  d_postsamp <- unlist(lapply(1:chains,function(x){exp(parallel_results[[x]]$reg_params_t[,1])}))
  alphaU_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,2]}))
  alphaR_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,3]}))
  names(tau_postsamp) <- names(d_postsamp) <- names(alphaU_postsamp) <- names(alphaR_postsamp) <- NULL
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp,alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  
  acc_reg_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_reg_t})))/((n_iter-burnin)*chains)
  acc_Y_rate <- Reduce('+', lapply(1:chains,function(x){parallel_results[[x]]$acc_Y_t}))/((n_iter-burnin)*chains)
  
  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              tau=tau_postsamp, 
              reg_params=reg_params_postsamp,
              d=d_postsamp, 
              Y=Y_postsamp,
              acc_reg_rate = acc_reg_rate, acc_Y_rate = acc_Y_rate))
}



#######################################################################
#######################################################################
# Simulations for Alg 5.1-5.3 -----------

## set populations and simulation parameters
{
  ## this section is calibrated to be relatively similar to Sierra Leone 2019 survey, in terms of sample sizes
  n_areas <- 20 #numbers of areas
  n_clusters_urban <- rnegbin(n_areas,200,5) #number of urban clusters in each area
  n_clusters_rural <- rnegbin(n_areas,300,5) #number of rural clusters in each area
  n_births_urban <- 30 #average number of births per urban cluster (0.1*number of households*3) <- assuming 3 year period
  n_births_rural <-  40 #average number of births per rural cluster (0.15*number of households*3) <- assuming 3 year period
  n_clusters_urban_samp <-  round(0.075*n_clusters_urban) #number of urban clusters sampled
  n_clusters_rural_samp <-  round(0.05*n_clusters_rural) #number of rural clusters sampled
  n_births_urban_samp <- 10 # average number of births sampled from each urban cluster (0.1*30*3) -- change to be a little smaller later
  n_births_rural_samp <- 15 # average number of births sampled from each rural cluster -- change to be a little smaller later
  
  # intercept (fixed effect) 
  alphaU <- rnorm(1,-3.5,0.25)
  alphaR <- rnorm(1,-3,0.25)
  
  # draw hyperparameters
  tau <- rgamma(1,10,.25)
  
  # draw overdispersion parameter
  d <- rlnorm(1,0,0.25)
  
  # draw random effects
  b <- as.vector(Rfast::rmvnorm(1,rep(0,n_areas),diag(1/tau,n_areas)))
}

## generate data 
{
  all_dat <- data.frame(
    # (true) number of births in each cluster
    N = round(c(rnorm(sum(n_clusters_urban),n_births_urban,n_births_urban/10),rnorm(sum(n_clusters_rural),n_births_rural,n_births_rural/10))),
    # urban or rural strata (U=1 if urban, otherwise 0)
    U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
    # admin area of each cluster
    A = c(unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_rural[i])}))))
  all_dat$cluster <- 1:nrow(all_dat)
  all_dat$Y <- sapply(1:nrow(all_dat),function(i){rnbinom(1,
                                                          size = 1/d*all_dat$N[i]*exp(alphaU*all_dat$U[i] + alphaR*(1-all_dat$U[i]) + b[all_dat$A[i]] + rnorm(1,0,0.05)),
                                                          prob=1/(1+d))})
  
  ### sample clusters
  obs_dat <- NULL
  for(area in 1:n_areas){
    # randomly select clusters from urban strata based on size of cluster
    obs_dat <- rbind(obs_dat,all_dat[all_dat$A==area & all_dat$U==1,][sample(1:n_clusters_urban[area],n_clusters_urban_samp[area],prob = all_dat[all_dat$U==1 & all_dat$A==area,]$N),],
                     # randomly select cluster from rural strata  based on size of cluster
                     all_dat[all_dat$A==area & all_dat$U==0,][sample(1:n_clusters_rural[area],n_clusters_rural_samp[area],prob = all_dat[all_dat$U==0& all_dat$A==area,]$N),])
  }
  
  # sample births
  obs_dat$n <- obs_dat$Z <- NA
  obs_dat[obs_dat$U==1,]$n <- round(min(rnorm(sum(obs_dat$U),n_births_urban_samp,n_births_urban_samp/10),obs_dat[obs_dat$U==1,]$N))
  obs_dat[obs_dat$U==0,]$n <- round(min(rnorm(sum(1-obs_dat$U),n_births_rural_samp,n_births_rural_samp/10),obs_dat[obs_dat$U==0,]$N))
  for(i in 1:nrow(obs_dat)){
    obs_dat$Z[i] <- sum(sample(c(rep(1,obs_dat$Y[i]),rep(0,(obs_dat$N[i]-obs_dat$Y[i]))), obs_dat$n[i]))
  }
  
  
  data_list <- list(obs_dat = obs_dat, #observed data
                    Yplus = sum(all_dat$Y), #ALL deaths
                    pop_strata_dat = all_dat %>% group_by(A,U) %>% summarise(N=sum(N))) # number of births per strata (area x urban)
}

## run model once

start.time <- Sys.time()
postsamp_5.1 <- postsamp_Alg5.1_smMALA_MCMC(eps=0.7, data_list = data_list, 
                                            n_iter=4000,chains = 4)
Sys.time() - start.time

postsamp_5.1$acc_reg_rate
quantile(postsamp_5.1$acc_Y_rate,seq(0,1,0.1))
sum(is.na(postsamp_5.1$reg_params[,1]))

pdf("Handcoded MCMC Simulations/Alg 5.1 smMALA 231218, 4 chains of 5k iterations posterior dist.pdf")
{
  plot(postsamp_5.1$reg_params[1:1000%%10==0,1],type='l')
  abline(h=log(d),col='red')
  plot(postsamp_5.1$reg_params[1:1000%%10==0,2],type='l')
  abline(h=alphaU,col='red')
  plot(postsamp_5.1$reg_params[1:1000%%10==0,3],type='l')
  abline(h=alphaR,col='red')

  par(mfrow=c(3,2))
  for(k in 4:(n_areas+3)){
    plot(postsamp_5.1$reg_params[1:1000%%10==0,k],type='l',main=k)
    abline(h=b[k-3],col='red')
  }

  par(mfrow=c(3,2))
  for(k in sample(1:nrow(data_list$obs_dat),30,replace = F)){
   plot(postsamp_5.1$Y[1:1000%%10==0,k],type='l',main=k)
    abline(h=obs_dat$Y[k],col='red')
  }

  hist(log(postsamp_5.1$tau))
  abline(v=log(tau),col='red')
  hist(log(postsamp_5.1$d))
  abline(v=log(d),col='red')
  hist(postsamp_5.1$alphaU)
  abline(v=alphaU,col='red')
  hist(postsamp_5.1$alphaR)
  abline(v=alphaR,col='red')

  results_r <- gather(data.frame(postsamp_5.1$r),'strata','r')
  results_r$strata <- as.numeric(str_remove(results_r$strata,'V'))
  true_r <- data.frame(strata = 1:(2*n_areas),r=exp(c(alphaU+b,alphaR+b)))
  true_r <- true_r[order(true_r$r),]
  true_r$r_rank <- 1:nrow(true_r)
  results_r <- left_join(results_r,true_r,by='strata')
  g <- results_r %>% filter(strata %in% 1:20) %>% ggplot() + geom_boxplot(aes(y=r.x)) +
    geom_hline(aes(yintercept = r.y),col='red') + facet_grid(~r_rank) + ggtitle('Posterior for NMR (for each urban strata)') 
  print(g)
  g <- results_r %>% filter(strata %in% 21:40) %>% ggplot() + geom_boxplot(aes(y=r.x)) +
    geom_hline(aes(yintercept = r.y),col='red') + facet_grid(~r_rank) + ggtitle('Posterior for NMR (for each rural strata)') 
  print(g)
  
  # results_Y <- gather(data.frame(postsamp_5.1$Y),'cluster','Y')
  # results_Y$cluster <- as.numeric(str_remove(results_Y$cluster,'X'))
  # true_Y <- data.frame(cluster=1:nrow(data_list$obs_dat), Y=data_list$obs_dat$Y)
  # true_Y <-true_Y[order(true_Y$Y),]
  # true_Y$rank <- 1:nrow(true_Y)
  # clusters_plot <- seq(1,nrow(true_Y),50)
  # results_Y_plot <- results_Y[results_Y$cluster %in% clusters_plot,]
  # results_Y_plot <- left_join(results_Y_plot,true_Y,by='cluster')
  # g <- results_Y_plot %>% ggplot() + geom_boxplot(aes(y=Y.x)) +
  #   geom_hline(aes(yintercept = Y.y),col='red') + facet_grid(~rank) #+ ggtitle('XX births per area') +
  # print(g)
}
dev.off()

## run model multiple times
n_runs <- 50
tau_postmeans <- acc_reg <- na_count <- rep(NA,n_runs)
reg_params_postmeans <- matrix(NA,n_runs,23)
r_postmeans <- matrix(NA,n_runs,2*n_areas)
Y_postmedians <- matrix(NA,n_runs,nrow(obs_dat))
acc_Y <- list()
start.time <- Sys.time()
for(q in 1:25){
  message(paste0('Starting simulation ',q,'\n'))
  ## generate data 
  {
    all_dat <- data.frame(
      # (true) number of births in each cluster
      N = round(c(rnorm(sum(n_clusters_urban),n_births_urban,n_births_urban/10),rnorm(sum(n_clusters_rural),n_births_rural,n_births_rural/10))),
      # urban or rural strata (U=1 if urban, otherwise 0)
      U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
      # admin area of each cluster
      A = c(unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_rural[i])}))))
    all_dat$cluster <- 1:nrow(all_dat)
    all_dat$Y <- sapply(1:nrow(all_dat),function(i){rnbinom(1,
                                                            size = 1/d*all_dat$N[i]*exp(alphaU*all_dat$U[i] + alphaR*(1-all_dat$U[i]) + b[all_dat$A[i]] + rnorm(1,0,0.05)),
                                                            prob=1/(1+d))})
    
    ### sample clusters
    obs_dat <- NULL
    for(area in 1:n_areas){
      # randomly select clusters from urban strata based on size of cluster
      obs_dat <- rbind(obs_dat,all_dat[all_dat$A==area & all_dat$U==1,][sample(1:n_clusters_urban[area],n_clusters_urban_samp[area],prob = all_dat[all_dat$U==1 & all_dat$A==area,]$N),],
                       # randomly select cluster from rural strata  based on size of cluster
                       all_dat[all_dat$A==area & all_dat$U==0,][sample(1:n_clusters_rural[area],n_clusters_rural_samp[area],prob = all_dat[all_dat$U==0& all_dat$A==area,]$N),])
    }
    
    # sample births
    obs_dat$n <- obs_dat$Z <- NA
    obs_dat[obs_dat$U==1,]$n <- round(min(rnorm(sum(obs_dat$U),n_births_urban_samp,n_births_urban_samp/10),obs_dat[obs_dat$U==1,]$N))
    obs_dat[obs_dat$U==0,]$n <- round(min(rnorm(sum(1-obs_dat$U),n_births_rural_samp,n_births_rural_samp/10),obs_dat[obs_dat$U==0,]$N))
    for(i in 1:nrow(obs_dat)){
      obs_dat$Z[i] <- sum(sample(c(rep(1,obs_dat$Y[i]),rep(0,(obs_dat$N[i]-obs_dat$Y[i]))), obs_dat$n[i]))
    }
    
    
    data_list <- list(obs_dat = obs_dat, #observed data
                      Yplus = sum(all_dat$Y), #ALL deaths
                      pop_strata_dat = all_dat %>% group_by(A,U) %>% summarise(N=sum(N))) # number of births per strata (area x urban)
  }
  
  start.time <- Sys.time()
  postsamp_5.1 <- postsamp_Alg5.1_smMALA_MCMC(eps=0.7,
                                              data_list = data_list, 
                                              n_iter=4000,chains = 4)
  Sys.time() - start.time
  
  #make sure things are going ok 
  message(paste0('Acceptance rate for regression parameters is ',postsamp_5.1$acc_reg_rate,'\n'))
  
  tau_postmeans[q] <- mean(postsamp_5.1$tau,na.rm=T)
  reg_params_postmeans[q,] <- colMeans(postsamp_5.1$reg_params,na.rm=T)
  r_postmeans[q,] <- colMeans(postsamp_5.1$r,na.rm=T)
  Y_postmedians[q,] <- colMedians(postsamp_5.1$Y,na.rm=T)
  acc_reg[q] <- postsamp_5.1$acc_reg_rate
  acc_Y[[q]] <- postsamp_5.1$acc_Y_rate
  na_count[q] <- sum(is.na(postsamp_5.1$reg_params[,1]))
}
Sys.time() - start.time

results=list(tau_postmeans=tau_postmeans,reg_params_postmeans=reg_params_postmeans,r_postmeans=r_postmeans,
             acc_reg=acc_reg,d=d,alphaU=alphaU,alphaR=alphaR,b=b,tau=tau)

good_runs <- which(na_count<1000)

save(file='Handcoded MCMC Simulations/Alg 5.1 smMALA 231216, 50 datasets 2k iterations.rda',results)

pdf("Handcoded MCMC Simulations/Alg 5.1 smMALA 231216, 50 datasets 2k iterations.pdf")
{
hist(log(tau_postmeans[good_runs]))
abline(v=log(tau),col='red')
hist(reg_params_postmeans[good_runs,1],main='log(d) Posterior means')
abline(v=log(d),col='red')
hist(reg_params_postmeans[good_runs,2],main='alphaU Posterior means')
abline(v=alphaU,col='red')
hist(reg_params_postmeans[good_runs,3],main='alphaR Posterior means')
abline(v=alphaR,col='red')

results_r <- gather(data.frame(r_postmeans[good_runs,]),'strata','r')
results_r$strata <- as.numeric(str_remove(results_r$strata,'X'))
true_r <- data.frame(strata = 1:(2*n_areas),r=exp(c(alphaU+b,alphaR+b)))
true_r <- true_r[order(true_r$r),]
true_r$r_rank <- 1:nrow(true_r)
results_r <- left_join(results_r,true_r,by='strata')
g <- results_r %>% filter(strata %in% 1:20) %>% ggplot() + geom_boxplot(aes(y=r.x)) +
  geom_hline(aes(yintercept = r.y),col='red') + facet_grid(~r_rank) + ggtitle('Posterior mean for NMR (for each urban strata)') 
print(g)
g <- results_r %>% filter(strata %in% 21:40) %>% ggplot() + geom_boxplot(aes(y=r.x)) +
  geom_hline(aes(yintercept = r.y),col='red') + facet_grid(~r_rank) + ggtitle('Posterior mean for NMR (for each rural strata)') 
print(g)
}
dev.off()

#######################################################################
#######################################################################
###### DEPRECATED CODE ##########

## Version of Alg 2.3 where each parameter is updated separately
postsamp_artauphi_Y_MCMC_DEP <- function(tau_init,phi_init,alpha_init, #initial values
                                         prop_sd_phi,
                                         data,Q_inv,n_iter,burnin=n_iter/2,chains=2){
  n_areas <- length(unique(data$A))
  #collapse data into areas only
  newdat <- data.frame(data %>% group_by(A) %>% summarise(Y=sum(Y),N=sum(N)))
  #hyperpriors
  a_tau <- b_tau <- 0.01
  a_phi <- 2
  b_phi <- 3
  mu_a <- -6
  sd_a <- 3
  tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
  phi_postsamp_t <- c(phi_init, rep(NA,n_iter-1))
  alpha_postsamp_t <- c(alpha_init, rep(NA,n_iter-1))
  acc_phi <- acc_eta <- 0
  
  eta_init <- Rfast::rmvnorm(1,rep(alpha_init,n_areas),1/tau_init*(W(phi_init,Q_inv))) 
  eta_postsamp_t <- rbind(eta_init, matrix(NA,n_iter-1,n_areas))
  tau_postsamp <- phi_postsamp <- alpha_postsamp <- b_postsamp <- eta_postsamp <- r_postsamp <- NULL
  for(chain in 1:chains){
    for(i in 2:n_iter){
      tau_current <- tau_postsamp_t[i-1]
      phi_current <- phi_postsamp_t[i-1]
      alpha_current <- alpha_postsamp_t[i-1]
      eta_current <- eta_postsamp_t[(i-1),]
      
      # Draw new tau from full conditional
      tau_current <- rgamma(1,a_tau + n_areas/2,
                            b_tau + 0.5*t(eta_current-alpha_current)%*%solve(W(phi_current,Q_inv))%*%(eta_current-alpha_current))
      
      # Draw proposal for phi
      phi_proposal <- expit(rnorm(1,logitlink(phi_current),prop_sd_phi)) 
      
      aprob_phi <- exp(-0.5*as.numeric(determinant(W(phi_proposal,Q_inv))$modulus) - #normalizing constant (maybe find faster way later)
                         (-0.5*as.numeric(determinant(W(phi_current,Q_inv))$modulus)) +
                         -0.5*tau_current*t(eta_current-alpha_current)%*%solve(W(phi_proposal,Q_inv))%*%(eta_current-alpha_current) - #likelihood eta-a|phi
                         (-0.5*tau_current*t(eta_current-alpha_current)%*%solve(W(phi_current,Q_inv))%*%(eta_current-alpha_current)) +
                         (a_phi-1)*(log(phi_proposal/phi_current)) + (b_phi-1)*(log((1-phi_proposal)/(1-phi_current))) + #prior
                         # pc_log_pdf_phi(phi_proposal,lambda_phi,Q_inv,gamma_til_phi) - pc_log_pdf_phi(phi_current,lambda_phi,Q_inv,gamma_til_phi) + #PC prior
                         (-log(phi_current)-log(1-phi_current) +log(phi_proposal) + log(1-phi_proposal)))   # proposal distribution
      if(runif(1) < aprob_phi){
        phi_current <- phi_proposal
        if(i>burnin){
          acc_phi <- acc_phi+1 
        }
      }
      
      # Draw new alpha from full conditional
      alpha_current <- rnorm(1, (mu_a + sd_a^2*tau_current*sum(t(eta_current)%*%solve(W(phi_current,Q_inv))))/(1+sd_a^2*tau_current*sum(solve(W(phi_current,Q_inv)))),
                             1/sqrt(tau_current*sum(solve(W(phi_current,Q_inv))) + 1/sd_a^2))
      
      # Draw proposal for eta for Taylor approximation of full conditional
      eta_proposal <- as.vector(Rfast::rmvnorm(1,prop_mean_bym2_eta(eta_current,tau_current,phi_current,alpha_current,newdat,Q_inv),
                                               prop_sigma_bym2_eta(eta_current,tau_current,phi_current,newdat,Q_inv)))
      
      Var_b <- (1/tau_current)*W(phi_current,Q_inv)
      aprob_eta <- exp(sum(sapply(1:n_areas,function(area){
        dpois(newdat$Y[area],newdat$N[area]*exp(eta_proposal[area]),log=T) - #likelihood
          dpois(newdat$Y[area],newdat$N[area]*exp(eta_current[area]),log=T)})) +
          Rfast::dmvnorm(eta_proposal,rep(alpha_current,n_areas),Var_b, logged = T) - #prior
          Rfast::dmvnorm(eta_current,rep(alpha_current,n_areas),Var_b, logged = T) +
          Rfast::dmvnorm(eta_current,prop_mean_bym2_eta(eta_proposal,tau_current,phi_current,alpha_current,newdat,Q_inv), # target distribution
                         prop_sigma_bym2_eta(eta_proposal,tau_current,phi_current,newdat,Q_inv),logged = T) - 
          Rfast::dmvnorm(eta_proposal,prop_mean_bym2_eta(eta_current,tau_current,phi_current,alpha_current,newdat,Q_inv), 
                         prop_sigma_bym2_eta(eta_current,tau_current,phi_current,newdat,Q_inv),logged = T))
      
      if(runif(1) < aprob_eta){
        eta_current <- eta_proposal
        if(i>burnin){
          acc_eta <- acc_eta+1
        }
      }
      
      # record state
      tau_postsamp_t[i] <- tau_current
      phi_postsamp_t[i] <- phi_current
      alpha_postsamp_t[i] <- alpha_current
      eta_postsamp_t[i,] <- eta_current
      
    } 
    #get rid of burn-in and add to other chains
    eta_postsamp <- rbind(eta_postsamp,eta_postsamp_t[(burnin+1):n_iter,])
    alpha_postsamp <- c(alpha_postsamp,alpha_postsamp_t[(burnin+1):n_iter])
    tau_postsamp <- c(tau_postsamp,tau_postsamp_t[(burnin+1):n_iter])
    phi_postsamp <- c(phi_postsamp,phi_postsamp_t[(burnin+1):n_iter])
  }
  
  b_postsamp <- eta_postsamp-alpha_postsamp
  r_postsamp <- exp(eta_postsamp)
  acc_rates <- matrix(c(acc_phi/((n_iter-burnin)*chains),acc_eta/((n_iter-burnin)*chains)),nrow=1)
  colnames(acc_rates) <- c('phi','eta') 
  
  return(list(alpha = alpha_postsamp,b=b_postsamp, eta=eta_postsamp,r=r_postsamp, tau=tau_postsamp,phi=phi_postsamp,acc_rates=acc_rates))
}

#######################################################################
#######################################################################