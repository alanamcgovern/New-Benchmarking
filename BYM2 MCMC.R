
## helper functions -----
# returns (1-phi)I + phiQ^(-1) (matrix component of BYM2 variance)
W <- function(phi,Q_inv){ 
  return(diag(1-phi,nrow(Q_inv)) + phi*Q_inv)
} 

#mean of prop dist for eta in BYM2 model
prop_mean_bym2_eta <- function(eta_star,tau,phi,alpha,data,Q_inv,Y=NULL){ 
  if(is.null(Y)){
    Y <- data$Y}
  out <- (Y + alpha*tau*colSums(solve(W(phi,Q_inv))) - data$N*exp(eta_star)*(1-eta_star))
  return(out%*%prop_sigma_bym2_eta(eta_star,tau,phi,data,Q_inv))
} 

#var of prop dist for eta in BYM2 model
prop_sigma_bym2_eta <- function(eta_star,tau,phi,data,Q_inv){ 
  out <- tau*solve(W(phi,Q_inv)) + diag(data$N*exp(eta_star))
  return(solve(out))} 

# generate spatial data according to BYM2
simdat_bym2sp <- function(n_areas, # number of areas (integer)
                          n_clusters=NULL, # number of clusters in each area (area-length vector)
                          n_births=NULL, # average number of births in each cluster (integer)
                          n_clusters_urban=NULL, # number of urban clusters in each area (area-length vector)
                          n_clusters_rural=NULL, # number of rural clusters in each area (area-length vector)
                          n_births_urban=NULL, # average number of births in each urban cluster (integer)
                          n_births_rural=NULL, # average number of births in each rural cluster (integer)
                          a=NULL, Var_b,b=NULL,
                          alphaU=NULL, alphaR=NULL, e=NULL){
  
  if(is.null(b)==T){
    #random effects
    b <- Rfast::rmvnorm(1,rep(0,n_areas),Var_b)
  }
  
  #generate without stratification
  if(!is.null(n_births)){
    #observed data
    dat <- data.frame(
      # (true) number of births in each cluster
      N = rpois(sum(n_clusters),n_births),
      # admin area of each cluster
      A = c(sapply(1:n_areas,function(i){rep(i,n_clusters[i])})))
    dat$Y <- sapply(1:nrow(dat),function(i){rpois(1,dat$N[i]*exp(a + b[dat$A[i]]+rnorm(1,0,0.05)))})
    dat$cluster <- 1:nrow(dat)
    return(list(dat=dat,b=b, Var_b=Var_b))
    
    #generate with stratification
  }else{
    dat <- data.frame(
      # (true) number of births in each cluster
      N = c(rpois(sum(n_clusters_urban),n_births_urban),rpois(sum(n_clusters_rural),n_births_rural)),
      # urban or rural strata (U=1 if urban, otherwise 0)
      U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
      # admin area of each cluster
      A = c(unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_areas,function(i){rep(i,n_clusters_rural[i])}))))
    dat$Y <- sapply(1:nrow(dat),function(i){rpois(1,dat$N[i]*exp(alphaU*dat$U[i] + alphaR*(1-dat$U[i]) + b[dat$A[i]] + e[i] +rnorm(1,0,0.05)))})
    dat$cluster <- 1:nrow(dat)
    return(list(dat=dat,b=b, Var_b=Var_b, e=e))
  }
}


## MCMC algorithm function ------
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

## define ground truth -----

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
  
  # intercept (fixed effect) 
  alpha <- rnorm(1,-6,0.25)
  
  # draw hyperparameters
  tau <- rgamma(1,10,.25)
  phi <- rbeta(1,2,3)

## simulate data ------
  Var_b <- (1/tau)*W(phi,Q_scaled_inv)
  b <- Rfast::rmvnorm(1,rep(0,n_areas),Var_b)
  dat2 <- simdat_bym2sp(n_areas=n_areas,n_clusters=n_clusters,n_births=n_births,a=alpha,Var_b=Var_b,b=b)$dat
  

  
## run MCMC for one dataset ----
postsamp_2.3 <- postsamp_artauphi_Y_MCMC(tau_init = 50,phi_init = 0.15,alpha_init = -6,
                                         prop_sd_phi=.25,prop_sd_tau = .25,
                                         data = dat2,Q_inv = Q_scaled_inv,10000,chains=2)

## Some diagnostic plots ------
  plot(log(postsamp_2.3$tau),type='l')
  abline(h=log(tau),col='red')
  plot(logitlink(postsamp_2.3$phi),type='l')
  abline(h=logitlink(phi),col='red')
  plot(postsamp_2.3$alpha,type='l')
  abline(h=alpha,col='red')
  plot(postsamp_2.3$b[,7],type='l')
  abline(h=b[7],col='red')
  plot(postsamp_2.3$b[,15],type='l')
  abline(h=b[15],col='red')
  plot(postsamp_2.3$b[,28],type='l')
  abline(h=b[28],col='red')
  plot(postsamp_2.3$b[,42],type='l')
  abline(h=b[42],col='red')
  
  hist(log(postsamp_2.3$tau))
  abline(v=log(tau),col='red')
  hist(logitlink(postsamp_2.3$phi))
  abline(v=logitlink(phi),col='red')
  hist(postsamp_2.3$alpha)
  abline(v=alpha,col='red')
  