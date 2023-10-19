
library(locfit)
library(LaplacesDemon)
library(Rfast)
library(MASS)
library(INLA)
library(SUMMER)
library(VGAM)

# Simulation 1: verify relationship between Zi and sum of Zis (TRUE counts of deaths (unobservable in practice)) -------
n_areas <- 5
nsim <- 1*10^3

#region specific intercepts, change to BYM2 later
b <- rnorm(n_areas,0,0.5)
#probability of death in each region

# assuming we can know the number of births (worldpop)
N <- rpois(n_areas,10000)

z_area <- z_area_plus <- matrix(NA,nsim,n_areas)
# z's for each area distributed poisson
for(j in 1:n_areas){
  z_area[,j] <- rpois(nsim,N[j]*p[j])
}

# z's for each area distributed multinomial depending on sum of z's
z_plus <- rpois(nsim,sum(N*p))
multinom_prob <- N*p/(sum(N*p))
for(i in 1:nsim){
  z_area_plus[i,] <- rmultinom(1,z_plus[i],multinom_prob)
  
}

# compare distributions
plot(density(z_area[,1]),col=2,lwd=2,xlim=c(0,max(z_area)),ylim=c(0,0.05))
lines(density(z_area_plus[,1]),col=2,lwd=2,lty=2)
for(area in 2:n_areas){
  lines(density(z_area[,area]),col=area+1,lwd=2)
  lines(density(z_area_plus[,area]),col=area+1,lwd=2,lty=2)
}

# Simulate data for MCMC (fixed spatial effect)  -------

#number of admin1 areas
n_areas <- 5
#total number of clusters in each area
n_clusters <- rep(15,n_areas)
#ground truth
b <- rnorm(n_areas,-2,0.25)

#observed data
dat <- data.frame(
  # (true) number of births in each cluster
  N = rpois(sum(n_clusters),500),
  # admin area of each cluster
  A = c(sapply(1:n_areas,function(i){rep(i,n_clusters[i])})))
dat$Y <- sapply(1:nrow(dat),function(i){rpois(1,dat$N[i]*expit(b[dat$A[i]]+rnorm(1,0,0.05)))})
dat$cluster <- 1:nrow(dat)

# MCMC to get P(b)|Ys (fixed spatial effect)  -------

postsamp_b_Y_MCMC <- function(b_init,prop_sd,data,n_iter){
  log_target_dist <- function(b,data){
    #likelihood
    log(dpois(sum(data$Y),sum(data$N)*expit(b))) + 
      #prior
      log(dnorm(b,-2,10))
  }
  n_areas <- length(b_init)
  b_postsamp <- rbind(b_init, matrix(NA,n_iter-1,n_areas))
  for(i in 2:n_iter){
    for(j in 1:n_areas){
      b_current <- b_postsamp[(i-1),j]
      b_proposed <- rnorm(1,b_current,prop_sd^2)
      a_prob <- exp(log_target_dist(b_proposed,data[data$A==j,]) - log_target_dist(b_current,data[data$A==j,]))
      u <- runif(1)
      if(u < a_prob){
        b_postsamp[i,j] <- b_proposed
      }else if(u>= a_prob){
        b_postsamp[i,j] <- b_current
      }
    }
  }  
  #get rid of burn-in
  b_postsamp <- b_postsamp[1001:n_iter,]
  return(list(b=b_postsamp,lambda=expit(b_postsamp)))
}

postsamp_p <- postsamp_b_Y_MCMC(rep(-1.5,n_areas),0.15,dat,10000)

plot(postsamp_p$b[,2],type='l')
plot(postsamp_p$lambda[,2],type='l')
plot(postsamp_p$lambda[,5],type='l')
hist(postsamp_p$lambda[,2])
  
  
# MCMC to get P(weighted sum of theta, q|Ys,Y+)  w priors in terms of b (fixed spatial effect) -------  
postsamp_thetaq_yyplus_MCMC <- function(b_init,prop_sd,data,n_iter){
  log_target_dist <- function(b,data){
                #add up # of deaths across clusters for each admin area
    dmultinom((data%>% group_by(A) %>% summarise(sum(Y)))$`sum(Y)`,
              #Nithetai/sum(Nithetai)
              prob=sapply(1:n_areas,function(j){sum(data[data$A==j,]$N)*expit(b[j])/sum(sapply(1:n_areas,function(i){sum(data[data$A==i,]$N)*expit(b[i])}))})
              ,log=T) +
    dpois(sum(data$Y),sum(sapply(1:n_areas,function(i){sum(data[data$A==i,]$N)*expit(b[i])})),log=T) +
      sum(dnorm(b,-2,10,log = T))
  }
  n_areas <- length(b_init)
  b_postsamp <- rbind(b_init, matrix(NA,n_iter-1,n_areas))
  for(i in 2:n_iter){
    b_current <- b_postsamp[(i-1),]
    for(j in 1:n_areas){
      b_proposed <- b_current
      b_proposed[j] <- rnorm(1,b_current[j],prop_sd^2)
      a_prob <- exp(log_target_dist(b_proposed,data) - log_target_dist(b_current,data))
      u <- runif(1)
      if(u < a_prob){
        b_current[j] <- b_proposed[j]
      }
    }
    b_postsamp[i,] <- b_current
  }
  #get rid of burn-in
  b_postsamp <- b_postsamp[1001:n_iter,]
  
  #get sum of Ni*thetai
  thetaw_postsamp <- sapply(1:nrow(b_postsamp),function(j){sum(sapply(1:n_areas,function(i){sum(data[data$A==i,]$N)*expit(b_postsamp[j,i])}))})
  q_postsamp <-  (sapply(1:n_areas,function(i){sum(data[data$A==i,]$N)})*(expit(b_postsamp)))/thetaw_postsamp
  return(list(b=b_postsamp,
         thetaw=thetaw_postsamp,
         q=q_postsamp))
}
mcmc2 <- postsamp_thetaq_yyplus_MCMC(rep(-1.5,n_areas),0.15,dat,2000)
plot(mcmc2$b[,1],type='l')
plot(mcmc2$b[,3],type='l')
plot(mcmc2$b[,5],type='l')
hist(mcmc2$thetaw)
hist(mcmc2$q[,1])
hist(mcmc2$q[,3])
hist(mcmc2$q[,5])

# Simulate data for MCMC (BYM2 spatial effect) -------  

#GROUND TRUTH
#number of admin1 areas
n_areas <- 5
#total number of clusters in each area
n_clusters <- rep(15,n_areas)
#intercept (fixed effect)
a <- rnorm(1,-5,0.25)
#hyperprior specification for BYM2 precision (tau_b)
pc_alpha_tau <- 0.01
pc_u_tau <- 1
#hyperprior specifications for BYM2 spatial variation (phi)
pc_alpha_phi <- 2/3
pc_u_phi <- 0.5
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

#function to reparamterized covariance matrix
get_Varw <- function(n_areas,Q,tau_b,phi){
  eigQ <- eigen(Q)
  Q_scaled <- inla.scale.model(Q, constr=list(A=t(eigQ$vectors[,eigQ$values<1e-10]), e=rep(0,sum(eigQ$values<1e-10))))
  Prec_w <- rbind(cbind(tau_b/(1-phi)*diag(1,n_areas), -sqrt(phi*tau_b)/(1-phi)*diag(1,n_areas)),
                  cbind(-sqrt(phi*tau_b)/(1-phi)*diag(1,n_areas), Q_scaled + phi/(1-phi)*diag(1,n_areas)))
  Var_w <- ginv(as.matrix(Prec_w))
  Var_w <- nearPD(Var_w,ensureSymmetry = T)$mat
  return(Var_w)
}

#simulate data
simdat <- function(n_areas,n_clusters,Q_scaled,pc_theta,a){
  #draw spatial precision from hyperprior
  tau_b <- INLA:::inla.pc.rprec(1,pc_u_tau,pc_alpha_tau)
  
  #draw spatial variation from uniform hyperprior
  phi <- runif(1)
  
  #reparametrized matrix for BYM2 parameters
  Var_w <- get_Varw(n_areas,Q,tau_b,phi)
  
  #random effects
  w_draw <- Rfast::rmvnorm(1,rep(0,2*n_areas),Var_w)
  b <- w_draw[,1:(n_areas)]
  #structured component
  #u <- w_draw[,(n_areas+1):(2*n_areas)]
  #unstructured component
  #v <- (sqrt(tau_b)*b - sqrt(phi)*u)/(sqrt(1-phi))
  
  #observed data
  dat <- data.frame(
    # (true) number of births in each cluster
    N = rpois(sum(n_clusters),500),
    # admin area of each cluster
    A = c(sapply(1:n_areas,function(i){rep(i,n_clusters[i])})))
  dat$Y <- sapply(1:nrow(dat),function(i){rpois(1,dat$N[i]*expit(a + b[dat$A[i]]+rnorm(1,0,0.05)))})
  dat$cluster <- 1:nrow(dat)
  return(dat)
}


# MCMC to get P(b)|Ys (BYM spatial effect)  -------

a_init <- -5
b_init <- rep(0,n_areas)
#hyperprior specification for BYM2 precision (tau_b)
pc_alpha_tau <- 0.01
pc_u_tau <- 1
#hyperprior specifications for BYM2 spatial variation (phi)
pc_alpha_phi <- 2/3
pc_u_phi <- 0.5
data <- dat
n_iter <- 10000

postsamp_b_y_bym2_MCMC <- function(a_init,b_init,pc_u_tau,pc_alpha_tau,pc_u_phi,pc_alpha_phi,data,Q,n_iter){
  log_target_dist <- function(a,b,tau,phi){
    #likelihood
    sum(sapply(1:n_areas,function(k){log(dpois(Y_sums[k],N_sums[k]*expit(a + b[k])))})) + 
      #prior on a
      log(dnorm(a,-5,10)) +
      #prior on b
      sum(dmvnorm(b,rep(0,n_areas),(1/sqrt(tau))*(sqrt(1-phi)*diag(n_areas)+phi*Q_scaled_inv),logged = T)) +
      #hyperprior on tau
      INLA:::inla.pc.dprec(tau,u=pc_u_tau,alpha=pc_alpha_tau,log = T) +
      #hyperprior on phi
      phi.log.prior(phi)
    
  }
  n_areas <- length(b_init)
  eigQ <- eigen(Q)
  Q_scaled <- inla.scale.model(Q, constr=list(A=t(eigQ$vectors[,eigQ$values<1e-10]), e=rep(0,sum(eigQ$values<1e-10))))
  Q_scaled_inv <- ginv(as.matrix(Q_scaled))
  Y_sums <- (data %>% group_by(A) %>% summarise(sum(Y)))$`sum(Y)`
  N_sums <- (data %>% group_by(A) %>% summarise(sum(N)))$`sum(N)`
  phi.log.prior <- INLA:::inla.pc.bym.phi(Q,rankdef = 1,u=pc_u_phi,alpha = pc_alpha_phi)
  
  b_postsamp <- rbind(b_init, matrix(NA,n_iter-1,n_areas))
  a_postsamp <- c(a_init,rep(NA,n_iter-1))
  tau_draw <- c(NA,rep(NA,n_iter-1))
  phi_draw <- c(NA,rep(NA,n_iter-1))
  for(i in 2:n_iter){
    #not accepting/rejecting, just taking draws (will fix in Stan)
    tau_draw[i] <- INLA:::inla.pc.rprec(1,pc_u_tau,pc_alpha_tau)
    phi_draw[i] <- runif(1)
    
    a_current <- a_postsamp[i-1]
    a_proposed <- rnorm(1,a_postsamp[i-1],1)
    a_prob <- exp(log_target_dist(a_proposed,b_postsamp[(i-1),],tau_draw[i],phi_draw[i]) - 
                    log_target_dist(a_current,b_postsamp[(i-1),],tau_draw[i],phi_draw[i]))
    u <- runif(1)
    if(u < a_prob){
      a_postsamp[i] <- a_proposed
    }else if(u>= a_prob){
      a_postsamp[i] <- a_current
    }
    for(j in 1:n_areas){
      if(j==1){
        b_current <- b_postsamp[(i-1),]
      }else{
        b_current <- c(b_postsamp[i,(1:(j-1))],b_postsamp[(i-1),(j:n_areas)])
      }
      b_proposed <- b_current
      b_proposed[j] <- rnorm(1,b_current[j],0.25)
      a_prob <- exp(log_target_dist(a_postsamp[i],b_proposed,tau_draw[i],phi_draw[i]) - 
                      log_target_dist(a_postsamp[i],b_current,tau_draw[i],phi_draw[i]))
      u <- runif(1)
      if(u < a_prob){
        b_postsamp[i,] <- b_proposed
      }else if(u>= a_prob){
        b_postsamp[i,] <- b_current
      }
    }
  }  
  #get rid of burn-in
  a_postsamp <- a_postsamp[1001:n_iter]
  b_postsamp <- b_postsamp[1001:n_iter,]
  return(list(b=b_postsamp,theta=expit(a_postsamp + b_postsamp)))
}

postsamp_b <- postsamp_b_y_bym2_MCMC(-5,rep(0,n_areas),1,0.01,0.5,2/3,dat,Q,10000)

plot(postsamp_b$b[,2],type='l')
plot(postsamp_b$theta[,2],type='l')
plot(postsamp_b$theta[,5],type='l')
sapply(1:n_areas,function(i){hist(postsamp_b$theta[,i])})


# Simulate data (Ys) from ground truth data (Zs) -------  

#for now assume each cluster is chosen from area with equal probability
# number of clusters to sample from each area
samp_cluster_num <- rep(5,n_areas)

for(area in 1:n_areas){
  dat_tmp <- dat[dat$A==area,]
  #sample clusters from that area
  dat_tmp <- dat_tmp[sample(1:nrow(dat_tmp),samp_cluster_num[area]),]
  #generate number of births sampled in each cluster
  n <- rpois(nrow(dat_tmp),100)
  y <- rep(NA,nrow(dat_tmp))
  for(i in 1:nrow(dat_tmp)){
    y[i] <- rbinom(1,n[i],dat_tmp$Z[i]/dat_tmp$N[i])
  }
  dat_obs_tmp <- data.frame(cluster=dat_tmp$cluster,A=dat_tmp$A,n=n,Y=y)
  if(area==1){
    dat_obs <- dat_obs_tmp
  }else{
    dat_obs <- rbind(dat_obs,dat_obs_tmp)
  }
}











