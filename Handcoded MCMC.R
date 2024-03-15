library(sf)
library(SUMMER)
library(tidyverse)
library(data.table)
library(parallel)
library(LaplacesDemon)
library(VGAM)
library(Matrix)
library(Rfast)
library(MASS)
library(locfit)

code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]
setwd(paste(code.path.splitted[1: (length(code.path.splitted)-1)], collapse = "/"))

## helper functions -----
source('NUTS helpers.R')
source('nuts.R')

W_fn <- function(phi,Q_inv){
  return(diag(1-phi,nrow(Q_inv)) + phi*Q_inv)
}

log_target_5.2a_reg_params_4NUTS <- function(params,cons_params){
  tau <- exp(params[1])
  phi <- expit(params[2])
  W <- W_fn(phi,cons_params$Q_scaled_inv)
  d <- exp(params[3])
  alphaU <- params[4]
  alphaR <- params[5]
  b <- params[6:(length(params)-1)]
  Yplus <- params[length(params)]
  
  Nr_obs <- cons_params$data$N*exp(cons_params$data$U*alphaU + (1-cons_params$data$U)*alphaR + b[cons_params$data$A])
  Nr_sum_obs <- sum(Nr_obs)
  Nr_sum <- sum(cons_params$pop_dat$N*exp(cons_params$pop_dat$U*alphaU + (1-cons_params$pop_dat$U)*alphaR + b[cons_params$pop_dat$A]))
  
  term1 <- Yplus*log(d/(1+d)) - log(1+d)*Nr_sum/d #Negbin likelihood
  term2 <- lgamma(Yplus - sum(cons_params$Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - lgamma((1/d)*(Nr_sum-Nr_sum_obs)) - lgamma(Yplus - sum(cons_params$Y) + 1)
  term3 <- sum(lgamma(cons_params$Y+1/d*Nr_obs) - lgamma(1/d*Nr_obs))
  term4 <- -0.5*tau*t(b)%*%solve(W)%*%b # BYM2
  term5 <- -1/(2*cons_params$logit_r_hat_var)*(logitlink(Yplus/sum(cons_params$pop_dat$N)) - cons_params$logit_r_hat)^2 
  term6 <- -(alphaU-cons_params$mu_U)^2/(2*cons_params$sd_U^2) - (alphaR-cons_params$mu_R)^2/(2*cons_params$sd_R^2) - (log(d)-cons_params$mu_d)^2/(2*cons_params$sd_d^2) #alpha and d priors
  term7 <- (length(b)-1)/2*log(tau) + log(cons_params$alpha_tau)/cons_params$U_tau*tau^(-0.5) - 0.5*logdet(W) + cons_params$alpha_phi*(log(phi)) + cons_params$beta_phi*log(1-phi) #tau and phi likelihood + priors
  return(term1+term2+term3+term4+term5+term6+term7)
}

grad_log_target_5.2a_4NUTS <- function(params,cons_params){
  tau <- exp(params[1])
  phi <- expit(params[2])
  W <- W_fn(phi,cons_params$Q_scaled_inv)
  d <- exp(params[3])
  alphaU <- params[4]
  alphaR <- params[5]
  b <- params[6:(length(params)-1)]
  Yplus <- params[length(params)]
  
  Nr_obs <- cons_params$data$N*exp(cons_params$data$U*alphaU + (1-cons_params$data$U)*alphaR + b[cons_params$data$A])
  Nr_obs_urban <- cons_params$data$U*Nr_obs
  Nr_obs_rural <- (1-cons_params$data$U)*Nr_obs
  
  Nr_sum_urban <- sum(cons_params$pop_dat$U*cons_params$pop_dat$N*exp(alphaU + b[cons_params$pop_dat$A]))
  Nr_sum_rural <- sum((1-cons_params$pop_dat$U)*cons_params$pop_dat$N*exp(alphaR + b[cons_params$pop_dat$A]))
  Nr_sum <- Nr_sum_urban + Nr_sum_rural
  
  #digamma terms to be included in gradients
  digamma_1 <- digamma(Yplus - sum(cons_params$Y) + 1/d*(Nr_sum - sum(Nr_obs))) - digamma(1/d*(Nr_sum - sum(Nr_obs)))
  digamma_2 <- digamma(cons_params$Y + 1/d*Nr_obs) - digamma(1/d*Nr_obs)
  
  grad_logtau <- 0.5*(length(b)-1) - 0.5*log(cons_params$alpha_tau)/cons_params$U_tau*tau^(-0.5) - 0.5*tau*t(b)%*%solve(W)%*%b
  
  grad_logitphi <- -0.5*phi*(1-phi)*tr(solve(W)%*%(cons_params$Q_scaled_inv - diag(1,length(b)))) + 
    0.5*tau*phi*(1-phi)*t(b)%*%solve(W)%*%(cons_params$Q_scaled_inv - diag(1,length(b)))%*%solve(W)%*%b +
    cons_params$alpha_phi*(1-phi) - cons_params$beta_phi*phi

  grad_logd <- Yplus*(1/(1+d)) - Nr_sum*(1/(1+d)-log(1+d)/d) - (1/d)*(Nr_sum - sum(Nr_obs))*digamma_1 + 
    (1/d)*sum(Nr_obs*(digamma(Nr_obs/d) - digamma(cons_params$Y + Nr_obs/d))) - (1/cons_params$sd_d^2)*(log(d) - cons_params$mu_d)
  
  grad_alphaU <- -log(1+d)/d*Nr_sum_urban + 1/d*(Nr_sum_urban - sum(Nr_obs_urban))*digamma_1 +1/d*sum(Nr_obs_urban*digamma_2) - (1/cons_params$sd_U^2)*(alphaU - cons_params$mu_U)
  
  grad_alphaR <- -log(1+d)/d*Nr_sum_rural + 1/d*(Nr_sum_rural - sum(Nr_obs_rural))*digamma_1 + 1/d*sum(Nr_obs_rural*digamma_2) - (1/cons_params$sd_R^2)*(alphaR - cons_params$mu_R)
  
  grad_b <- sapply(1:length(b),function(area){
    Nr_obs_area <- I(cons_params$data$A==area)*Nr_obs
    Nr_sum_area <- sum(I(cons_params$pop_dat$A==area)*cons_params$pop_dat$N*exp(cons_params$pop_dat$U*alphaU + (1-cons_params$pop_dat$U)*alphaR + b[area]))
    grad_tmp <- -log(1+d)/d*Nr_sum_area + 1/d*(Nr_sum_area - sum(Nr_obs_area))*digamma_1 + 1/d*sum(Nr_obs_area*digamma_2)
    return(grad_tmp)
  })
  grad_b <- grad_b -tau*solve(W)%*%b
  
  grad_Yplus <- log(d/(1+d)) + digamma(Yplus - sum(cons_params$Y) + (1/d)*(Nr_sum-sum(Nr_obs))) - digamma(Yplus - sum(cons_params$Y) + 1) -
    1/cons_params$logit_r_hat_var*(logitlink(Yplus/sum(cons_params$pop_dat$N)) - cons_params$logit_r_hat)*sum(cons_params$pop_dat$N)/(Yplus*(sum(cons_params$pop_dat$N)-Yplus))
 
  return(c(grad_logtau,grad_logitphi,grad_logd,grad_alphaU,grad_alphaR,grad_b,grad_Yplus))
}

log_target_cond_5.2_y <- function(y,data,Yplus,pop_dat,Y,i,Nr_sum,Nr_obs,d){
  term1 <- lgamma(Yplus - sum(Y[-i]) - y + 1/d*(Nr_sum - sum(Nr_obs))) - lgamma(Yplus - sum(Y[-i]) -y + 1)
  term2 <- lgamma(data$N[i]-y+1) - lgamma(y-data$Z[i]+1) - lgamma(data$N[i]-y-data$n[i]+data$Z[i] +1)
  term3 <- lgamma(y + 1/d*Nr_obs[i])
  
  return(term1 + term2 + term3)
}

log_target_5.3a_reg_params_4NUTS <- function(params,cons_params){
  tau <- exp(params[1])
  phi <- expit(params[2])
  W <- W_fn(phi,cons_params$Q_inv)
  W_inv <- solve(W)
  d <- exp(params[3])
  alphaU <- params[4]
  alphaR <- params[5]
  b <- params[6:(length(params))]
  #b <- params[6:(length(params)-length(cons_params$Y_sum_obs))]
  #Yplus <- tail(params,length(cons_params$Y_sum_obs))
  
  # add Nr_obs to data
  cons_params$data_table[, ':='(Nr_obs=N*exp(U*alphaU + (1-U)*alphaR + b[A2]))]
  
  # get grouped sums faster w data.table package
  N_sum <- (cons_params$pop_dt[, sum(N), by=A1][order(A1)])$V1
  Nr_sum <- (cons_params$pop_dt[, sum(N*exp(U*alphaU + (1-U)*alphaR + b[A2])), by=A1][order(A1)])$V1
  Nr_sum_obs <- (cons_params$data_table[, sum(Nr_obs),by=A1][order(A1)])$V1
  
  out <- sum(cons_params$Yplus)*log(d/(1+d)) - sum(log(1+d)*Nr_sum/d) -
    sum(lgamma(cons_params$Yplus - cons_params$Y_sum_obs + 1)) +
    sum(lgamma(cons_params$Yplus - cons_params$Y_sum_obs + (1/d)*(Nr_sum-Nr_sum_obs))- lgamma((1/d)*(Nr_sum-Nr_sum_obs))) +
    sum(lgamma(cons_params$data_table$Y+1/d*cons_params$data_table$Nr_obs) - lgamma(1/d*cons_params$data_table$Nr_obs)) +
    -0.5*sum(1/cons_params$logit_r_hat_var*(logitlink(cons_params$Yplus/N_sum) - cons_params$logit_r_hat)^2) +
    -0.5*tau*t(b)%*%W_inv%*%b -
    (alphaU-cons_params$mu_U)^2/(2*cons_params$sd_U^2) - (alphaR-cons_params$mu_R)^2/(2*cons_params$sd_R^2) -
    (log(d)-cons_params$mu_logd)^2/(2*cons_params$sd_logd^2) +
    (length(b)-1)/2*log(tau) + log(cons_params$alpha_tau)/cons_params$U_tau*tau^(-0.5) - 0.5*logdet(W) + cons_params$alpha_phi*(log(phi)) + cons_params$beta_phi*log(1-phi)

   return(out)
}

grad_log_target_5.3a_4NUTS <- function(params,cons_params){
  tau <- exp(params[1])
  phi <- expit(params[2])
  W <- W_fn(phi,cons_params$Q_inv)
  W_inv <- solve(W)
  d <- exp(params[3])
  alphaU <- params[4]
  alphaR <- params[5]
  b <- params[6:(length(params))]
  #b <- params[6:(length(params)-length(cons_params$Y_sum_obs))]
  #Yplus <- tail(params,length(cons_params$Y_sum_obs))
  
  # add Nr_obs to data
  cons_params$data_table[, ':='(Nr_obs_urban=U*N*exp(U*alphaU + b[A2]))]
  cons_params$data_table[, ':='(Nr_obs_rural=(1-U)*N*exp((1-U)*alphaR + b[A2]))]
  Nr_sum_obs_urban <- (cons_params$data_table[, sum(Nr_obs_urban),by=A1][order(A1)])$V1
  Nr_sum_obs_rural <- (cons_params$data_table[, sum(Nr_obs_rural),by=A1][order(A1)])$V1
  Nr_sum_obs <- (cons_params$data_table[, sum(Nr_obs_urban+Nr_obs_rural),by=A1][order(A1)])$V1
  
  # get grouped sums faster w data.table package
  N_sum <- (cons_params$pop_dt[, sum(N), by=A1][order(A1)])$V1
  
  Nr_sum_urban <- (cons_params$pop_dt[U==1, sum(N*exp(U*alphaU + (1-U)*alphaR + b[A2])), by=A1][order(A1)])$V1
  Nr_sum_rural <- (cons_params$pop_dt[U==0, sum(N*exp(U*alphaU + (1-U)*alphaR + b[A2])), by=A1][order(A1)])$V1
  Nr_sum <- (cons_params$pop_dt[, sum(N*exp(U*alphaU + (1-U)*alphaR + b[A2])), by=A1][order(A1)])$V1
  
  Nr_sum_area <- cons_params$pop_dt[,sum(N*exp(U*alphaU + (1-U)*alphaR + b[A2])), by=A2][order(A2)]$V1
  Nr_sum_obs_area <- (cons_params$data_table[, sum(Nr_obs_urban+Nr_obs_rural),by=A2][order(A2)])$V1
  
  # digamma terms to be included in gradients
  digamma_1 <- digamma(cons_params$Yplus - cons_params$Y_sum_obs + 1/d*(Nr_sum - Nr_sum_obs)) - digamma(1/d*(Nr_sum - Nr_sum_obs))
  digamma_2 <- digamma(cons_params$data_table$Y + 1/d*(cons_params$data_table$Nr_obs_urban+cons_params$data_table$Nr_obs_rural)) - 
    digamma(1/d*(cons_params$data_table$Nr_obs_urban+cons_params$data_table$Nr_obs_rural))
  cons_params$data_table$digamma_2 <- digamma_2
  
  grad_logtau <- 0.5*(length(b)-1) - 0.5*log(cons_params$alpha_tau)/cons_params$U_tau*tau^(-0.5) - 0.5*tau*t(b)%*%solve(W)%*%b
  
   grad_logitphi <- -0.5*phi*(1-phi)*tr(solve(W)%*%(cons_params$Q_inv - diag(1,length(b)))) + 
     0.5*tau*phi*(1-phi)*t(b)%*%solve(W)%*%(cons_params$Q_inv - diag(1,length(b)))%*%solve(W)%*%b +
     cons_params$alpha_phi*(1-phi) - cons_params$beta_phi*phi
  
  grad_logd <- sum(cons_params$Yplus)*(1/(1+d)) - sum(Nr_sum)*(1/(1+d)-log(1+d)/d) - (1/d)*sum((Nr_sum - Nr_sum_obs)*digamma_1) +
    (1/d)*sum((cons_params$data_table$Nr_obs_urban + cons_params$data_table$Nr_obs_rural)*(digamma((cons_params$data_table$Nr_obs_urban + cons_params$data_table$Nr_obs_rural)/d) - digamma(cons_params$data_table$Y + (cons_params$data_table$Nr_obs_urban + cons_params$data_table$Nr_obs_rural)/d))) - 
    (1/cons_params$sd_logd^2)*(log(d) - cons_params$mu_logd)

  grad_alphaU <- sum(-log(1+d)/d*Nr_sum_urban + 1/d*(Nr_sum_urban - Nr_sum_obs_urban)*digamma_1) +
    1/d*sum(cons_params$data_table$Nr_obs_urban*digamma_2) - 
    (1/cons_params$sd_U^2)*(alphaU - cons_params$mu_U)

  grad_alphaR <- sum(-log(1+d)/d*Nr_sum_rural + 1/d*(Nr_sum_rural - Nr_sum_obs_rural)*digamma_1) +
    1/d*sum(cons_params$data_table$Nr_obs_rural*digamma_2) - 
    (1/cons_params$sd_R^2)*(alphaR - cons_params$mu_R)
  
  grad_b <- -log(1+d)/d*Nr_sum_area +
    1/d*(Nr_sum_area - Nr_sum_obs_area)*digamma_1[unique(cons_params$pop_dt[order(A2),c('A1','A2')])$A1] +
    1/d*(cons_params$data_table[, sum(digamma_2*(Nr_obs_urban+Nr_obs_rural)),by=A2][order(A2)])$V1 - 
    tau*solve(W)%*%b
  
  return(c(grad_logtau,grad_logitphi,grad_logd,grad_alphaU,grad_alphaR,grad_b))
}

log_target_cond_5.3_Yplus <- function(Yplus,pop_dat,Y,Nr_sum,Nr_obs,d,logit_r_hat,logit_r_hat_var){
  term1 <- -1/(2*logit_r_hat_var)*(logitlink(Yplus/sum(pop_dat$N)) - logit_r_hat)^2 
  term2 <- Yplus*log(d/(1+d)) #negative linear
  term3 <- lgamma(Yplus - sum(Y) + (1/d)*(Nr_sum-sum(Nr_obs))) - lgamma(Yplus - sum(Y) + 1)
  
  return(term1 + term2 + term3)
}

log_target_cond_5.3_y <- function(y,data,Yplus,pop_dat,Y,i,Nr_sum,Nr_obs,d){
  term1 <- lgamma(Yplus - sum(Y[-i]) - y + 1/d*(Nr_sum - sum(Nr_obs))) - lgamma(Yplus - sum(Y[-i]) -y + 1)
  term2 <- lgamma(data$N[i]-y+1) - lgamma(y-data$Z[i]+1) - lgamma(data$N[i]-y-data$n[i]+data$Z[i] +1)
  term3 <- lgamma(y + 1/d*Nr_obs[i])
  
  return(term1 + term2 + term3)
}


# more helper functions ----------
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
log_target_5.2_reg_params <- function(params,data,Yplus,pop_dat,Y,tau,d,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R){
  alphaU <- params[1]
  alphaR <- params[2]
  b <- params[3:length(params)]
  
  Nr_obs <- data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A])
  Nr_sum_obs <- sum(Nr_obs)
  Nr_sum <- sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  term1 <- Yplus*log(d/(1+d)) - log(1+d)*Nr_sum/d
  term2 <- lgamma(Yplus - sum(Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - lgamma((1/d)*(Nr_sum-Nr_sum_obs)) 
  term3 <- sum(lgamma(Y+1/d*Nr_obs) - lgamma(1/d*Nr_obs))
  term4 <- log(d) - tau*sum(b^2)/2 - (alphaU-mu_U)^2/(2*sd_U^2) - (alphaR-mu_R)^2/(2*sd_R^2) - (log(d)-mu_d)^2/(2*sd_d^2)
  
  return(term1 + term2 + term3 + term4)
}

log_target_5.2a_reg_params <- function(params,data,Yplus,pop_dat,Y,tau,W_inv,mu_d,sd_d,mu_U,sd_U,mu_R,sd_R){
  d <- exp(params[1])
  alphaU <- params[2]
  alphaR <- params[3]
  b <- params[4:length(params)]
  
  Nr_obs <- data$N*exp(data$U*alphaU + (1-data$U)*alphaR + b[data$A])
  Nr_sum_obs <- sum(Nr_obs)
  Nr_sum <- sum(pop_dat$N*exp(pop_dat$U*alphaU + (1-pop_dat$U)*alphaR + b[pop_dat$A]))
  
  term1 <- Yplus*log(d/(1+d)) - log(1+d)*Nr_sum/d
  term2 <- lgamma(Yplus - sum(Y) + (1/d)*(Nr_sum-Nr_sum_obs)) - lgamma((1/d)*(Nr_sum-Nr_sum_obs)) 
  term3 <- sum(lgamma(Y+1/d*Nr_obs) - lgamma(1/d*Nr_obs))
  term4 <- log(d) - 0.5*tau*t(b)%*%W_inv%*%b - (alphaU-mu_U)^2/(2*sd_U^2) - (alphaR-mu_R)^2/(2*sd_R^2) - (log(d)-mu_d)^2/(2*sd_d^2)
  
  return(term1 + term2 + term3 + term4)
}

log_target_5.2_logd <- function(d,Yplus,Nr_sum,Nr_obs,Y,mu_logd,sd_logd){
  
  term1 <- Yplus*log(d/(1+d)) - log(1+d)*Nr_sum/d
  term2 <- lgamma(Yplus - sum(Y) + (1/d)*(Nr_sum-sum(Nr_obs))) - lgamma((1/d)*(Nr_sum-sum(Nr_obs))) 
  term3 <- sum(lgamma(Y+1/d*Nr_obs) - lgamma(1/d*Nr_obs))
  term4 <- log(d) - (log(d)-mu_d)^2/(2*sd_d^2)
  
  return(term1 + term2 + term3 + term4)
}

log_target_cond_5.2_Yplus <- function(Yplus,pop_dat,Y,Nr_sum,Nr_obs,d,logit_r_hat,logit_r_hat_var){
  term1 <- -1/(2*logit_r_hat_var)*(logitlink(Yplus/sum(pop_dat$N)) - logit_r_hat)^2 
  term2 <- Yplus*log(d/(1+d)) #negative linear
  term3 <- lgamma(Yplus - sum(Y) + (1/d)*(Nr_sum-sum(Nr_obs))) - lgamma(Yplus - sum(Y) + 1)
  
  return(term1 + term2 + term3)
}

log_target_5.3_reg_params <- function(params,data_table,Yplus,pop_dt,d,Y_sum_obs,tau,mu_U,sd_U,mu_R,sd_R){
  alphaU <- params[1]
  alphaR <- params[2]
  b <- params[3:length(params)]
  
  # add Nr_obs to data
  data_table[, ':='(Nr_obs=N*exp(U*alphaU + (1-U)*alphaR + b[A2]))]
  
  # get grouped sums faster w data.table package
  Nr_sum <- (pop_dt[, sum(N*exp(U*alphaU + (1-U)*alphaR + b[A2])), by=A1])$V1
  Nr_sum_obs <- (data_table[, sum(Nr_obs),by=A1])$V1
  
  out <- sum(-log(1+d)*Nr_sum/d) +
    sum(lgamma(Yplus - Y_sum_obs + (1/d)*(Nr_sum-Nr_sum_obs))- lgamma((1/d)*(Nr_sum-Nr_sum_obs))) + 
    sum(lgamma(data_table$Y+1/d*data_table$Nr_obs) - lgamma(1/d*data_table$Nr_obs)) + 
    (-1)*tau*sum(b^2)/2 - (alphaU-mu_U)^2/(2*sd_U^2) - (alphaR-mu_R)^2/(2*sd_R^2)
  
  return(out)
}

log_target_5.3_logd <- function(log_d,alphaU,alphaR,b,data,Yplus,pop_dat,Y,mu_logd,sd_logd){
  d <- exp(log_d)
  
  likelihood <- sapply(unique(data$A1),function(i){
    sub_data <- data[data$A1==i,]
    sub_pop_dat <- pop_dat[pop_dat$A1==i,]
    
    Nr_obs <- sub_data$N*exp(sub_data$U*alphaU + (1-sub_data$U)*alphaR + b[sub_data$A2])
    Nr_sum_obs <- sum(Nr_obs)
    Nr_sum <- sum(sub_pop_dat$N*exp(sub_pop_dat$U*alphaU + (1-sub_pop_dat$U)*alphaR + b[sub_pop_dat$A2]))
    
    term1 <- Yplus[i]*log(d/(1+d)) - log(1+d)*Nr_sum/d
    term2 <- lgamma(Yplus[i] - sum(Y[[i]]) + (1/d)*(Nr_sum-Nr_sum_obs)) - lgamma((1/d)*(Nr_sum-Nr_sum_obs)) 
    term3 <- sum(lgamma(Y[[i]]+1/d*Nr_obs) - lgamma(1/d*Nr_obs))
    
    return(term1 + term2 + term3)
  })
  
  priors <- -(log_d-mu_logd)^2/(2*sd_logd^2)
  
  return(sum(likelihood) + priors)
}

log_target_cond_5.3_Yplus <- function(Yplus,pop_dat,Y,Nr_sum,Nr_obs,d,logit_r_hat,logit_r_hat_var){
  term1 <- -1/(2*logit_r_hat_var)*(logitlink(Yplus/sum(pop_dat$N)) - logit_r_hat)^2 
  term2 <- Yplus*log(d/(1+d)) #negative linear
  term3 <- lgamma(Yplus - sum(Y) + (1/d)*(Nr_sum-sum(Nr_obs))) - lgamma(Yplus - sum(Y) + 1)
  
  return(term1 + term2 + term3)
}

log_target_5.3a_reg_params <- function(params,data_table,Yplus,pop_dt,d,Y_sum_obs,tau,W_inv,mu_U,sd_U,mu_R,sd_R){
  alphaU <- params[1]
  alphaR <- params[2]
  b <- params[3:length(params)]
  
  # add Nr_obs to data
  data_table[, ':='(Nr_obs=N*exp(U*alphaU + (1-U)*alphaR + b[A2]))]
  
  # get grouped sums faster w data.table package
  Nr_sum <- (pop_dt[, sum(N*exp(U*alphaU + (1-U)*alphaR + b[A2])), by=A1])$V1
  Nr_sum_obs <- (data_table[, sum(Nr_obs),by=A1])$V1
  
  out <- sum(-log(1+d)*Nr_sum/d) +
    sum(lgamma(Yplus - Y_sum_obs + (1/d)*(Nr_sum-Nr_sum_obs))- lgamma((1/d)*(Nr_sum-Nr_sum_obs))) + 
    sum(lgamma(data_table$Y+1/d*data_table$Nr_obs) - lgamma(1/d*data_table$Nr_obs)) + 
    (-1/2)*tau*t(b)%*%W_inv%*%b - (alphaU-mu_U)^2/(2*sd_U^2) - (alphaR-mu_R)^2/(2*sd_R^2)
  
  return(out)
}


#######################################################################
#######################################################################
# Alg 5.2 (Conditioning on Natl, IID): P(alphaU,alphaR,b,tau,d,Y^s,Y+|r+_est,Z,delta) ------

# tau_init=40
# d_init = 1
# alphaU_init = seq(-4,-3,by=0.1)
# alphaR_init= seq(-3.5,-2.5,by=0.1) #initial values
# a_tau = 0.01
# b_tau = 0.01
# mu_U = -3.5
# sd_U = 3
# mu_R = -3.2
# sd_R = 3
# mu_logd = 0
# sd_logd = 1#hyperpriors  
# 
# eps=0.6
# prop_sd_logd=0.3
# logit_r_hat=dir.est$logit.est
# logit_r_hat_var = dir.est$var.est
# n_iter=100
# chains = 4
# burnin=00

#preconditioned MALA as described in "Adaptive step size selection for Hessian-based manifold Langevian samplers" Tore Selland Kleppe
postsamp_Alg5.2_MALA_MCMC <- function(tau_init=40, d_init = 1, alphaU_init = seq(-4,-3,by=0.1), alphaR_init= seq(-3.5,-2.5,by=0.1), #initial values
                                      a_tau = 0.01,b_tau = 0.01,mu_U = -3.5,sd_U = 3,mu_R = -3.2,sd_R = 3,mu_logd = 0,sd_logd = 1,#hyperpriors  
                                      eps, prop_sd_logd, #tuning params for regression param and Yplus proposals
                                      logit_r_hat,logit_r_hat_var,data_list, 
                                      n_iter, burnin=n_iter/2,chains=5){
  
  data <- data_list$obs_dat
  data <- data[order(data$A),]
  pop_dat <- data_list$pop_strata_dat
  n_areas <- length(unique(data_list$pop_strata_dat$A))
  
  parallel_results <- mclapply(1:chains, function(j){
    
    b_init <- rnorm(n_areas,0,1/sqrt(tau_init))
    Yplus_init <- round(sum(pop_dat$N)*expit(logit_r_hat))
    
    ## ensure we are starting with a PD negative hessian 
    reg_params_init_grid <- as.matrix(expand.grid(alphaU_init,alphaR_init))
    search_order <- sample(1:nrow(reg_params_init_grid),nrow(reg_params_init_grid))
    
    for(j in 1:nrow(reg_params_init_grid)){
      reg_params_init <-  c(reg_params_init_grid[search_order[j],],b_init)
      Nr_sum_init <- sum(pop_dat$N*exp(pop_dat$U*reg_params_init[1] + (1-pop_dat$U)*reg_params_init[2] + reg_params_init[2+pop_dat$A]))
      Nr_obs_init <- data$N*exp(data$U*reg_params_init[1] + (1-data$U)*reg_params_init[2] + reg_params_init[2+data$A])
      Y_init <- rmultinom(1,Yplus_init,prob=c(Nr_obs_init,Nr_sum_init-sum(Nr_obs_init)))
      Y_init <- Y_init[-length(Y_init)]
      Y_init[Y_init<data$Z] <- data_list$obs_dat$Z[Y_init<data$Z]
      
      hess_init <- pracma::hessian(log_target_5.2_reg_params,x=reg_params_init,data=data,Yplus = Yplus_init,pop_dat=pop_dat,
                                   Y=Y_init,tau=tau_init,d=d_init,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      if(is.positive.definite(-hess_init)){
        break
      }
    }
    message('Initial values obtained.')
    
    ## initialize data frames 
    tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
    d_postsamp_t <- c(d_init, rep(NA,n_iter-1))
    reg_params_postsamp_t <- rbind(reg_params_init, matrix(NA,n_iter-1,length(reg_params_init)))
    Y_postsamp_t <- rbind(Y_init, matrix(NA,n_iter-1,nrow(data)))
    Yplus_postsamp_t <- c(Yplus_init, rep(NA,n_iter-1))
    acc_reg <- acc_d <- 0
    
    # set current states
    tau_current <- tau_init
    d_current <- d_init
    alphaU_current <- reg_params_init[1]
    alphaR_current <- reg_params_init[2]
    b_current <- b_init
    reg_params_current <- reg_params_init
    Y_current <- Y_init
    Yplus_current <- Yplus_init
    
    #start.time <- Sys.time()
    for(i in 2:n_iter){
      
      # draw new tau from full conditional ----
      tau_current <- rgamma(1,n_areas/2+a_tau,0.5*sum(b_current^2)+b_tau)
      
      # draw proposal for regression params ----
      grad_current_reg_params <- pracma::grad(log_target_5.2_reg_params,x=reg_params_current,data=data,Yplus = Yplus_current, pop_dat=pop_dat,
                                              Y=Y_current,tau=tau_current,d=d_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      
      hess_current_reg_params <- pracma::hessian(log_target_5.2_reg_params,x=reg_params_current,data=data,Yplus = Yplus_current,pop_dat=pop_dat,
                                                 Y=Y_current,tau=tau_current,d=d_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      if(is.positive.definite(-hess_current_reg_params)){
        cov_current_reg_params <- matlib::inv(-hess_current_reg_params)
      }else{
        message(paste0('current Hessian is not PD at iteration ',i))
        L_current_reg_params <- get.GMW.chol(-hess_current_reg_params,0.01)
        cov_current_reg_params <- matlib::inv(L_current_reg_params%*%t(L_current_reg_params))
      }
      
      
      reg_params_proposal <- as.vector(reg_params_current + eps^2/2*cov_current_reg_params%*%grad_current_reg_params +
                                         t(Rfast::rmvnorm(1,rep(0,length(reg_params_current)),eps^2*cov_current_reg_params))) 
      
      
      grad_proposal_reg_params <- pracma::grad(log_target_5.2_reg_params,x=reg_params_proposal,data=data,Yplus = Yplus_current,pop_dat=pop_dat,
                                               Y=Y_current,tau=tau_current,d=d_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      
      hess_proposal_reg_params <- pracma::hessian(log_target_5.2_reg_params,x=reg_params_proposal,data=data,Yplus = Yplus_current,pop_dat=pop_dat,
                                                  Y=Y_current,tau=tau_current,d=d_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      if(is.symmetric(hess_proposal_reg_params)){
        if(is.positive.definite(-hess_proposal_reg_params)){
          cov_proposal_reg_params <- matlib::inv(-hess_proposal_reg_params)
          skip.it <- F
          
        }else{
          message(paste0('Proposed Hessian not PD at iteration ',i,'\n'))
          skip.it <- T}
      }else{
        message(paste0('Proposed Hessian not symmetric at iteration ',i,'\n'))
        skip.it <- T}
      
      
      if(!skip.it){
        # accept or reject proposal for (alpha, b, d)
        
        a_prob_reg <- min(1,exp({
          log_target_5.2_reg_params(reg_params_proposal,data=data,Yplus = Yplus_current,pop_dat=pop_dat,Y=Y_current,tau=tau_current,d=d_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R) +
            Rfast::dmvnorm(reg_params_current,reg_params_proposal + eps^2/2*cov_proposal_reg_params%*%grad_proposal_reg_params,eps^2*cov_proposal_reg_params, log=T) - #optimize
            log_target_5.2_reg_params(reg_params_current,data=data,Yplus = Yplus_current,pop_dat=pop_dat,Y=Y_current,tau=tau_current,d=d_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R) -
            Rfast::dmvnorm(reg_params_proposal,reg_params_current + eps^2/2*cov_current_reg_params%*%grad_current_reg_params,eps^2*cov_current_reg_params, log=T) #optimize
        }))
        
        if(a_prob_reg>runif(1)){
          reg_params_current <- reg_params_proposal
          alphaU_current <- reg_params_proposal[1]
          alphaR_current <- reg_params_proposal[2]
          b_current <- reg_params_proposal[3:length(reg_params_proposal)]
          if(i>burnin){acc_reg <- acc_reg + 1}
        }
      }
      
      Nr_sum_current <- sum(pop_dat$N*exp(pop_dat$U*alphaU_current + (1-pop_dat$U)*alphaR_current + b_current[pop_dat$A]))
      Nr_obs_current <- data$N*exp(data$U*alphaU_current + (1-data$U)*alphaR_current + b_current[data$A])
      
      #draw proposal for d -------
      logd_proposal <- rnorm(1,log(d_current),prop_sd_logd)
      
      a_prob_d <- exp({log_target_5.2_logd(exp(logd_proposal),Yplus_current,Nr_sum_current,Nr_obs_current,Y_current,mu_logd,sd_logd) -
          log_target_5.2_logd(d_current,Yplus_current,Nr_sum_current,Nr_obs_current,Y_current,mu_logd,sd_logd)})
      
      if(a_prob_d > runif(1)){
        d_current <- exp(logd_proposal)
        if(i>burnin){
          acc_d <- acc_d + 1
        }
      }
      # draw Yplus ----
      
      domain <- round(0.5*Yplus_current):(2*Yplus_current)
      f <- log_target_cond_5.2_Yplus(domain,pop_dat,Y_current,Nr_sum_current,Nr_obs_current,
                                     d_current,logit_r_hat,logit_r_hat_var)
      
      del_f <- f-mean(f) #normalize
      #restrict domain to densities >0
      zero_dens <- which(exp(del_f)==0)
      while(length(zero_dens)>0){
        domain <- domain[-zero_dens]
        f <- f[-zero_dens]
        del_f <- f-mean(f)
        zero_dens <- which(exp(del_f)==0)
      }
      
      target_density <- (exp(del_f)/sum(exp(del_f)))
      target_F <- cumsum(target_density)
      
      #take draw from inverse CDF
      U <- runif(1)
      Yplus_current <- (min(domain[target_F>U]))
      
      # draw Ys ----
      for(k in 1:nrow(data)){
        # derive CDF
        domain <- data$Z[k]:round(0.5*data$N[k])
        f <- log_target_cond_5.2_y(domain,data,Yplus_current,pop_dat,
                                   Y_current,k,Nr_sum_current,Nr_obs_current,d_current)
        domain <- domain[f>-Inf]
        f <- f[f>-Inf]
        del_f <- f-mean(f) #normalize
        #restrict domain to densities >0
        zero_dens <- which(exp(del_f)==0)
        while(length(zero_dens)>0){
          domain <- domain[-zero_dens]
          f <- f[-zero_dens]
          del_f <- f-mean(f)
          zero_dens <- which(exp(del_f)==0)
        }
        
        target_density <- (exp(del_f)/sum(exp(del_f)))
        target_F <- cumsum(target_density)
        #take draw from inverse CDF
        U <- runif(1)
        if(is.na(min(domain[target_F>U]))){
          stop()
        }
        Y_current[k] <- min(domain[target_F>U])
      }
      
      
      # record current state ----
      tau_postsamp_t[i] <- tau_current
      d_postsamp_t[i] <- d_current
      if(!skip.it){
        reg_params_postsamp_t[i,] <- reg_params_current
      }
      Y_postsamp_t[i,] <- Y_current
      Yplus_postsamp_t[i] <- Yplus_current
    } 
    #Sys.time() - start.time
    
    results <- list(reg_params_t=reg_params_postsamp_t[(burnin+1):n_iter,], 
                    d_t=d_postsamp_t[(burnin+1):n_iter], 
                    tau_t=tau_postsamp_t[(burnin+1):n_iter], 
                    Y_t=Y_postsamp_t[(burnin+1):n_iter,],
                    Yplus_t=Yplus_postsamp_t[(burnin+1):n_iter], 
                    acc_reg_t = acc_reg,acc_d_t = acc_d)
    
    return(results)
  },mc.cores = chains)
  
  #combine chains
  b_postsamp <- reg_params_postsamp <- Y_postsamp <- NULL
  for(j in 1:chains){
    b_postsamp <- rbind(b_postsamp,parallel_results[[j]]$reg_params_t[,3:(2+n_areas)])
    reg_params_postsamp <- rbind(reg_params_postsamp,parallel_results[[j]]$reg_params_t)
    Y_postsamp <- rbind(Y_postsamp,parallel_results[[j]]$Y_t)
  }
  tau_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$tau_t}))
  d_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$d_t}))
  alphaU_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,1]}))
  alphaR_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,2]}))
  Yplus_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$Yplus_t}))
  names(tau_postsamp) <- names(d_postsamp) <- names(alphaU_postsamp) <- names(alphaR_postsamp) <- NULL
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp,alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  
  acc_reg_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_reg_t})))/((n_iter-burnin)*chains)
  acc_d_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_d_t})))/((n_iter-burnin)*chains)
  
  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              tau=tau_postsamp, 
              reg_params=reg_params_postsamp,
              d=d_postsamp, 
              Y=Y_postsamp,
              Yplus=Yplus_postsamp,
              acc_reg_rate = acc_reg_rate,
              acc_d_rate = acc_d_rate))
}

#######################################################################
#######################################################################
# Simulations for Alg 5.2 -----------

## set populations and simulation parameters
{
  ## this section is calibrated to be relatively similar to Sierra Leone 2019 survey, in terms of sample sizes
  n_areas <- 20 #numbers of areas
  n_clusters_urban <- rnegbin(n_areas,200,5) #number of urban clusters in each area
  n_clusters_rural <- rnegbin(n_areas,300,5) #number of rural clusters in each area
  n_births_urban <- 75 #average number of births per urban cluster (0.1*number of households*3) <- assuming 3 year period
  n_births_rural <-  100 #average number of births per rural cluster (0.15*number of households*3) <- assuming 3 year period
  n_clusters_urban_samp <-  round(0.075*n_clusters_urban) #number of urban clusters sampled in each area
  n_clusters_rural_samp <-  round(0.05*n_clusters_rural) #number of rural clusters sampled in each area
  n_births_urban_samp <- 20 # average number of births sampled from each urban cluster (0.1*30*3) -- change to be a little smaller later
  n_births_rural_samp <- 30 # average number of births sampled from each rural cluster -- change to be a little smaller later
  
  # intercept (fixed effect) 
  alphaU <- rnorm(1,-3.5,0.25)
  alphaR <- rnorm(1,-3.2,0.25)
  
  # draw hyperparameters
  tau <- rgamma(1,10,.25)
  
  # draw overdispersion parameter
  d <- rlnorm(1,0,0.25)
  
  # draw random effects
  b <- as.vector(Rfast::rmvnorm(1,rep(0,n_areas),diag(1/tau,n_areas)))
  
  true_params <- c(log(d),alphaU,alphaR,b)
  true_r <- exp(c(alphaU +b,alphaR +b))
}

## generate data 
{
  all_dat <- data.frame(
    # (true) number of births in each cluster
    N = round(c(rnorm(sum(n_clusters_urban),n_births_urban,10),rnorm(sum(n_clusters_rural),n_births_rural,12))),
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
  
  # assign sampling weights
  #births in strata/(# cluster sampled * #births sampled in cluster)
  pop_dat <- all_dat %>% group_by(A,U) %>% summarise(N=sum(N))
  pop_dat$num_clust_samp <- pop_dat$num_births_samp <- NA
  pop_dat$num_clust_samp[pop_dat$U==0] <- n_clusters_rural_samp
  pop_dat$num_clust_samp[pop_dat$U==1] <- n_clusters_urban_samp
  pop_dat$num_births_samp[pop_dat$U==0] <- n_births_rural_samp
  pop_dat$num_births_samp[pop_dat$U==1] <- n_births_urban_samp
  pop_dat$wt <- pop_dat$N/(pop_dat$num_births_samp*pop_dat$num_clust_samp)
  pop_dat$wt <- pop_dat$wt/100
  obs_dat <- merge(obs_dat,pop_dat[c('A','U','wt')])
  obs_dat <- obs_dat[order(obs_dat$A),]
  
  data_list <- list(obs_dat = obs_dat, #observed data
                    Yplus = sum(all_dat$Y), #ALL deaths
                    pop_strata_dat = pop_dat) # number of births per strata (area x urban)
  
  #what is the implied national urban and rural rates
  pop_dat <- pop_dat %>% mutate(alpha=ifelse(U==1,alphaU,alphaR))
  pop_dat$b <- b[pop_dat$A]
  pop_dat$r <- exp(pop_dat$alpha+pop_dat$b)
  natl_r <- as.numeric(pop_dat %>% ungroup() %>% summarise(r = sum(r*N)/sum(N)))
  
}

## get national direct estimate
{
  #obtain a direct estimate and variance for national rate
  data_for_direct <- data_list$obs_dat
  data_for_direct$years <- 'All'
  data_for_direct$died <- data_for_direct$Z
  data_for_direct$strata <- paste0(data_for_direct$A,'.',c('rural','urban')[data_for_direct$U+1])
  data_for_direct$age <- '0'
  
  dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                               regionVar = "A",timeVar = "years", clusterVar =  "~cluster",
                               Ntrials = "n",weightsVar = "wt",national.only = T)
}

## run Alg 5.2 model once
{
  data_list5.2 <- data_list
  data_list5.2$Yplus <- NULL
  data_list5.2$obs_dat <- data_list$obs_dat %>% dplyr::select(-c(Y))
  
  #run model
  start.time <- Sys.time()
  postsamp_5.2 <- postsamp_Alg5.2_MALA_MCMC(eps=0.65,prop_sd_logd = 0.15,
                                            data_list = data_list5.2, 
                                            logit_r_hat=dir.est$logit.est, 
                                            logit_r_hat_var = dir.est$var.est,
                                            n_iter=3000,chains = 4, burnin=1000)
  Sys.time() - start.time
  
  postsamp_5.2$acc_reg_rate
  postsamp_5.2$acc_d_rate
  sum(is.na(postsamp_5.2$reg_params[,1]))
  
  # test using average cluster size
  data_list5.2a <- data_list
  data_list5.2a$Yplus <- NULL
  #changing cluster sizes to be estimates
  data_list5.2a$obs_dat <- merge(data_list5.2a$obs_dat,data_list5.2a$obs_dat %>% group_by(A,U) %>% 
                                   summarise(N = round(sum(N)/n())),by=(c('A','U')))
  data_list5.2a$obs_dat <- data_list5.2a$obs_dat %>% select(-c(Y,N.x)) %>% rename(N=N.y)
  data_list5.2a$obs_dat <- data_list5.2a$obs_dat[order(data_list5.2a$obs_dat$A),]
  
  start.time <- Sys.time()
  postsamp_5.2a <- postsamp_Alg5.2_MALA_MCMC(eps=0.6,  k_Yplus = 500,
                                             data_list = data_list5.2a, 
                                             logit_r_hat=dir.est$logit.est, 
                                             #logit_r_hat=logitlink(data_list5.2$Yplus/sum(pop_dat$N)), # to compare with Alg 5.2
                                             logit_r_hat_var = dir.est$var.est,
                                             n_iter=4000,chains = 5, burnin=1000)
  Sys.time() - start.time
  
}

pdf("Handcoded MCMC Simulations/Alg 5.2 MALA 240117, 5 chains of 5k iterations posterior dist.pdf")
{
  hist(log(postsamp_5.2$tau))
  abline(v=log(tau),col='red')
  hist(log(postsamp_5.2$d))
  abline(v=log(d),col='red')
  hist(postsamp_5.2$alphaU)
  abline(v=alphaU,col='red')
  hist(postsamp_5.2$alphaR)
  abline(v=alphaR,col='red')
  
  results_r <- gather(data.frame(postsamp_5.2$r),'strata','r')
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
  
  plot(log(postsamp_5.2$d),type='l')
  abline(h=log(d),col='red')
  plot(postsamp_5.2$reg_params[,1],type='l')
  abline(h=alphaU,col='red')
  plot(postsamp_5.2$reg_params[,2],type='l')
  abline(h=alphaR,col='red')
  
  par(mfrow=c(3,2))
  for(k in 3:(n_areas+2)){
    plot(postsamp_5.2$reg_params[,k],type='l',main=k-2)
    abline(h=b[k-2],col='red')
  }
  
  each_chain <- 2000
  plot(cumsum(postsamp_5.2$Yplus[1:each_chain])/(1:each_chain),type='l',col=1,ylab='',
       ylim=c(min(postsamp_5.2$Yplus),max(postsamp_5.2$Yplus)),main = paste0('Running mean of Yplus'))
  for(c in 2:chains){
    lines(cumsum(postsamp_5.2$Yplus[((c-1)*each_chain+1):(c*each_chain)])/(1:each_chain),col=c)
  }
  
}
dev.off()

pdf("Handcoded MCMC Simulations/Compare Alg 5.2 w known v estimated cluster sizes (SD=10,12) 240118.pdf")
{
  # code to compare 2 runs (with same data)
  plot(colMeans(postsamp_5.2$r),colMeans(postsamp_5.2v2$r,na.rm = T))
  abline(0,1)
  plot((colMeans(postsamp_5.2$r) - colMeans(postsamp_5.2v2$r,na.rm = T))/colMeans(postsamp_5.2$r),ylab='Percentage difference in rate estimate')
  
  par(mfrow=c(3,2))
  for(i in 1:ncol(postsamp_5.2$r)){
    hist(postsamp_5.2$r[,i],probability = T,col='lightblue')
    hist(postsamp_5.2a$r[,i],add=T,probability = T,col=rgb(1,0,0,0.25))
  }
}
dev.off()

#compare to Negbin INLA model
{
  prior.hyper <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01)))
  prior.fixed <- list(mean =list("factor(U)0"=-3, "factor(U)1"=-3.5, default = 0), prec=list(default = 1/9))
  inla.fit <- INLA::inla(Z ~ factor(U) -1 +
                           f(A, model = "iid",hyper=prior.hyper),
                         data=data_list$obs_dat, family='nbinomial', E=n,
                         control.predictor = list(compute = F, link = 1),
                         control.fixed = prior.fixed,
                         control.family = list(link='log'),control.compute = list(config=T))
  # sample from posterior
  cs <- inla.fit$misc$configs$contents$tag
  cs <- cs[cs != "Predictor"]
  select <- list()
  for (i in 1:length(cs)) {
    select[[i]] <- 0
    names(select)[i] <- cs[i]
  }
  
  sampFull <- INLA::inla.posterior.sample(n = 1000, result = inla.fit, intern = TRUE, selection = select)
  sampFull.draws <- matrix(NA,1000,length(sampFull[[1]]$latent))
  for(i in 1:1000){
    sampFull.draws[i,] <- sampFull[[i]]$latent
  }
  fields <- colnames(sampFull.draws) <- row.names(sampFull[[1]]$latent)
  
  inla.eta.res <- cbind(sampFull.draws[,n_areas+2] + sampFull.draws[,1:n_areas],sampFull.draws[,n_areas+1] + sampFull.draws[,1:n_areas])
  inla.r.res <- exp(inla.eta.res)
  
  pdf('Compare Alg 5.2 to INLA 240103, simulated dataset.pdf')
  {
    
    #compare urban and rural intercepts
    hist(postsamp_5.2$alphaU,probability = T,col='lightblue',main='Urban intercept',border=F)
    hist((sampFull.draws[,22]),add=T,probability = T,col=rgb(0,0,0,0),border='grey50')
    abline(v=alphaU,col='red')
    
    hist(postsamp_5.2$alphaR,probability = T,col='lightblue',main='Rural intercept',border=F)
    hist((sampFull.draws[,21]),add=T,probability = T,col=rgb(0,0,0,0),border='grey50')
    abline(v=alphaR,col='red')
    
    #compare random effect estimates
    par(mfrow=c(3,2))
    for(i in 1:(n_areas)){
      hist(postsamp_5.2$b[,i],probability = T,col='lightblue',main=paste0('Area ',i),border = F)
      hist((sampFull.draws[,i]),add=T,probability = T,col=rgb(0,0,0,0),border='grey50')
      abline(v=b[i],col='red')
    }
    par(mfrow=c(2,1))
    
    #compare rates
    plot(1:20,colmeans(postsamp_5.2$r)[1:20],col='blue',pch=20,ylim=c(0.015,0.06),main='Urban NMR estimate for each area')
    points(1:20,colmeans(inla.r.res)[1:20],col='grey50',pch=20)
    points(1:20,exp(alphaU+b),col='red',pch=3)
    abline(h=natl_r_strat$r[2],col='red',lty=2)
    legend('topleft',legend = c('Alg 5.2','INLA','truth'),col=c('blue','grey50','red'),pch=c(20,20,3),cex=0.75)
    legend('topright',legend = c('true National rate'),col=c('red'),lty=3,cex=0.5)
    
    plot(21:40,colmeans(postsamp_5.2$r)[21:40],col='blue',pch=20,ylim=c(0.015,0.07),main='Rural NMR estimate for each area')
    points(21:40,colmeans(inla.r.res)[21:40],col='grey50',pch=20)
    points(21:40,exp(alphaR+b),col='red',pch=3)
    abline(h=natl_r_strat$r[1],col='red',lty=2)
    legend('bottomleft',legend = c('Alg 5.2','INLA','truth'),col=c('blue','grey50','red'),pch=c(20,20,3),cex=0.75)
    legend('bottomright',legend = c('true National rate'),col=c('red'),lty=3,cex=0.5)
    
    par(mfrow=c(3,2))
    for(i in 1:(2*n_areas)){
      hist(postsamp_5.2$r[,i],probability = T,col='lightblue',border = F,main=paste0('Stratum ',i),xlim=c(0.015,0.1))
      hist(inla.r.res[,i],add=T,probability = T,col=rgb(0,0,0,0),border='grey50')
      abline(v=exp(c(alphaU+b,alphaR + b)[i]),col='red')
    }
  }
  dev.off()
}

#######################################################################
#######################################################################
# Alg 5.2a (Conditioning on Natl, BYM2): P(alphaU,alphaR,b,phi,tau,d,Y^s,Y+|r+_est,Z,delta) ------

# tau_init=40
# phi_init = 0.5
# d_init = seq(0.8,1.2,by=0.1) 
# alphaU_init = seq(-4,-3,by=0.1)
# alphaR_init= seq(-3.5,-2.5,by=0.1)  #initial values
# U_tau = 0.5
# alpha_tau = 2/3
# mu_U = -3.5
# mu_R = -3
# sd_U<- sd_R <- 9
# mu_d = 0
# sd_d = 1 #hyperpriors
# 
# eps=0.65
# prop_sd_tau <- 0.3
# prop_sd_phi <- 0.3
# logit_r_hat=dir.est$logit.est
# logit_r_hat_var = dir.est$var.est
# n_iter=4000
# chains = 4
# burnin=1000

#preconditioned MALA as described in "Adaptive step size selection for Hessian-based manifold Langevian samplers" Tore Selland Kleppe
postsamp_Alg5.2a_MALA_MCMC <- function(tau_init=40, phi_init=0.5, d_init = seq(0.8,1.2,by=0.1), alphaU_init = seq(-4,-3,by=0.1), alphaR_init= seq(-3.5,-2.5,by=0.1), #initial values
                                       U_tau = 0.5,alpha_tau = 2/3,mu_U = -3.5,sd_U = 3,mu_R = -3.2,sd_R = 3,mu_d = 0,sd_d = 1,#hyperpriors  
                                       eps, prop_sd_tau, prop_sd_phi, #tuning params for regression param and Yplus proposals
                                       logit_r_hat,logit_r_hat_var,data_list, 
                                       n_iter, burnin=n_iter/2,chains=5){
  
  data <- data_list$obs_dat
  data <- data[order(data$A),]
  pop_dat <- data_list$pop_strata_dat
  Q_scaled_inv <- data_list$Q_scaled_inv
  n_areas <- length(unique(data_list$pop_strata_dat$A))
  
  parallel_results <- mclapply(1:chains, function(j){
    
    b_init <- rnorm(n_areas,0,1/sqrt(tau_init))
    Yplus_init <- round(sum(pop_dat$N)*expit(logit_r_hat))
    W_inv_init <- solve(W(phi_init,Q_scaled_inv))
    
    ## ensure we are starting with a PD negative hessian 
    reg_params_init_grid <- as.matrix(expand.grid(log(d_init),alphaU_init,alphaR_init))
    search_order <- sample(1:nrow(reg_params_init_grid),nrow(reg_params_init_grid))
    
    for(j in 1:nrow(reg_params_init_grid)){
      reg_params_init <-  c(reg_params_init_grid[search_order[j],],b_init)
      Nr_sum_init <- sum(pop_dat$N*exp(pop_dat$U*reg_params_init[2] + (1-pop_dat$U)*reg_params_init[3] + reg_params_init[3+pop_dat$A]))
      Nr_obs_init <- data$N*exp(data$U*reg_params_init[2] + (1-data$U)*reg_params_init[3] + reg_params_init[3+data$A])
      Y_init <- rmultinom(1,Yplus_init,prob=c(Nr_obs_init,Nr_sum_init-sum(Nr_obs_init)))
      Y_init <- Y_init[-length(Y_init)]
      Y_init[Y_init<data$Z] <- data_list$obs_dat$Z[Y_init<data$Z]
      
      hess_init <- pracma::hessian(log_target_5.2a_reg_params,x=reg_params_init,data=data,Yplus = Yplus_init,pop_dat=pop_dat,
                                   Y=Y_init,tau=tau_init,W_inv=W_inv_init,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      if(is.positive.definite(-hess_init)){
        break
      }
    }
    
    ## initialize data frames 
    tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
    phi_postsamp_t <- c(phi_init, rep(NA,n_iter-1))
    reg_params_postsamp_t <- rbind(reg_params_init, matrix(NA,n_iter-1,length(reg_params_init)))
    Y_postsamp_t <- rbind(Y_init, matrix(NA,n_iter-1,nrow(data)))
    Yplus_postsamp_t <- c(Yplus_init, rep(NA,n_iter-1))
    acc_reg <- acc_hyper <- 0
    
    # set current states
    tau_current <- tau_init
    phi_current <- phi_init
    d_current <- exp(reg_params_init[1])
    alphaU_current <- reg_params_init[2]
    alphaR_current <- reg_params_init[3]
    b_current <- b_init
    reg_params_current <- reg_params_init
    Y_current <- Y_init
    Yplus_current <- Yplus_init
    W_inv_current <- W_inv_init
    
    #start.time <- Sys.time()
    for(i in 2:n_iter){
      
      # draw new tau and phi----
      
      log_tau_proposal <- rnorm(1,log(tau_current),prop_sd_tau)
      logit_phi_proposal <- rnorm(1,logitlink(phi_current),prop_sd_phi)
      W_inv_proposal <- solve(W(expit(logit_phi_proposal),Q_scaled_inv))
      
      a_prob_reg <- exp({
        #likelihood
        (n_areas-1)/2*log_tau_proposal - 0.5*exp(log_tau_proposal)*t(b_current)%*%W_inv_proposal%*%b_current + log(alpha_tau)/U_tau*exp(-0.5*log_tau_proposal) -
          ((n_areas-1)/2*log(tau_current) - 0.5*exp(log(tau_current))*t(b_current)%*%W_inv_current%*%b_current + log(alpha_tau)/U_tau*exp(-0.5*log(tau_current))) +
          #proposal dist (jacobian for logit phi)
          log(expit(logit_phi_proposal)*(1-expit(logit_phi_proposal))) - log(phi_current*(1-phi_current))
      })
      
      if(a_prob_reg>runif(1)){
        tau_current <- exp(log_tau_proposal)
        phi_current <- expit(logit_phi_proposal)
        W_inv_current <- solve(W(phi_current,Q_scaled_inv))
        if(i>burnin)
          acc_hyper <- acc_hyper + 1
      }
      
      # draw proposal for regression params ----
      grad_current <- pracma::grad(log_target_5.2a_reg_params,x=reg_params_current,data=data,Yplus = Yplus_current, pop_dat=pop_dat,
                                   Y=Y_current,tau=tau_current,W_inv=W_inv_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      
      hess_current <- pracma::hessian(log_target_5.2a_reg_params,x=reg_params_current,data=data,Yplus = Yplus_current,pop_dat=pop_dat,
                                      Y=Y_current,tau=tau_current,W_inv=W_inv_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      if(is.positive.definite(-hess_current)){
        cov_current <- matlib::inv(-hess_current)
      }else{
        message(paste0('current Hessian is not PD at iteration ',i))
        hess_chol <- get.GMW.chol(-hess_current,0.01)
        cov_current <- matlib::inv(hess_chol%*%t(hess_chol))
      }
      L_current <- t(chol(cov_current))
      
      r <- rnorm(length(reg_params_current),0,1)
      reg_params_proposal <- reg_params_current + eps^2/2*L_current%*%t(L_current)%*%grad_current + eps*L_current%*%r
      
      grad_proposal <- pracma::grad(log_target_5.2a_reg_params,x=reg_params_proposal,data=data,Yplus = Yplus_current,pop_dat=pop_dat,
                                    Y=Y_current,tau=tau_current,W_inv=W_inv_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      hess_proposal <- pracma::hessian(log_target_5.2a_reg_params,x=reg_params_proposal,data=data,Yplus = Yplus_current,pop_dat=pop_dat,
                                       Y=Y_current,tau=tau_current,W_inv=W_inv_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      if(is.symmetric(hess_proposal)){
        if(is.positive.definite(-hess_proposal)){
          cov_proposal <- matlib::inv(-hess_proposal)
          if(is.positive.definite(cov_proposal)){
            skip.it <- F
          }else{
            message(paste0('Proposed Hessian not PD at iteration ',i,'\n'))
            skip.it <- T
          }
        }else{
          message(paste0('Proposed Hessian not PD at iteration ',i,'\n'))
          skip.it <- T}
      }else{
        message(paste0('Proposed Hessian not symmetric at iteration ',i,'\n'))
        skip.it <- T}
      
      
      if(!skip.it){
        # accept or reject proposal for (alpha, b, d)
        L_proposal <- t(chol(cov_proposal))
        
        a_prob_reg <- min(1,exp({
          log_target_5.2a_reg_params(reg_params_proposal,data=data,Yplus = Yplus_current,pop_dat=pop_dat,Y=Y_current,tau=tau_current,W_inv=W_inv_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R) +
            (-0.5*sum((solve(L_proposal)%*%(reg_params_current-reg_params_proposal)/eps - 0.5*eps*t(L_proposal)%*%grad_proposal)^2)) - 
            log_target_5.2a_reg_params(reg_params_current,data=data,Yplus = Yplus_current,pop_dat=pop_dat,Y=Y_current,tau=tau_current,W_inv=W_inv_current,mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R) -
            (-0.5*sum(r^2))
        }))
        
        if(a_prob_reg>runif(1)){
          reg_params_current <- reg_params_proposal
          d_current <- exp(reg_params_proposal[1])
          alphaU_current <- reg_params_proposal[2]
          alphaR_current <- reg_params_proposal[3]
          b_current <- reg_params_proposal[4:length(reg_params_proposal)]
          if(i>burnin){acc_reg <- acc_reg + 1}
        }
      }
      
      Nr_sum_current <- sum(pop_dat$N*exp(pop_dat$U*alphaU_current + (1-pop_dat$U)*alphaR_current + b_current[pop_dat$A]))
      Nr_obs_current <- data$N*exp(data$U*alphaU_current + (1-data$U)*alphaR_current + b_current[data$A])
      # draw proposal for Yplus ----
      
      domain <- round(0.5*Yplus_current):(2*Yplus_current)
      f <- log_target_cond_5.2_Yplus(domain,pop_dat,Y_current,Nr_sum_current,Nr_obs_current,
                                     d_current,logit_r_hat,logit_r_hat_var)
      
      del_f <- f-mean(f) #normalize
      #restrict domain to densities >0
      zero_dens <- which(exp(del_f)==0)
      while(length(zero_dens)>0){
        domain <- domain[-zero_dens]
        f <- f[-zero_dens]
        del_f <- f-mean(f)
        zero_dens <- which(exp(del_f)==0)
      }
      
      target_density <- (exp(del_f)/sum(exp(del_f)))
      target_F <- cumsum(target_density)
      
      #take draw from inverse CDF
      U <- runif(1)
      Yplus_current <- (min(domain[target_F>U]))
      
      # draw proposal for Ys ----
      for(k in 1:nrow(data)){
        # derive CDF
        domain <- data$Z[k]:round(0.5*data$N[k])
        f <- log_target_cond_5.2_y(domain,data,Yplus_current,pop_dat,
                                   Y_current,k,Nr_sum_current,Nr_obs_current,d_current)
        del_f <- f-mean(f) #normalize
        #restrict domain to densities >0
        zero_dens <- which(exp(del_f)==0)
        while(length(zero_dens)>0){
          domain <- domain[-zero_dens]
          f <- f[-zero_dens]
          del_f <- f-mean(f)
          zero_dens <- which(exp(del_f)==0)
        }
        
        target_density <- (exp(del_f)/sum(exp(del_f)))
        target_F <- cumsum(target_density)
        #take draw from inverse CDF
        U <- runif(1)
        Y_current[k] <- min(domain[target_F>U])
      }
      
      
      # record current state ----
      tau_postsamp_t[i] <- tau_current
      phi_postsamp_t[i] <- phi_current
      if(!skip.it){
        reg_params_postsamp_t[i,] <- reg_params_current
      }
      Y_postsamp_t[i,] <- Y_current
      Yplus_postsamp_t[i] <- Yplus_current
    } 
    #Sys.time() - start.time
    
    results <- list(reg_params_t=reg_params_postsamp_t[(burnin+1):n_iter,], 
                    tau_t=tau_postsamp_t[(burnin+1):n_iter], 
                    phi_t=phi_postsamp_t[(burnin+1):n_iter], 
                    Y_t=Y_postsamp_t[(burnin+1):n_iter,],
                    Yplus_t=Yplus_postsamp_t[(burnin+1):n_iter], 
                    acc_reg_t = acc_reg,acc_hyper_t = acc_hyper)
    
    return(results)
  },mc.cores = chains)
  
  #combine chains
  b_postsamp <- reg_params_postsamp <- Y_postsamp <- NULL
  for(j in 1:chains){
    b_postsamp <- rbind(b_postsamp,parallel_results[[j]]$reg_params_t[,4:(3+n_areas)])
    reg_params_postsamp <- rbind(reg_params_postsamp,parallel_results[[j]]$reg_params_t)
    Y_postsamp <- rbind(Y_postsamp,parallel_results[[j]]$Y_t)
  }
  tau_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$tau_t}))
  phi_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$phi_t}))
  d_postsamp <- unlist(lapply(1:chains,function(x){exp(parallel_results[[x]]$reg_params_t[,1])}))
  alphaU_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,2]}))
  alphaR_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,3]}))
  Yplus_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$Yplus_t}))
  names(tau_postsamp) <-names(phi_postsamp) <- names(d_postsamp) <- names(alphaU_postsamp) <- names(alphaR_postsamp) <- NULL
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp,alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  
  acc_reg_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_reg_t})))/((n_iter-burnin)*chains)
  acc_hyper_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_hyper_t})))/((n_iter-burnin)*chains)
  
  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              tau=tau_postsamp, phi=phi_postsamp, 
              reg_params=reg_params_postsamp,
              d=d_postsamp, 
              Y=Y_postsamp,
              Yplus=Yplus_postsamp,
              acc_reg_rate = acc_reg_rate,
              acc_hyper_rate = acc_hyper_rate))
}

#######################################################################
#######################################################################
# Applying Alg 5.2/5.2a to real DHS data (kind of a mess bc Nigeria was not good example) ------
country_t <- 'Nigeria'
year_range <- 2014:2018
survey_year <- 2018

load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/",country_t,"_cluster_dat.rda"))
mod.dat <- as.data.frame(nmr.dat)
mod.dat$years <- max(year_range)
mod.dat$survey <- survey_year
mod.dat$age <- 0
# stick to definition of urban/rural at time of sampling frame
mod.dat[mod.dat$admin1.name=='Lagos',]$urban <- 'urban'
mod.dat[mod.dat$admin1.name=='Lagos',]$strata <- 'Lagos south west.urban'

load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/",country_t,"_Amat.rda"))

load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/worldpop/adm1_weights_u1.rda"))
admin_key <- readxl::read_excel("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Admin1_Admin2_Key.xlsx") %>% filter(country==country_t) %>% dplyr::select(-country)
admin_key$A1 <- as.numeric(sapply(1:nrow(admin_key),function(i){str_split(admin_key$admin1.char[i],'_')[[1]][2]}))
admin_key <- admin_key[order(admin_key$A1),c('A1','admin1.char','admin1.name')]

# prepare direct estimates
{
  data_for_direct <- mod.dat %>% filter(years %in% year_range, survey==survey_year,age==0)
  data_for_direct$years <- 'All'
  data_for_direct$died <- data_for_direct$Y
  data_for_direct$v005 <- data_for_direct$v005/1e6
  dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                               regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                               Ntrials = "total",weightsVar = "v005",national.only = T)
}

# input strata level DHS survey info 
{
  dhs_dat <- NULL
  total_births <- data.frame(year=year_range, births = c(7189000,7273000,7384000,7487000,7590000)) # total births in each year (UN estimate)
  for(year_t in year_range){
    adm1_UR_weights_tmp <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_NGA/admin1_u1_",year_t,"_urban_frac.rds"))
    #subset year
    adj_pop <- weight.adm1.u1[,c(1,4,5)] %>% filter(year==year_t)
    #rescale number of births in each area based on national total
    adj_pop$admin_pop_adj <- adj_pop$proportion*total_births[total_births$year==year_t,]$births
    pop_tmp <- merge(adj_pop,adm1_UR_weights_tmp[,c(2,3)],by.x='admin1.char',by.y='adm_idx')
    pop_tmp <- merge(expand_grid(admin1.char = unique(admin_key$admin1.char),U=c(0,1)),pop_tmp)
    pop_tmp$frac <- ifelse(pop_tmp$U==1,pop_tmp$urb_frac,1-pop_tmp$urb_frac)
    pop_tmp$N <- pop_tmp$admin_pop_adj*pop_tmp$frac # assuming fertility per population is same across urban and rural
    dhs_dat <- rbind(dhs_dat,pop_tmp[,c(1,2,4,8)])
  }
  
  pop_strata_dat <- dhs_dat %>% group_by(admin1.char,U) %>% summarise(N=sum(N))
  pop_strata_dat <- merge(admin_key,pop_strata_dat)
  pop_strata_dat <- pop_strata_dat[order(pop_strata_dat$A1,1-pop_strata_dat$U),]
  
  # total number of clusters in strata (from DHS)
  pop_strata_dat$n_clusters <- c(2106,9463,2820,9988,908,16205,18409,3498,2761,17124,2628,6379,2006,20850,7798,16288,1410,14912,9008,9201,
                                 11911,1977,7964,4829,9438,2123,9774,4223,2452,1138,1955,7539,10006,9567,2293,18900,9529,12263,16957,19402,
                                 6874,26442,2621,14020,5492,10354,11715,4556,25424,0,2008,7211,5126,18319,7085,7408,8588,10667,19810,6097,
                                 22405,8701,3949,11930,12480,12381,2548,10231,1657,8943,3053,11870,3090,13942)
  #pop_strata_dat$n_HH_samp <- 30 #number of HH sampled in each cluster
  pop_strata_dat$N_bar <- pop_strata_dat$N/pop_strata_dat$n_clusters #average number of births per cluster
}

# organize survey data
{
  alg_dat <- mod.dat %>% filter(years %in% year_range, survey==survey_year,age==0) %>% mutate(U=ifelse(urban=='urban',1,0)) %>% 
    dplyr::select(cluster,total,Y,U,admin1,admin1.char,admin1.name) %>% rename(Z=Y,n=total,A1=admin1)
  alg_dat <- merge(alg_dat,admin_key)
}

# estimate average number of births in SAMPLED cluster
{
  frame.info <- merge(pop_strata_dat,alg_dat %>% group_by(admin1.char,U) %>% summarise(n_clusters_samp = n())) %>% 
    dplyr::select(A1,admin1.char,U,N_bar,n_clusters,n_clusters_samp)
  frame.info$sd_N <- 0.75*frame.info$N_bar
  
  n_iter <- 100
  N_obs_bar <- rep(NA,n_iter)
  frame.info$N <- NA
  
  for(i in 1:nrow(frame.info)){
    for(k in 1:n_iter){
      N = round(rnorm(frame.info$n_clusters[i],frame.info$N_bar[i],frame.info$sd_N[i]))
      N[N<0] <- 0
      N_obs <- N[sample(1:frame.info$n_clusters[i],frame.info$n_clusters_samp[i],prob = N)]
      N_obs_bar[k] <- mean(N_obs)
    }
    frame.info$N[i] <- round(mean(N_obs_bar))
  }
}

# add estimated number of total births per cluster
alg_dat <- merge(alg_dat,frame.info) %>% dplyr::select(-N_bar,-n_clusters,-n_clusters_samp) 

# organize population and survey info
pop_strata_dat <- pop_strata_dat[,c('A1','U','N')]
pop_strata_dat$N <- round(pop_strata_dat$N)

# get generalized inverse of scaled Q for BYM2
Q_scaled <- INLA::inla.scale.model(admin1.mat, 
                                   constr=list(A=t(eigen(admin1.mat)$vectors[,eigen(admin1.mat)$values<1e-10]), e=rep(0,sum(eigen(admin1.mat)$values<1e-10))))
Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))

data_list <- list(obs_dat = alg_dat,
                  pop_strata_dat = pop_strata_dat,
                  Q_scaled_inv = Q_scaled_inv)


start.time <- Sys.time()
postsamp_5.2 <- postsamp_Alg5.2_MALA_MCMC(eps=0.15, prop_sd_logd = 0.075,
                                          data_list = data_list, 
                                          logit_r_hat=dir.est$logit.est, 
                                          logit_r_hat_var = dir.est$var.est,
                                          n_iter=1000,chains = 4,burnin=00)
Sys.time() - start.time

postsamp_5.2$acc_reg_rate
postsamp_5.2$acc_d_rate
sum(is.na(postsamp_5.2$alphaU))
plot(postsamp_5.2$alphaR,type='l')

#######################################################################
#######################################################################
# Alg 5.3 (Conditioning on Admin1, IID): P(alphaU,alphaR,b,tau,d,Y^s,Y+|r+_est,Z,delta) ------

# tau_init=40
# d_init = 1
# alphaU_init = seq(-4,-3,by=0.25)
# alphaR_init= seq(-3.5,-2.5,by=0.25)  #initial values
# a_tau <- b_tau <- 0.01
# mu_U = -3.5
# mu_R = -3
# sd_U<- sd_R <- 9
# mu_logd = 0
# sd_logd = 1 #hyperpriors
# 
# eps=0.65
# prop_sd_logd=0.3
# data_list = data_list5.3
# logit_r_hat=dir.est$logit.est
# logit_r_hat_var = dir.est$var.est
# n_iter=100
# burnin=50
# chains = 4

#preconditioned MALA as described in "Adaptive step size selection for Hessian-based manifold Langevian samplers" Tore Selland Kleppe
postsamp_Alg5.3_MALA_MCMC <- function(tau_init=40, d_init = 1, alphaU_init = seq(-4,-3,by=0.25), alphaR_init= seq(-3.5,-2.5,by=0.25), #initial values
                                      a_tau = 0.01,b_tau = 0.01,mu_U = -3.5,sd_U = 3,mu_R = -3.2,sd_R = 3,mu_logd = 0,sd_logd = 1,#hyperpriors  
                                      eps, prop_sd_logd, #tuning params for regression param and d
                                      logit_r_hat,logit_r_hat_var,data_list, 
                                      n_iter, burnin=n_iter/2,chains=4){
  
  data <- data_list$obs_dat
  data <- data[order(data$A2),]
  pop_dat <- data_list$pop_strata_dat
  #make data tables
  data_table <- as.data.table(data)
  pop_dt <- as.data.table(pop_dat)
  
  admin1_pop_dat <- pop_dat %>% group_by(A1) %>% summarise(N = sum(N))
  n_admin1 <- length(unique(data_list$pop_strata_dat$A1))
  n_admin2 <- length(unique(data_list$pop_strata_dat$A2))
  
  parallel_results <- mclapply(1:chains, function(j){
    
    b_init <- rnorm(n_admin2,0,1/sqrt(tau_init))
    Yplus_init <- round(admin1_pop_dat$N*expit(logit_r_hat))
    
    ## ensure we are starting with a PD negative hessian 
    reg_params_init_grid <- as.matrix(expand.grid(alphaU_init,alphaR_init))
    search_order <- sample(1:nrow(reg_params_init_grid),nrow(reg_params_init_grid))
    
    for(j in 1:nrow(reg_params_init_grid)){
      reg_params_init <-  c(reg_params_init_grid[search_order[j],],b_init)
      
      Y_init <- sapply(1:n_admin1,function(i){
        subdat <- data %>% filter(A1==i)
        subpop_dat <- pop_dat %>% filter(A1==i)
        
        Nr_obs_init <- subdat$N*exp(subdat$U*reg_params_init[1] + (1-subdat$U)*reg_params_init[2] + reg_params_init[2+subdat$A2])
        Nr_sum_init <- sum(subpop_dat$N*exp(subpop_dat$U*reg_params_init[1] + (1-subpop_dat$U)*reg_params_init[2] + reg_params_init[2+subpop_dat$A2]))
        
        Y_init_t <- rmultinom(1,Yplus_init[i],prob=c(Nr_obs_init,Nr_sum_init-sum(Nr_obs_init)))
        Y_init_t <- Y_init_t[-length(Y_init_t)]
        Y_init_t[Y_init_t<subdat$Z] <- subdat$Z[Y_init_t<subdat$Z]
        
        return(Y_init_t)
      })
      data_table[, ':='(Y=unlist(Y_init))]
      Y_sum_obs_init <- (data_table[, sum(Y),by=A1])$V1
      
      hess_init <- pracma::hessian(log_target_5.3_reg_params,x=reg_params_init,data_table=data_table,Y_sum_obs=Y_sum_obs_init,
                                   Yplus = Yplus_init,pop_dt=pop_dt,d=d_init,
                                   tau=tau_init,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      
      if(is.positive.definite(-hess_init)){
        L_init <- t(chol(-hess_init))
        break
      }
      print(j)
    }
    
    ## initialize data frames 
    tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
    d_postsamp_t <- c(d_init, rep(NA,n_iter-1))
    reg_params_postsamp_t <- rbind(reg_params_init, matrix(NA,n_iter-1,length(reg_params_init)))
    Yplus_postsamp_t <- rbind(Yplus_init, matrix(NA,n_iter-1,n_admin1))
    acc_reg <- acc_d <- 0
    
    # set current states
    tau_current <- tau_init
    d_current <- d_init
    alphaU_current <- reg_params_init[1]
    alphaR_current <- reg_params_init[2]
    b_current <- b_init
    reg_params_current <- reg_params_init
    Y_current <- Y_init
    Yplus_current <- Yplus_init
    Y_sum_obs_current <- Y_sum_obs_init
    
    for(i in 2:n_iter){
      # draw new tau from full conditional ----
      tau_current <- rgamma(1,n_admin2/2+a_tau,0.5*sum(b_current^2)+b_tau)
      
      # draw proposal for regression params ----
      grad_current <- pracma::grad(log_target_5.3_reg_params,x=reg_params_current,data_table=data_table,d=d_current,Yplus = Yplus_current, pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,
                                   tau=tau_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      hess_current <- pracma::hessian(log_target_5.3_reg_params,x=reg_params_current,data_table=data_table,d=d_current,Yplus = Yplus_current,pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,
                                      tau=tau_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      if(is.positive.definite(-hess_current)){
        cov_current <- matlib::inv(-hess_current)
      }else{
        message(paste0('current Hessian is not PD at iteration ',i))
        hess_chol <- get.GMW.chol(-hess_current,0.01)
        cov_current <- matlib::inv(hess_chol%*%t(hess_chol))
      }
      L_current <- t(chol(cov_current))
      
      r <- rnorm(length(reg_params_current),0,1)
      reg_params_proposal <- reg_params_current + eps^2/2*L_current%*%t(L_current)%*%grad_current + eps*L_current%*%r
      
      grad_proposal <- pracma::grad(log_target_5.3_reg_params,x=reg_params_proposal,data_table=data_table,d=d_current,Yplus = Yplus_current,pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,
                                    tau=tau_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      hess_proposal <- pracma::hessian(log_target_5.3_reg_params,x=reg_params_proposal,data_table=data_table,d=d_current,Yplus = Yplus_current,pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,
                                       tau=tau_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      #use taylors approximation -- not working
      # y <- grad_proposal - grad_current
      # s <- reg_params_proposal - reg_params_current
      # B_mat <- diag(1,length(reg_params_current))
      # 
      # B <- optim(par=B_mat[lower.tri(B_mat,diag=T)],fn=function(B){
      #   B_mat[lower.tri(B_mat,diag=T)] <- B
      #   sum((y+t(s)%*%B_mat%*%t(B_mat))^2)},method = 'BFGS')
      # 
      # B_mat[lower.tri(B_mat,diag=T)] <- B$par
      # L_proposal <- t(solve(B_mat)) #LLt=cov_proposal
      
      if(is.symmetric(hess_proposal)){
        if(is.positive.definite(-hess_proposal)){
          cov_proposal <- matlib::inv(-hess_proposal)
          if(is.positive.definite(cov_proposal)){
            skip.it <- F
          }else{
            message(paste0('Proposed Hessian not PD at iteration ',i,'\n'))
            skip.it <- T
          }
        }else{
          message(paste0('Proposed Hessian not PD at iteration ',i,'\n'))
          skip.it <- T}
      }else{
        message(paste0('Proposed Hessian not symmetric at iteration ',i,'\n'))
        skip.it <- T}
      
      if(!skip.it){
        # accept or reject proposal for (alpha, b, d)
        L_proposal <- t(chol(cov_proposal))
        
        a_prob_reg <- min(1,exp({
          log_target_5.3_reg_params(reg_params_proposal,data_table=data_table,Yplus = Yplus_current,pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,d=d_current,tau=tau_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R) +
            (-0.5*sum((solve(L_proposal)%*%(reg_params_current-reg_params_proposal)/eps - 0.5*eps*t(L_proposal)%*%grad_proposal)^2)) - #optimize
            log_target_5.3_reg_params(reg_params_current,data_table=data_table,Yplus = Yplus_current,pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,d=d_current,tau=tau_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R) -
            (-0.5*sum(r^2)) # proposed given current
        }))
        
        if(a_prob_reg>runif(1)){
          reg_params_current <- reg_params_proposal
          alphaU_current <- reg_params_proposal[1]
          alphaR_current <- reg_params_proposal[2]
          b_current <- reg_params_proposal[3:length(reg_params_proposal)]
          if(i>burnin){acc_reg <- acc_reg + 1}
        }
      }
      
      Nr_sum_current <- sapply(1:n_admin1,function(i){
        subpop_dat <- pop_dat %>% filter(A1==i)
        Nr_sum_current_t <- sum(subpop_dat$N*exp(subpop_dat$U*alphaU_current + (1-subpop_dat$U)*alphaR_current + b_current[subpop_dat$A2]))
        return(Nr_sum_current_t)
      })
      
      Nr_obs_current <- sapply(1:n_admin1,function(i){
        subdat <- data %>% filter(A1==i)
        Nr_obs_t <- subdat$N*exp(subdat$U*alphaU_current + (1-subdat$U)*alphaR_current + b_current[subdat$A2])
        return(Nr_obs_t)
      })
      
      # draw proposal for d -------
      logd_proposal <- rnorm(1,log(d_current),prop_sd_logd)
      
      a_prob_d <- exp({log_target_5.3_logd(logd_proposal,alphaU_current,alphaR_current,b_current,data,
                                           Yplus_current,pop_dat,Y_current,mu_logd,sd_logd) -
          log_target_5.3_logd(log(d_current),alphaU_current,alphaR_current,b_current,data,
                              Yplus_current,pop_dat,Y_current,mu_logd,sd_logd)})
      
      if(a_prob_d > runif(1)){
        d_current <- exp(logd_proposal)
        if(i>burnin){
          acc_d <- acc_d + 1
        }
      }
      
      # draw for Yplus's from conditional distribution (ITS) ----
      
      Yplus_current <- sapply(1:n_admin1,function(k){
        # calculate Poisson parameter for domain centered on current Y
        domain <- round(0.5*Yplus_current[k]):(2*Yplus_current[k])
        f <- log_target_cond_5.3_Yplus(domain,pop_dat[pop_dat$A1==k,],Y_current[[k]],Nr_sum_current[k],Nr_obs_current[[k]],
                                       d_current,logit_r_hat[k],logit_r_hat_var[k])
        
        del_f <- f-mean(f) #normalize
        #restrict domain to densities >0
        zero_dens <- which(exp(del_f)==0)
        while(length(zero_dens)>0){
          domain <- domain[-zero_dens]
          f <- f[-zero_dens]
          del_f <- f-mean(f)
          zero_dens <- which(exp(del_f)==0)
        }
        
        target_density <- (exp(del_f)/sum(exp(del_f)))
        target_F <- cumsum(target_density)
        
        #take draw from inverse CDF
        U <- runif(1)
        return(min(domain[target_F>U]))
      })
      
      # draw new Ys from condtional distribution (ITS) ----
      for(k in 1:nrow(data)){
        area_id <- data$A1[k]
        within_area_id <- k - min(which(data$A1==area_id)) + 1
        # derive CDF
        domain <- data$Z[k]:round(0.5*data$N[k])
        f <- log_target_cond_5.3_y(domain,data[data$A1==area_id,],Yplus_current[area_id],pop_dat[pop_dat$A1==area_id,],
                                   Y_current[[area_id]],within_area_id,Nr_sum_current[area_id],Nr_obs_current[[area_id]],d_current)
        del_f <- f-mean(f) #normalize
        #restrict domain to densities >0
        zero_dens <- which(exp(del_f)==0)
        while(length(zero_dens)>0){
          domain <- domain[-zero_dens]
          f <- f[-zero_dens]
          del_f <- f-mean(f)
          zero_dens <- which(exp(del_f)==0)
        }
        
        target_density <- (exp(del_f)/sum(exp(del_f)))
        target_F <- cumsum(target_density)
        #take draw from inverse CDF
        U <- runif(1)
        Y_current[[area_id]][within_area_id] <- min(domain[target_F>U])
      }
      
      # record current state ----
      tau_postsamp_t[i] <- tau_current
      d_postsamp_t[i] <- d_current
      if(!skip.it){
        reg_params_postsamp_t[i,] <- reg_params_current
      }
      Yplus_postsamp_t[i,] <- Yplus_current
      data_table[, ':='(Y=unlist(Y_current))]
      Y_sum_obs_current <- (data_table[, sum(Y),by=A1])$V1
    } 
    
    results <- list(reg_params_t=reg_params_postsamp_t[(burnin+1):n_iter,], 
                    d_t=d_postsamp_t[(burnin+1):n_iter], 
                    tau_t=tau_postsamp_t[(burnin+1):n_iter], 
                    Yplus_t=Yplus_postsamp_t[(burnin+1):n_iter,], 
                    acc_reg_t = acc_reg, acc_d_t = acc_d)
    
    return(results)
  },mc.cores = chains)
  
  #combine chains
  b_postsamp <- reg_params_postsamp <- Yplus_postsamp <- NULL
  for(j in 1:chains){
    b_postsamp <- rbind(b_postsamp,parallel_results[[j]]$reg_params_t[,3:(2+n_admin2)])
    reg_params_postsamp <- rbind(reg_params_postsamp,parallel_results[[j]]$reg_params_t)
    Yplus_postsamp <- rbind(Yplus_postsamp,parallel_results[[j]]$Yplus_t)
  }
  tau_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$tau_t}))
  d_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$d_t}))
  alphaU_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,1]}))
  alphaR_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,2]}))
  names(tau_postsamp) <- names(d_postsamp) <- names(alphaU_postsamp) <- names(alphaR_postsamp) <- NULL
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp,alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  
  acc_reg_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_reg_t})))/((n_iter-burnin)*chains)
  acc_d_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_d_t})))/((n_iter-burnin)*chains)
  
  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              tau=tau_postsamp, 
              reg_params=reg_params_postsamp,
              d=d_postsamp, 
              Yplus=Yplus_postsamp,
              acc_reg_rate = acc_reg_rate, 
              acc_d_rate = acc_d_rate))
}

#######################################################################
#######################################################################
# Simulations for 5.3 ------

## set populations and simulation parameters
{
  ## this section is calibrated to be relatively similar to Sierra Leone 2019 survey, in terms of sample sizes
  n_admin1 <- 5 #numbers of admin1 areas
  n_admin2_per_admin1 <- c(3,4,3,4,2) #number of admin2 areas in each admin1 area
  n_admin2 <- sum(n_admin2_per_admin1)
  n_clusters_urban <- rnegbin(n_admin2,200,5) #number of urban clusters in each admin2 area
  n_clusters_rural <- rnegbin(n_admin2,300,5) #number of rural clusters in each admin2 area
  n_births_urban <- 75 #average number of births per urban cluster (0.3*number of households) 
  n_births_rural <-  100 #average number of births per rural cluster (0.45*number of households) 
  n_clusters_urban_samp <-  round(0.075*n_clusters_urban) #number of urban clusters sampled in each admin2 area
  n_clusters_rural_samp <-  round(0.05*n_clusters_rural) #number of rural clusters sampled in each admin2 area
  n_births_urban_samp <- 20 # average number of births sampled from each urban cluster 
  n_births_rural_samp <- 30 # average number of births sampled from each rural cluster 
  
  # intercept (fixed effect) 
  alphaU <- rnorm(1,-3.5,0.25)
  alphaR <- rnorm(1,-3.2,0.25)
  
  # draw hyperparameters
  tau <- rgamma(1,10,.25)
  
  # draw overdispersion parameter
  d <- rlnorm(1,0,0.25)
  
  # draw random effects
  b <- as.vector(Rfast::rmvnorm(1,rep(0,n_admin2),diag(1/tau,n_admin2)))
  
  true_params <- c(log(d),alphaU,alphaR,b)
  true_r <- exp(c(alphaU +b,alphaR +b))
}

## generate data 
{
  admin_key <- data.frame(A1 = unlist(sapply(1:n_admin1,function(x){rep(x,n_admin2_per_admin1[x])})), A2=1:n_admin2)
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
  pop_dat <- all_dat %>% group_by(A1,A2,U) %>% summarise(N=sum(N))
  pop_dat$num_clust_samp <- pop_dat$num_births_samp <- NA
  pop_dat$num_clust_samp[pop_dat$U==0] <- n_clusters_rural_samp
  pop_dat$num_clust_samp[pop_dat$U==1] <- n_clusters_urban_samp
  pop_dat$num_births_samp[pop_dat$U==0] <- n_births_rural_samp
  pop_dat$num_births_samp[pop_dat$U==1] <- n_births_urban_samp
  pop_dat$wt <- pop_dat$N/(pop_dat$num_births_samp*pop_dat$num_clust_samp)
  pop_dat$wt <- pop_dat$wt/100
  obs_dat <- merge(obs_dat,pop_dat[c('A2','U','wt')])
  obs_dat <- obs_dat[order(obs_dat$A2),]
  
  data_list <- list(obs_dat = obs_dat, #observed data
                    Yplus = all_dat %>% group_by(A1) %>% summarise(Yplus = sum(Y)), #ALL deaths in each admin1 area
                    pop_strata_dat = pop_dat) # number of births per strata (admin2 x urban)
  
  #what is the implied national urban and rural rates
  pop_dat <- pop_dat %>% mutate(alpha=ifelse(U==1,alphaU,alphaR))
  pop_dat$b <- b[pop_dat$A2]
  pop_dat$r <- exp(pop_dat$alpha+pop_dat$b)
  natl_r <- as.numeric(pop_dat %>% ungroup() %>% summarise(r = sum(r*N)/sum(N)))
  
}

## get admin1 direct estimates
{
  #obtain a direct estimate and variance for national rate
  data_for_direct <- data_list$obs_dat
  data_for_direct$years <- 'All'
  data_for_direct$died <- data_for_direct$Z
  data_for_direct$strata <- paste0(data_for_direct$A2,'.',c('rural','urban')[data_for_direct$U+1])
  data_for_direct$age <- '0'
  
  dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                               regionVar = "A1",timeVar = "years", clusterVar =  "~cluster",
                               Ntrials = "n",weightsVar = "wt")
  dir.est <- dir.est[dir.est$region!='All',]
}

## run model with known cluster size
{
  data_list5.3 <- data_list
  data_list5.3$obs_dat <- data_list$obs_dat %>% dplyr::select(-c(Y))
  
  #run model
  start.time <- Sys.time()
  postsamp_5.3 <- postsamp_Alg5.3_MALA_MCMC(eps=0.65, prop_sd_logd=0.25,
                                            data_list = data_list5.3, 
                                            logit_r_hat=dir.est$logit.est, 
                                            logit_r_hat_var = dir.est$var.est,
                                            n_iter=5000,chains = 4,burnin=0)
  Sys.time() - start.time
  
}

## estimate cluster sizes (need to choose an SD)
{
  frame.info <- all_dat %>% group_by(A2,U) %>% summarise(Mbar = mean(N))
  frame.info$n_clusters <- frame.info$n_clusters_samp <- NA
  frame.info[frame.info$U==1,]$n_clusters <- n_clusters_urban
  frame.info[frame.info$U==1,]$n_clusters_samp <- n_clusters_urban_samp
  frame.info[frame.info$U==0,]$n_clusters <- n_clusters_rural
  frame.info[frame.info$U==0,]$n_clusters_samp <- n_clusters_rural_samp
  
  frame.info$sd_M <- 0.05*frame.info$Mbar
  
  n_iter <- 1000
  M_obs_bar <- rep(NA,n_iter)
  frame.info$N <- NA
  
  for(i in 1:nrow(frame.info)){
    for(k in 1:n_iter){
      M = round(rnorm(frame.info$n_clusters[i],frame.info$Mbar[i],frame.info$sd_M[i]))
      M_obs <- M[sample(1:frame.info$n_clusters[i],frame.info$n_clusters_samp[i],prob = M)]
      M_obs_bar[k] <- mean(M_obs)
    }
    frame.info$N[i] <- round(mean(M_obs_bar))
  }
  data_list5.3$obs_dat <- merge(data_list5.3$obs_dat,frame.info[,c('A2','U','N')],by=c('A2','U')) %>% rename(N=N.y)
}

## run model with estimated cluster size
{
  #run model
  start.time <- Sys.time()
  postsamp_5.3_SD.05 <- postsamp_Alg5.3_MALA_MCMC(eps=0.75, k_Yplus = 500, prop_sd_logd=0.25,
                                                  data_list = data_list5.3, 
                                                  logit_r_hat=dir.est$logit.est[2:(n_admin1+1)], 
                                                  logit_r_hat_var = dir.est$var.est[2:(n_admin1+1)],
                                                  n_iter=2000,chains = 5)
  Sys.time() - start.time
  
}

postsamp_5.3$acc_reg_rate
postsamp_5.3$acc_d_rate
sum(is.na(postsamp_5.3$reg_params[,1]))

#just look at some set of iterations for each chain
iter_include <- 1000:5000
n_iter <- 5000 
chains <- 4
#iter_id <- 1:(chains*(n_iter))
iter_id <- as.vector(sapply(1:chains,function(x){iter_include+(x-1)*n_iter}))
iter_id <- iter_id[iter_id%%10==0]

par(mfrow=c(2,3))
plot(log(postsamp_5.3$d[iter_id]),type='l',ylab='',xlab='',main='log(d)')
#abline(h=(d),col='red')
plot(postsamp_5.3$reg_params[iter_id,1],type='l',ylab='',xlab='',main='alphaU')
#abline(h=alphaU,col='red')
plot(postsamp_5.3$reg_params[iter_id,2],type='l',ylab='',xlab='',main='alphaR')

plot(postsamp_5.3$reg_params[iter_id,3],type='l',ylab='',xlab='',main='b_3')
plot(postsamp_5.3$reg_params[iter_id,8],type='l',ylab='',xlab='',main='b_8')
plot(postsamp_5.3$reg_params[iter_id,15],type='l',ylab='',xlab='',main='b_15')

each_chain <- nrow(postsamp_5.3$Yplus)/chains
par(mfrow=c(2,2))
for(k in 1:ncol(postsamp_5.3$Yplus)){
  plot(cumsum(postsamp_5.3$Yplus[1:each_chain,k])/(1:each_chain),type='l',col=1,ylab='',
       ylim=c(min(postsamp_5.3$Yplus[,k]),max(postsamp_5.3$Yplus[,k])),main = paste0('Running mean of Yplus_',k))
  for(c in 2:chains){
    lines(cumsum(postsamp_5.3$Yplus[((c-1)*each_chain+1):(c*each_chain),k])/(1:each_chain),col=c)
  }
}

postsamp_5.3$true_d <- d
postsamp_5.3$true_alphaU <- alphaU
postsamp_5.3$true_alphaR <- alphaR
postsamp_5.3$true_b <- b
postsamp_5.3$true_Yplus <- data_list5.3$Yplus$Yplus
#save(postsamp_5.3,file="Handcoded MCMC Simulations/Alg 5.3 MALA 240205, 4 chains of 5k iterations posterior dist.rda")


#######################################################################
#######################################################################
# Alg 5.3a (Conditioning on Admin1, BYM2): P(alphaU,alphaR,b,phi,tau,d,Y^s,Y+|r+_est,Z,delta) ------

tau_init=40
phi_init <- 0.7
d_init = 1
alphaU_init = seq(-4,-3,by=0.25)
alphaR_init= seq(-3.5,-2.5,by=0.25)  #initial values
a_tau <- b_tau <- 0.01
mu_U = -3.5
mu_R = -3
sd_U<- sd_R <- 9
mu_logd = 0
sd_logd = 1 #hyperpriors
U_tau <- 0.5
alpha_tau <- 2/3

eps=0.65
prop_sd_tau=0.4
prop_sd_phi=0.4
prop_sd_logd=0.3
logit_r_hat = dir.est$logit.est[-1]
logit_r_hat_var = dir.est$var.est[-1]
n_iter=10
burnin=0

#preconditioned MALA as described in "Adaptive step size selection for Hessian-based manifold Langevian samplers" Tore Selland Kleppe
postsamp_Alg5.3a_MALA_MCMC <- function(tau_init=40, phi_init=0.7, d_init = 1, alphaU_init = seq(-4,-3,by=0.25), alphaR_init= seq(-3.5,-2.5,by=0.25), #initial values
                                       alpha_tau = 2/3,U_tau = 0.5,mu_U = -3.5,sd_U = 3,mu_R = -3.2,sd_R = 3,mu_logd = 0,sd_logd = 1,#hyperpriors  
                                       eps, prop_sd_tau, prop_sd_phi, prop_sd_logd,  #tuning params for regression param,d, and Yplus proposals
                                       logit_r_hat,logit_r_hat_var,data_list, 
                                       n_iter, burnin=n_iter/2,chains=4){
  
  data <- data_list$obs_dat
  data <- data[order(data$A2),]
  pop_dat <- data_list$pop_strata_dat
  #make data tables
  data_table <- as.data.table(data)
  pop_dt <- as.data.table(pop_dat)
  
  admin1_pop_dat <- pop_dat %>% group_by(A1) %>% summarise(N = sum(N))
  Yplus_init <- round(admin1_pop_dat$N*expit(logit_r_hat))
  W_inv_init <- solve(W(phi_init,Q_scaled_inv))
  
  n_admin1 <- length(unique(data_list$pop_strata_dat$A1))
  n_admin2 <- nrow(data_list$Q_scaled_inv)
  
  parallel_results <- mclapply(1:chains, function(j){
    
    b_init <- rnorm(n_admin2,0,1/sqrt(tau_init))
    
    ## ensure we are starting with a PD negative hessian 
    reg_params_init_grid <- as.matrix(expand.grid(alphaU_init,alphaR_init))
    search_order <- sample(1:nrow(reg_params_init_grid),nrow(reg_params_init_grid))
    
    for(j in 1:nrow(reg_params_init_grid)){
      reg_params_init <-  c(reg_params_init_grid[search_order[j],],b_init)
      
      Y_init <- sapply(1:n_admin1,function(i){
        
        subdat <- data %>% filter(A1==i)
        subpop_dat <- pop_dat %>% filter(A1==i)
        
        Nr_obs_init <- subdat$N*exp(subdat$U*reg_params_init[1] + (1-subdat$U)*reg_params_init[2] + reg_params_init[2+subdat$A2])
        Nr_sum_init <- sum(subpop_dat$N*exp(subpop_dat$U*reg_params_init[1] + (1-subpop_dat$U)*reg_params_init[2] + reg_params_init[2+subpop_dat$A2]))
        
        Y_init_t <- rmultinom(1,Yplus_init[i],prob=c(Nr_obs_init,Nr_sum_init-sum(Nr_obs_init)))
        Y_init_t <- Y_init_t[-length(Y_init_t)]
        Y_init_t[Y_init_t<subdat$Z] <- subdat$Z[Y_init_t<subdat$Z]
        
        return(Y_init_t)
      })
      data_table[, ':='(Y=unlist(Y_init))]
      Y_sum_obs_init <- (data_table[, sum(Y),by=A1])$V1
      
      hess_init <- pracma::hessian(log_target_5.3a_reg_params,x=reg_params_init,data_table=data_table,Y_sum_obs=Y_sum_obs_init,
                                   Yplus = Yplus_init,pop_dt=pop_dt,d=d_init,
                                   tau=tau_init,W_inv=W_inv_init,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      
      if(is.positive.definite(-hess_init)){
        break
      }
      print(j)
    }
    
    ## initialize data frames 
    tau_postsamp_t <- c(tau_init, rep(NA,n_iter-1))
    phi_postsamp_t <- c(phi_init, rep(NA,n_iter-1))
    d_postsamp_t <- c(d_init, rep(NA,n_iter-1))
    reg_params_postsamp_t <- rbind(reg_params_init, matrix(NA,n_iter-1,length(reg_params_init)))
    Yplus_postsamp_t <- rbind(Yplus_init, matrix(NA,n_iter-1,n_admin1))
    acc_reg <- acc_d <- acc_hyper <- 0
    
    # set current states
    tau_current <- tau_init
    phi_current <- phi_init
    W_inv_current <- solve(W(phi_current,Q_scaled_inv))
    d_current <- d_init
    alphaU_current <- reg_params_init[1]
    alphaR_current <- reg_params_init[2]
    b_current <- b_init
    reg_params_current <- reg_params_init
    Y_current <- Y_init
    Yplus_current <- Yplus_init
    Y_sum_obs_current <- Y_sum_obs_init
    
    for(i in 2:n_iter){
     # start.time <- Sys.time()
      # draw proposal for tau and phi ----
      log_tau_proposal <- rnorm(1,log(tau_current),prop_sd_tau)
      logit_phi_proposal <- rnorm(1,logitlink(phi_current),prop_sd_phi)
      W_inv_proposal <- solve(W(expit(logit_phi_proposal),Q_scaled_inv))
      
      a_prob_reg <- exp({
        #likelihood
        (n_admin2-1)/2*log_tau_proposal - 0.5*exp(log_tau_proposal)*t(b_current)%*%W_inv_proposal%*%b_current + log(alpha_tau)/U_tau*exp(-0.5*log_tau_proposal) -
          ((n_admin2-1)/2*log(tau_current) - 0.5*exp(log(tau_current))*t(b_current)%*%W_inv_current%*%b_current + log(alpha_tau)/U_tau*exp(-0.5*log(tau_current))) +
          #proposal dist (jacobian for logit phi)
          log(expit(logit_phi_proposal)*(1-expit(logit_phi_proposal))) - log(phi_current*(1-phi_current))
      })
      
      if(a_prob_reg>runif(1)){
        tau_current <- exp(log_tau_proposal)
        phi_current <- expit(logit_phi_proposal)
        W_inv_current <- solve(W(phi_current,Q_scaled_inv))
        if(i>burnin)
          acc_hyper <- acc_hyper + 1
      }
      
      # draw proposal for regression params ----
      grad_current <- pracma::grad(log_target_5.3a_reg_params,x=reg_params_current,data_table=data_table,d=d_current,Yplus = Yplus_current, pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,
                                   tau=tau_current,W_inv=W_inv_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      hess_current <- pracma::hessian(log_target_5.3a_reg_params,x=reg_params_current,data_table=data_table,d=d_current,Yplus = Yplus_current,pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,
                                      tau=tau_current,W_inv=W_inv_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      if(is.positive.definite(-hess_current)){
        cov_current <- matlib::inv(-hess_current)
      }else{
        message(paste0('current Hessian is not PD at iteration ',i))
        hess_chol <- get.GMW.chol(-hess_current,0.01)
        cov_current <- matlib::inv(hess_chol%*%t(hess_chol))
      }
      L_current <- t(chol(cov_current))
      
      r <- rnorm(length(reg_params_current),0,1)
      reg_params_proposal <- reg_params_current + eps^2/2*L_current%*%t(L_current)%*%grad_current + eps*L_current%*%r
      
      grad_proposal <- pracma::grad(log_target_5.3a_reg_params,x=reg_params_proposal,data_table=data_table,d=d_current,Yplus = Yplus_current,pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,
                                    tau=tau_current,W_inv=W_inv_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      hess_proposal <- pracma::hessian(log_target_5.3a_reg_params,x=reg_params_proposal,data_table=data_table,d=d_current,Yplus = Yplus_current,pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,
                                       tau=tau_current,W_inv=W_inv_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R)
      
      if(isSymmetric(hess_proposal)){
        if(is.positive.definite(-hess_proposal)){
          cov_proposal <- matlib::inv(-hess_proposal)
          if(is.positive.definite(cov_proposal)){
            skip.it <- F
          }else{
            message(paste0('Proposed Hessian not PD at iteration ',i,'\n'))
            skip.it <- T
          }
        }else{
          message(paste0('Proposed Hessian not PD at iteration ',i,'\n'))
          skip.it <- T}
      }else{
        message(paste0('Proposed Hessian not symmetric at iteration ',i,'\n'))
        skip.it <- T}
      
      # skip.it <- T
      # 
      # log_target_current <- log_target_5.3a_reg_params(reg_params_current,data_table,Yplus_current,pop_dt,d_current,Y_sum_obs_current,tau_current,W_inv_current,mu_U,sd_U,mu_R,sd_R)
      # log_target_proposal <- log_target_5.3a_reg_params(reg_params_proposal,data_table,Yplus_current,pop_dt,d_current,Y_sum_obs_current,tau_current,W_inv_current,mu_U,sd_U,mu_R,sd_R)
      # 
      # L_vals <- optim(par=L_current[lower.tri(L_current,diag=T)],taylor.expand,method ="BFGS",
      #                 delta=reg_params_proposal-reg_params_current,
      #                 delta_log_target=log_target_proposal - log_target_current,grad=grad_proposal)$par
      # L_proposal_inv <- matrix(0,length(reg_params_current),length(reg_params_current))
      # L_proposal_inv[lower.tri(L_proposal_inv,diag=T)] <- L_vals
      # 
      # L_proposal <- t(chol(solve(L_proposal_inv%*%t(L_proposal_inv))))
      
      if(!skip.it){
        # accept or reject proposal for (alpha, b, d)
        L_proposal <- t(chol(cov_proposal))
        
        a_prob_reg <- min(1,exp({
          log_target_5.3a_reg_params(reg_params_proposal,data_table=data_table,Yplus = Yplus_current,pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,d=d_current,tau=tau_current,W_inv=W_inv_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R) +
            (-0.5*sum((solve(L_proposal)%*%(reg_params_current-reg_params_proposal)/eps - 0.5*eps*t(L_proposal)%*%grad_proposal)^2)) - 
            log_target_5.3a_reg_params(reg_params_current,data_table=data_table,Yplus = Yplus_current,pop_dt=pop_dt,Y_sum_obs=Y_sum_obs_current,d=d_current,tau=tau_current,W_inv=W_inv_current,mu_U=mu_U,sd_U=sd_U,mu_R=mu_R,sd_R=sd_R) -
            (-0.5*sum(r^2))
        }))
        
        if(a_prob_reg>runif(1)){
          reg_params_current <- reg_params_proposal
          alphaU_current <- reg_params_proposal[1]
          alphaR_current <- reg_params_proposal[2]
          b_current <- reg_params_proposal[3:length(reg_params_proposal)]
          if(i>burnin){acc_reg <- acc_reg + 1}
        }
      }
      
      Nr_sum_current <- sapply(1:n_admin1,function(i){
        subpop_dat <- pop_dat %>% filter(A1==i)
        Nr_sum_current_t <- sum(subpop_dat$N*exp(subpop_dat$U*alphaU_current + (1-subpop_dat$U)*alphaR_current + b_current[subpop_dat$A2]))
        return(Nr_sum_current_t)
      })
      
      Nr_obs_current <- sapply(1:n_admin1,function(i){
        subdat <- data %>% filter(A1==i)
        Nr_obs_t <- subdat$N*exp(subdat$U*alphaU_current + (1-subdat$U)*alphaR_current + b_current[subdat$A2])
        return(Nr_obs_t)
      })
      
      # draw proposal for d -------
      logd_proposal <- rnorm(1,log(d_current),prop_sd_logd)
      
      a_prob_d <- exp({log_target_5.3_logd(logd_proposal,alphaU_current,alphaR_current,b_current,data,
                                           Yplus_current,pop_dat,Y_current,mu_logd,sd_logd) -
          log_target_5.3_logd(log(d_current),alphaU_current,alphaR_current,b_current,data,
                              Yplus_current,pop_dat,Y_current,mu_logd,sd_logd)})
      
      if(a_prob_d > runif(1)){
        d_current <- exp(logd_proposal)
        if(i>burnin){
          acc_d <- acc_d + 1
        }
      }
      
      # draw for Yplus's from conditional distribution (ITS) ----
      
      Yplus_current <- sapply(1:n_admin1,function(k){
        # calculate Poisson parameter for domain centered on current Y
        domain <- round(0.5*Yplus_current[k]):(2*Yplus_current[k])
        f <- log_target_cond_5.3_Yplus(domain,pop_dat[pop_dat$A1==k,],Y_current[[k]],Nr_sum_current[k],Nr_obs_current[[k]],
                                       d_current,logit_r_hat[k],logit_r_hat_var[k])
        domain <- domain[f>-Inf]
        f <- f[f>-Inf]
        del_f <- f-mean(f) #normalize
        #restrict domain to densities >0
        zero_dens <- which(exp(del_f)==0)
        while(length(zero_dens)>0){
          domain <- domain[-zero_dens]
          f <- f[-zero_dens]
          del_f <- f-mean(f)
          zero_dens <- which(exp(del_f)==0)
        }
        
        target_density <- (exp(del_f)/sum(exp(del_f)))
        target_F <- cumsum(target_density)
        
        #take draw from inverse CDF
        U <- runif(1)
        if(is.na(min(domain[target_F>U]))){
          stop('New Yplus is NA')
        }
        return(min(domain[target_F>U]))
      })
      
      # draw new Ys from condtional distribution (ITS) ----
      for(k in 1:nrow(data)){
        area_id <- data$A1[k]
        within_area_id <- k - min(which(data$A1==area_id)) + 1
        
        # derive CDF
        domain <- data$Z[k]:round(0.5*data$N[k])
        f <- log_target_cond_5.3_y(domain,data[data$A1==area_id,],Yplus_current[area_id],
                                   Y_current[[area_id]],within_area_id,Nr_sum_current[area_id],Nr_obs_current[[area_id]],d_current)
        domain <- domain[f>-Inf]
        f <- f[f>-Inf]
        del_f <- f-mean(f) #normalize
        #restrict domain to densities >0
        zero_dens <- which(exp(del_f)==0)
        while(length(zero_dens)>0){
          domain <- domain[-zero_dens]
          f <- f[-zero_dens]
          del_f <- f-mean(f)
          zero_dens <- which(exp(del_f)==0)
        }
        
        target_density <- (exp(del_f)/sum(exp(del_f)))
        target_F <- cumsum(target_density)
        #take draw from inverse CDF
        U <- runif(1)
        Y_current[[area_id]][within_area_id] <- min(domain[target_F>U])
      }
      
      Sys.time() - start.time
      # record current state ----
      tau_postsamp_t[i] <- tau_current
      phi_postsamp_t[i] <- phi_current
      d_postsamp_t[i] <- d_current
      if(!skip.it){
        reg_params_postsamp_t[i,] <- reg_params_current
      }
      Yplus_postsamp_t[i,] <- Yplus_current
      data_table[, ':='(Y=unlist(Y_current))]
      Y_sum_obs_current <- (data_table[, sum(Y),by=A1])$V1
      
    } 
    
    results <- list(reg_params_t=reg_params_postsamp_t[(burnin+1):n_iter,], 
                    d_t=d_postsamp_t[(burnin+1):n_iter], 
                    tau_t=tau_postsamp_t[(burnin+1):n_iter], 
                    phi_t=phi_postsamp_t[(burnin+1):n_iter], 
                    Yplus_t=Yplus_postsamp_t[(burnin+1):n_iter,], 
                    acc_reg_t = acc_reg, acc_d_t = acc_d,acc_hyper_t = acc_hyper)
    
    return(results)
  },mc.cores = chains)
  
  #combine chains
  b_postsamp <- reg_params_postsamp <- Yplus_postsamp <- NULL
  for(j in 1:chains){
    b_postsamp <- rbind(b_postsamp,parallel_results[[j]]$reg_params_t[,3:(2+n_admin2)])
    reg_params_postsamp <- rbind(reg_params_postsamp,parallel_results[[j]]$reg_params_t)
    Yplus_postsamp <- rbind(Yplus_postsamp,parallel_results[[j]]$Yplus_t)
  }
  tau_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$tau_t}))
  phi_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$phi_t}))
  d_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$d_t}))
  alphaU_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,1]}))
  alphaR_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,2]}))
  names(tau_postsamp) <- names(d_postsamp) <- names(alphaU_postsamp) <- names(alphaR_postsamp) <- NULL
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp,alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  
  acc_reg_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_reg_t})))/((n_iter-burnin)*chains)
  acc_d_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_d_t})))/((n_iter-burnin)*chains)
  acc_hyper_rate <- sum(unlist(lapply(1:chains,function(x){parallel_results[[x]]$acc_hyper_t})))/((n_iter-burnin)*chains)
  
  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              tau=tau_postsamp, 
              phi=phi_postsamp,
              reg_params=reg_params_postsamp,
              d=d_postsamp, 
              Yplus=Yplus_postsamp,
              acc_reg_rate = acc_reg_rate, 
              acc_d_rate = acc_d_rate,
              acc_hyper_rate = acc_hyper_rate))
}

#######################################################################
#######################################################################
# Applying Alg 5.3/5.3a to real DHS data ------
country_t <- 'Sierra_Leone'
year_range <- 2017:2019
survey_year <- 2019

country_t <- 'Malawi'
year_range <- 2013:2015
survey_year <- 2015

country_t <- 'Mauritania'
year_range <- 2016:2020
survey_year <- 2018

# old countries
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/",country_t,"_cluster_dat_1frame.rda"))
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/shapeFiles_gadm/",country_t,"_Amat.rda"))
#load("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Sierra_Leone/shapeFiles/Sierra_Leone_Amat.rda")

#Mauritania
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/",country_t,"_cluster_dat.rda"))
mod.dat <- data_for_direct <- nmr.dat
mod.dat$years <- max(year_range)
mod.dat$survey <- survey_year
mod.dat$age <- 0

load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/shapeFiles/gadm41_MRT_shp/",country_t,"_Amat.rda"))

load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/worldpop/adm2_weights_u1.rda"))
admin_key <- readxl::read_excel("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Admin1_Admin2_Key.xlsx") %>% filter(country==country_t) %>% dplyr::select(-country)
admin_key$A1 <- as.numeric(sapply(1:nrow(admin_key),function(i){str_split(admin_key$admin1.char[i],'_')[[1]][2]}))
admin_key <- admin_key[order(admin_key$A1),]
admin_key$A2 <-  as.numeric(sapply(1:nrow(admin_key),function(i){str_split(admin_key$admin2.char[i],'_')[[1]][2]}))

# get generalized inverse of scaled Q for BYM2
{
  Q.admin2 <- -admin2.mat
  Q.admin2 <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})*Q.admin2
  diag(Q.admin2) <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})
  diag(Q.admin2)[diag(Q.admin2)==0] <- 1
  Q_scaled <- INLA::inla.scale.model(Q.admin2, constr=list(A=t(eigen(Q.admin2)$vectors[,eigen(Q.admin2)$values<1e-10]), e=rep(0,sum(eigen(Q.admin2)$values<1e-10))))
  Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))
}

# prepare direct estimates
{
  data_for_direct <- mod.dat %>% filter(years %in% year_range, survey==survey_year,age==0)
  data_for_direct$years <- 'All'
  data_for_direct$died <- data_for_direct$Y
  dir.est <- SUMMER::getDirect(data_for_direct, years='All',
                               regionVar = "admin1.char",timeVar = "years", clusterVar =  "~cluster",
                               Ntrials = "total",weightsVar = "v005")
  natl.dir <- dir.est[dir.est$region=='All',]
  dir.est <- dir.est[dir.est$region!='All',]
  dir.est <- dir.est[order(as.numeric(sapply(1:nrow(dir.est),function(i){str_split(dir.est$region[i],'_')[[1]][2]}))),]
}

# input strata level DHS survey info 
{
  dhs_dat <- rbind(admin_key,admin_key)
  dhs_dat$U <- c(rep(1,nrow(admin_key)),rep(0,nrow(admin_key)))
  dhs_dat$A1 <- as.numeric(sapply(1:nrow(dhs_dat),function(i){str_split(dhs_dat$admin1.char[i],'_')[[1]][2]}))
  dhs_dat <- dhs_dat[order(dhs_dat$A1,dhs_dat$A2),]
  if(country_t=='Sierra_Leone'){
    # total number of HH
    dhs_dat$M <- c(24189,61930,48343,63391,26032,57316,16256,37570,3055,43038,29319,67872,
                   8530,20844,30545,46499,1634,25100,18393,65854,6214,26324,34311,68412,4489,57391,4416,47098,82386,8898,229951,0)
    # average number of HH per cluster
    dhs_dat$M_bar <- c(120,106,110,93,94,93,81,100,118,105,100,96,69,77,115,100,68,76,89,79,88,67,106,97,121,99,134,86,130,137,108,0)
    # total number of clusters in strata
    dhs_dat$n_clusters <- c(201,586,441,678,276,615,200,376,26,409,294,706,123,271,266,464,24,330,207,834,71,390,323,708,37,579,33,549,635,65,2139,0)
    dhs_dat$n_HH_samp <- 24 #number of HH sampled in each cluster
  }else if(country_t=='Malawi'){
    # total number of HH
    dhs_dat$M <- c(2924,34856,8574,49234,299,1721,31061,138777,2276,39993,3847,32190,4489,141389,4479,117405,8964,118301,153717,275194,3570,93639,
                   5010,57458,3306,110485,1555,45873,6089,71442,5037,70619,153578,80879,2830,95205,592,70968,5303,109833,8473,177442,3243,124174,
                   3445,18573,366,25049,4227,48373,1117,75562,2405,139634,19041,142394)
    
    # total number of clusters in strata
    dhs_dat$n_clusters <- c(11,205,37,370,4,18,122,825,12,229,12,156,15,486,18,450,29,486,458,1173,12,374,16,177,11,468,6,204,22,416,17,275,412,381,
                            16,380,2,334,19,436,25,614,17,658,9,80,3,157,14,241,3,316,12,674,79,584)
    # average number of HH per cluster
    dhs_dat$M_bar <- round(dhs_dat$M/dhs_dat$n_clusters)
    dhs_dat <- dhs_dat %>% mutate(n_HH_samp=ifelse(U==1,30,33))
  }else if(country_t=='Mauritania'){
    # strata info only at admin1 level
    adm1_dhs_dat <- dhs_dat %>% dplyr::select(admin1.name,A1,U) %>% unique()
    # total number of HH in frame
    adm1_dhs_dat$M <- c(1261,1588,3512,10988,3108,9569,6103,292,3480,8243,1644,4659,3771,14383,2255,10531,1026,178,
                        61022,0,920,2561,2539,87,3456,7505)
    # total number of clusters in strata (from MICS)
    adm1_dhs_dat$n_clusters <- c(82,83,186,657,202,563,387,18,193,445,91,268,202,751,124,581,67,10,3634,0,
                                 51,148,164,6,220,412)
    # average number of HH per cluster
    adm1_dhs_dat$M_bar <- round(adm1_dhs_dat$M/adm1_dhs_dat$n_clusters)
    adm1_dhs_dat$n_HH_samp <- 10
  }
  
}

# organize survey data
{
  alg_dat <- mod.dat %>% filter(years %in% year_range, survey==survey_year,age==0) %>% mutate(U=ifelse(urban=='urban',1,0)) %>% 
    dplyr::select(cluster,total,Y,U,admin1,admin1.char,admin1.name,admin2.char,admin2.name) %>% 
    group_by(cluster,U,admin1,admin1.char,admin1.name,admin2.char,admin2.name) %>% 
    reframe(Z=sum(Y),n=sum(total)) %>% rename(A1=admin1)
  alg_dat <- merge(alg_dat,admin_key)
}

# estimate average number of HH in SAMPLED cluster
{
  if(country_t!='Mauritania'){ # if country has info at admin 2 level
    frame.info <- merge(dhs_dat,alg_dat %>% group_by(admin2.char,U) %>% summarise(n_clusters_samp = n())) %>% 
      dplyr::select(A2,admin2.char,U,M_bar,n_clusters,n_clusters_samp,n_HH_samp)
    sd_M <- 0.1*frame.info$M_bar
    
    n_iter <- 1000
    M_obs_bar <- rep(NA,n_iter)
    frame.info$M_bar_est <- NA
    
    for(i in 1:nrow(frame.info)){
      for(k in 1:n_iter){
        M = round(rnorm(frame.info$n_clusters[i],frame.info$M_bar[i],sd_M[i]))
        M_obs <- M[sample(1:frame.info$n_clusters[i],frame.info$n_clusters_samp[i],prob = M)]
        M_obs_bar[k] <- mean(M_obs)
      }
      frame.info$M_bar_est[i] <- mean(M_obs_bar)
    }
  }else{ #if country has info at admin1 level
    frame.info <- merge(adm1_dhs_dat,alg_dat %>% group_by(A1,U) %>% summarise(n_clusters_samp = n())) %>% 
      dplyr::select(A1,admin1.name,U,M_bar,n_clusters,n_clusters_samp,n_HH_samp)
    sd_M <- 0.1*frame.info$M_bar
    
    n_iter <- 1000
    M_obs_bar <- rep(NA,n_iter)
    frame.info$M_bar_est <- NA
    
    for(i in 1:nrow(frame.info)){
      if(frame.info$admin1.name[i] %in% c('Adrar','Tagant','Inchiri','Guidimaka')){
        frame.info$M_bar_est[i] <- frame.info$M_bar[i]
      }else{
        for(k in 1:n_iter){
          M = round(rnorm(frame.info$n_clusters[i],frame.info$M_bar[i],sd_M[i]))
          M_obs <- M[sample(1:frame.info$n_clusters[i],frame.info$n_clusters_samp[i],prob = M)]
          M_obs_bar[k] <- mean(M_obs)
        }
        frame.info$M_bar_est[i] <- mean(M_obs_bar)
      }
    }
  }
}

# add estimated number of total births per cluster
alg_dat <- merge(alg_dat,frame.info) %>% dplyr::select(-M_bar,-n_clusters,-n_clusters_samp) %>%
  mutate(N=round(M_bar_est*n/n_HH_samp))

# organize population and survey info
{
  if(country_t!='Mauritania'){
    pop_strata_dat <- merge(dhs_dat %>% dplyr::select(A1,A2,admin1.name,admin2.name,U,M),
                            alg_dat %>% group_by(A2,U) %>% summarise(birth_rate = sum(n)/sum(n_HH_samp)))
    pop_strata_dat$N <- round(pop_strata_dat$M*pop_strata_dat$birth_rate)
  }else{
    # get U1 population at admin2 level
    pop_est <- NULL
    for(year in year_range){
      adm2_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/UR Fractions/U1_fraction_MRT/admin2_u1_",year,"_urban_frac.rds"))
      adm2_UR_weights <- adm2_UR_weights[order(as.numeric(str_remove(adm2_UR_weights$adm_idx,'admin2_'))),]
      pop_est_tmp <- data.frame(admin2.char = rep(admin_key$admin2.char,each=2),U=rep(c(0,1),max(admin_key$A2)))
      pop_est_tmp <- merge(pop_est_tmp,weight.adm2.u1[weight.adm2.u1$year==year,c(1,3)])
      pop_est_tmp <- merge(pop_est_tmp,adm2_UR_weights[,2:3],by.y='adm_idx',by.x='admin2.char')
      pop_est_tmp <- pop_est_tmp %>% mutate(N=ifelse(U==1,admin_pop*urb_frac,admin_pop*(1-urb_frac)))
      pop_est_tmp$year <- year
      
      pop_est <- rbind(pop_est,pop_est_tmp[,c(1,2,5,6)])
    }
    pop_strata_dat <- pop_est %>% group_by(admin2.char,U) %>% summarise(N = round(sum(N)))
    pop_strata_dat <- merge(admin_key,pop_strata_dat)
    pop_strata_dat <- pop_strata_dat[order(pop_strata_dat$A2,pop_strata_dat$U),]
  }
  
  
  data_list <- list(obs_dat = alg_dat,
                    pop_strata_dat = as.data.frame(pop_strata_dat),
                    Q_scaled_inv = Q_scaled_inv)
  }

start.time <- Sys.time()
postsamp_5.3 <- postsamp_Alg5.3_MALA_MCMC(eps=0.6, prop_sd_logd=0.2,
                                          data_list = data_list, 
                                          logit_r_hat=dir.est$logit.est, 
                                          logit_r_hat_var = dir.est$var.est,
                                          n_iter=2000,chains = 5,burnin=00)
Sys.time() - start.time

postsamp_5.3$acc_reg_rate
postsamp_5.3$acc_d_rate
postsamp_5.3$acc_Yplus_rate

plot_id <- c(1000:2000,6000:7000,11000:12000,16000:17000,21000:22000)

plot(postsamp_5.3$alphaU[plot_id],type='l')
plot(postsamp_5.3$alphaR[plot_id],type='l')
plot(postsamp_5.3$r[plot_id,1],type='l')
plot(postsamp_5.3$r[plot_id,5],type='l')
plot(postsamp_5.3$r[plot_id,12],type='l')
plot(postsamp_5.3$Yplus[plot_id,1],type='l')
plot(postsamp_5.3$Yplus[plot_id,2],type='l')
plot(postsamp_5.3$Yplus[plot_id,3],type='l')
plot(postsamp_5.3$Yplus[plot_id,4],type='l')

hist(postsamp_5.3$alphaU[plot_id])
hist(postsamp_5.3$alphaR[plot_id])

## aggregate results
strat.res.dat <- dhs_dat[order(dhs_dat$U,decreasing=T),]
strat.res.dat$mean <- colMeans(postsamp_5.3$r[plot_id,],na.rm = T)
strat.res.dat$median <- colMedians(postsamp_5.3$r[plot_id,],na.rm = T)

adm2.res.dat <- merge(strat.res.dat,adm2_UR_weights[,2:3],by.x='admin2.char',by.y='adm_idx') %>% mutate(urb_frac = ifelse(U==0,1-urb_frac,urb_frac)) %>% 
  group_by(admin2.char,admin1.char) %>% summarise(mean = sum(urb_frac*mean),median = sum(urb_frac*median))

adm1.res.dat <- merge(adm2.res.dat,adm2_weights[,1:2],by.x='admin2.char',by.y='region') %>% group_by(admin1.char) %>%
  reframe(mean=sum(proportion*mean)/sum(proportion),median=sum(proportion*median)/sum(proportion))

ggplot() +# geom_point(data=adm2.res.dat,aes(median,admin1.char,col='2',pch='2'),size=4) +
  geom_point(data=dir.est,aes(mean,region,col='1',pch='1'),size=4) + 
  geom_point(data=adm1.res.dat,aes(median,admin1.char,col='3',pch='3'),size=4) +
  geom_point(data=all_res,aes(adm1_median_agg,admin1.char))+
  ggtitle('Sierra Leone 2017-2019') + theme_minimal() +
  ylab('') + xlab('NMR') + theme(legend.position = 'bottom') +
  scale_colour_manual(name = '', values =c('1'='blue','3'='red','2'='orange'),
                      labels = c('Admin 1 Direct','Admin 2 IID','Admin 2 IID, Aggregated')) +
  scale_shape_manual(name = '', values =c('1'=17,'2'=1,'3'=2),
                     labels = c('Admin 1 Direct','Admin 2 IID','Admin 2 IID, Aggregated'))

start.time <- Sys.time()
postsamp_5.3a <- postsamp_Alg5.3a_MALA_MCMC(eps=0.65, 
                                            prop_sd_logd=0.25, prop_sd_tau = 0.3, prop_sd_phi = 0.3,
                                            data_list = data_list, 
                                            logit_r_hat=dir.est$logit.est, 
                                            logit_r_hat_var = dir.est$var.est,
                                            n_iter=100,chains = 4,burnin=00)
Sys.time() - start.time

sum(is.na(postsamp_5.3a$reg_params[,1]))
postsamp_5.3a$acc_hyper_rate
postsamp_5.3a$acc_reg_rate
popostsamp_5.3a$acc_d_rate

#just look at some set of iterations for each chain
iter_include <- 1000:2000
n_iter <- 2000 
chains <- 4
#iter_id <- 1:(chains*(n_iter))
iter_id <- as.vector(sapply(1:chains,function(x){iter_include+(x-1)*n_iter}))
iter_id <- iter_id[iter_id%%10==0]

par(mfrow=c(1,3))
plot(log(postsamp_5.3a$d[iter_id]),type='l',ylab='',xlab='',main='log(d)')
plot(log(postsamp_5.3a$tau[iter_id]),type='l',ylab='',xlab='',main='log(tau)')
plot(logitlink(postsamp_5.3a$phi[iter_id]),type='l',ylab='',xlab='',main='logit(phi)')
#abline(h=(d),col='red')
par(mfrow=c(2,3))
plot(postsamp_5.3a$reg_params[iter_id,1],type='l',ylab='',xlab='',main='alphaU')
#abline(h=alphaU,col='red')
plot(postsamp_5.3a$reg_params[iter_id,2],type='l',ylab='',xlab='',main='alphaR')

plot(postsamp_5.3a$reg_params[iter_id,3],type='l',ylab='',xlab='',main='b_3')
plot(postsamp_5.3a$reg_params[iter_id,8],type='l',ylab='',xlab='',main='b_8')
plot(postsamp_5.3a$reg_params[iter_id,15],type='l',ylab='',xlab='',main='b_15')
plot(postsamp_5.3a$reg_params[iter_id,21],type='l',ylab='',xlab='',main='b_21')

each_chain <- nrow(postsamp_5.3a$Yplus)/chains
par(mfrow=c(1,3))
for(k in 1:ncol(postsamp_5.3a$Yplus)){
  plot(cumsum(postsamp_5.3a$Yplus[1:each_chain,k])/(1:each_chain),type='l',col=1,ylab='',
       ylim=c(min(postsamp_5.3a$Yplus[,k]),max(postsamp_5.3a$Yplus[,k])),main = paste0('Running mean of Yplus_',k))
  for(c in 2:chains){
    lines(cumsum(postsamp_5.3a$Yplus[((c-1)*each_chain+1):(c*each_chain),k])/(1:each_chain),col=c)
  }
}

save(postsamp_5.3a,file="Handcoded MCMC Simulations/Alg 5.3 BYM2 MALA 240206, MWI 2013-2015, 4 chains of 5k iterations posterior dist.rda")

adm2_UR_weights <- readRDS(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Results/",
                                  country_t,"/UR/U1_fraction/admin2_u1_2014_urban_frac.rds"))
adm2_UR_weights <- adm2_UR_weights[order(as.numeric(str_remove(adm2_UR_weights$adm_idx,'admin2_'))),]
load(paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/",country_t,"/worldpop/adm2_weights_u1.rda"))
adm2_weights <- weight.adm2.u1[weight.adm2.u1$years==2014,]

strat.res.dat <- dhs_dat[order(dhs_dat$U,decreasing=T),]
strat.res.dat$mean <- colMeans(postsamp_5.3a$r[iter_id,],na.rm = T)
strat.res.dat$median <- colMedians(postsamp_5.3a$r[iter_id,],na.rm = T)

adm2.res.dat <- merge(strat.res.dat,adm2_UR_weights[,2:3],by.x='admin2.char',by.y='adm_idx') %>% mutate(urb_frac = ifelse(U==0,1-urb_frac,urb_frac)) %>% 
  group_by(admin2.char,admin1.char) %>% summarise(mean = sum(urb_frac*mean),median = sum(urb_frac*median))

adm1.res.dat <- merge(adm2.res.dat,adm2_weights[,1:2],by.x='admin2.char',by.y='region') %>% group_by(admin1.char) %>%
  reframe(mean=sum(proportion*mean)/sum(proportion),median=sum(proportion*median)/sum(proportion))

ggplot() + geom_point(data=adm2.res.dat,aes(median,admin1.char,col='2',pch='2'),size=4) +
  geom_point(data=dir.est,aes(mean,region,col='1',pch='1'),size=4) + 
  geom_point(data=adm1.res.dat,aes(median,admin1.char,col='3',pch='3'),size=4) +
  geom_point(data=all_res,aes(adm1_median_agg,admin1.char,col='4',pch='4'),size=4) +
  ggtitle('Malawi 2013-2015') + theme_minimal() +
  ylab('') + xlab('NMR') + theme(legend.position = 'bottom') +
  scale_colour_manual(name = '', values =c('1'='blue','3'='red','2'='orange','4'='red'),
                      labels = c('Admin 1 Direct','Admin 2','Admin 2, Aggregated','Admin 2 INLA, Aggregated')) +
  scale_shape_manual(name = '', values =c('1'=17,'2'=1,'3'=2,'4'=17),
                     labels = c('Admin 1 Direct','Admin 2','Admin 2, Aggregated','Admin 2 INLA, Aggregated'))


#######################################################################
#######################################################################
# Alg 5.2a with NUTS ------

# after running data for Malawi
data_list$pop_strata_dat <- data_list$pop_strata_dat[order(1-data_list$pop_strata_dat$U,data_list$pop_strata_dat$A2),]
data_list$obs_dat$A <- data_list$obs_dat$A2
data_list$pop_strata_dat$A <- data_list$pop_strata_dat$A2

## find lambda to properly calibrate PC prior on phi
{
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

data_list$gamma_til <- gamma_til
data_list$lambda_opt <- lambda_opt
}

tau_init=40
phi_init=0.5
d_init = 1
alphaU_init = -3.5
alphaR_init= -3.2 #initial values
U_tau = 1
alpha_tau = 0.01
alpha_phi = 0.5
beta_phi =1
mu_U = -3.5
sd_U = 3
mu_R = -3.2
sd_R = 3
mu_d = 0
sd_d = 1

delta=0.5
max_treedepth = 10
eps0=0.01
kappa = 0.9
gamma = 0.5
t0=5
#M_diag = 1/c(1,2,1,rep(0.5,30),200),
M_diag = 1/c(0.25,0.5,0.5,rep(0.5,30),500)
logit_r_hat=natl.dir$logit.est
logit_r_hat_var=natl.dir$var.est
n_iter=100
warmup=50
chains=4

postsamp_Alg5.2a_NUTS <- function(tau_init=40, phi_init=0.5, d_init = 1, alphaU_init = -3.5, alphaR_init= -3.2, #initial values
                                  U_tau = 1,alpha_tau = 0.01, alpha_phi = 0.5, beta_phi =1,
                                  mu_U = -3.5,sd_U = 3,mu_R = -3.2,sd_R = 3,mu_d = 0,sd_d = 1,#hyperpriors  
                                  max_treedepth = 10, eps0=0.1, M_diag,
                                 # kappa = 0.9, t0=5, gamma=0.35, delta=0.5,
                                  logit_r_hat,logit_r_hat_var,data_list, 
                                  n_iter, warmup=500,chains=4){
  
  data <- data_list$obs_dat
  data <- data[order(data$A),]
  pop_dat <- data_list$pop_strata_dat
  Q_scaled_inv <- data_list$Q_scaled_inv
  n_areas <- length(unique(data_list$pop_strata_dat$A))
  
  #start.time <- Sys.time()
  parallel_results <- mclapply(1:chains, function(j){

    b_init <- rnorm(n_areas,0,1/sqrt(tau_init))
    Yplus_init <- round(sum(pop_dat$N)*expit(logit_r_hat))
    
    reg_params_init <-  c(log(tau_init),logitlink(phi_init),log(d_init),alphaU_init,alphaR_init,b_init)
    Nr_sum_init <- sum(pop_dat$N*exp(pop_dat$U*reg_params_init[4] + (1-pop_dat$U)*reg_params_init[5] + reg_params_init[5+pop_dat$A]))
    Nr_obs_init <- data$N*exp(data$U*reg_params_init[4] + (1-data$U)*reg_params_init[5] + reg_params_init[5+data$A])
    Y_init <- rmultinom(1,Yplus_init,prob=c(Nr_obs_init,Nr_sum_init-sum(Nr_obs_init)))
    Y_init <- Y_init[-length(Y_init)]
    Y_init[Y_init<data$Z] <- data_list$obs_dat$Z[Y_init<data$Z]
    
    ## initialize data frames 
    reg_params_postsamp_t <- rbind(reg_params_init, matrix(NA,n_iter,length(reg_params_init)))
    Y_postsamp_t <- rbind(Y_init, matrix(NA,n_iter,nrow(data)))
    Yplus_postsamp_t <- c(Yplus_init, rep(NA,n_iter))
    tree_depth_t <- eps_t <- rep(NA,n_iter)
    
    # set current states
    reg_params_current <- reg_params_init
    Y_current <- Y_init
    Yplus_current <- Yplus_init
    
    par_list <- list(eps=eps0,#M_adapt = warmup, 
                     M_diag=M_diag)
    
    for(i in 1:n_iter){
     # draw regression parameters using NUTS -----
    
      cons_params <- list(data=data,pop_dat=pop_dat,Q_scaled_inv=Q_scaled_inv,Y=Y_current,
                          logit_r_hat = logit_r_hat, logit_r_hat_var=logit_r_hat_var,
                          mu_d=mu_d,sd_d=sd_d,mu_U=mu_U,mu_R=mu_R,sd_U=sd_U,sd_R=sd_R,alpha_tau=alpha_tau,U_tau=U_tau,alpha_phi=alpha_phi,beta_phi=beta_phi)
     
      # if(i>20 & i<=warmup){
      #   par_list$M_diag <- -pracma::hessdiag(log_target_5.2a_reg_params_4NUTS,x=c(colMeans(reg_params_postsamp_t,na.rm = T),
      #                                                                             mean(Yplus_postsamp_t,na.rm=T)),
      #                                        cons_params=cons_params)
      # }
      
      nuts <- NUTS_one_step(c(reg_params_current,Yplus_current), i, f = log_target_5.2a_reg_params_4NUTS, grad_f=grad_log_target_5.2a_4NUTS, cons_params, 
                           par_list= par_list, max_treedepth = 10, verbose = F)
      
      reg_params_current <- nuts$theta[-length(nuts$theta)]
      Yplus_current <- round(nuts$theta[length(nuts$theta)])
      par_list <- nuts$pars
      tree_depth_t[i] <- nuts$tree_depth
      #eps_t[i] <- nuts$pars$eps_bar
      reg_params_postsamp_t[i+1, ] <- reg_params_current
      Yplus_postsamp_t[i + 1] <- Yplus_current
     # print(paste0(i, '--', nuts$tree_depth))
      
      Nr_sum_current <- sum(pop_dat$N*exp(pop_dat$U*reg_params_current[4] + (1-pop_dat$U)*reg_params_current[5] + reg_params_current[5+pop_dat$A]))
      Nr_obs_current <- data$N*exp(data$U*reg_params_current[4] + (1-data$U)*reg_params_current[5] + reg_params_current[5+data$A])
      d_current <- exp(reg_params_current[3])

     # draw Yplus ----
      
      # domain <- round(0.5*Yplus_current):(2*Yplus_current)
      # f <- log_target_cond_5.2_Yplus(domain,pop_dat,Y_current,Nr_sum_current,Nr_obs_current,
      #                                d_current,logit_r_hat,logit_r_hat_var)
      # domain <- domain[f>-Inf]
      # f <- f[f>-Inf]
      # del_f <- f-mean(f) #normalize
      # #restrict domain to densities >0
      # zero_dens <- which(exp(del_f)==0)
      # while(length(zero_dens)>0){
      #   domain <- domain[-zero_dens]
      #   f <- f[-zero_dens]
      #   del_f <- f-mean(f)
      #   zero_dens <- which(exp(del_f)==0)
      # }
      # while(sum(exp(del_f)==Inf)>0){
      #   domain <- domain[del_f>0]
      #   f <- f[del_f>0]
      #   del_f <- f-mean(f) #normalize
      # }
      # 
      # target_density <- (exp(del_f)/sum(exp(del_f)))
      # target_F <- cumsum(target_density)
      # 
      # #take draw from inverse CDF
      # U <- runif(1)
      # Yplus_current <- (min(domain[target_F>U]))
      # Yplus_postsamp_t[i + 1] <- Yplus_current
      
     # draw Ys ----
      for(k in 1:nrow(data)){
        # derive CDF
        domain <- data$Z[k]:round(0.5*data$N[k])
        f <- log_target_cond_5.2_y(domain,data,Yplus_current,pop_dat,
                                   Y_current,k,Nr_sum_current,Nr_obs_current,d_current)
        domain <- domain[f>-Inf]
        f <- f[f>-Inf]
        del_f <- f-mean(f) #normalize
        #restrict domain to densities >0
        zero_dens <- which(exp(del_f)==0)
        while(length(zero_dens)>0){
          domain <- domain[-zero_dens]
          f <- f[-zero_dens]
          del_f <- f-mean(f)
          zero_dens <- which(exp(del_f)==0)
        }
        
        target_density <- (exp(del_f)/sum(exp(del_f)))
        target_F <- cumsum(target_density)
        #take draw from inverse CDF
        U <- runif(1)
        Y_current[k] <- min(domain[target_F>U])
        if(is.na(Y_current[k])){
          stop()
        }
      }
      Y_postsamp_t[i + 1,] <- Y_current
      
    } 
    
    results <- list(reg_params_t=reg_params_postsamp_t, 
                    tree_depth_t=tree_depth_t,
                   # eps_t=eps_t,
                    Y_t=Y_postsamp_t,
                    Yplus_t=Yplus_postsamp_t)

    return(results)
  },mc.cores = chains)
  #Sys.time() - start.time
  
  #combine chains
  b_postsamp <- reg_params_postsamp <- Y_postsamp <- NULL
  for(j in 1:chains){
    b_postsamp <- rbind(b_postsamp,parallel_results[[j]]$reg_params_t[,6:(5+n_areas)])
    reg_params_postsamp <- rbind(reg_params_postsamp,parallel_results[[j]]$reg_params_t)
    Y_postsamp <- rbind(Y_postsamp,parallel_results[[j]]$Y_t)
  }
  tau_postsamp <- unlist(lapply(1:chains,function(x){exp(parallel_results[[x]]$reg_params_t[,1])}))
  phi_postsamp <- unlist(lapply(1:chains,function(x){expit(parallel_results[[x]]$reg_params_t[,2])}))
  d_postsamp <- unlist(lapply(1:chains,function(x){exp(parallel_results[[x]]$reg_params_t[,3])}))
  alphaU_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,4]}))
  alphaR_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,5]}))
  Yplus_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$Yplus_t}))
  tree_depth <- lapply(1:chains,function(x){parallel_results[[x]]$tree_depth_t})
  #eps <- lapply(1:chains,function(x){parallel_results[[x]]$eps_t})
  names(tau_postsamp) <-names(phi_postsamp) <- names(d_postsamp) <- names(alphaU_postsamp) <- names(alphaR_postsamp) <- NULL
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp,alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  
  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              tau=tau_postsamp, phi=phi_postsamp, 
              reg_params=reg_params_postsamp,
              d=d_postsamp, 
              Y=Y_postsamp,
              Yplus=Yplus_postsamp,
              #eps=eps,
              tree_depth=tree_depth))
}

start.time <- Sys.time()
NUTS_5.2a <- postsamp_Alg5.2a_NUTS(eps0=0.01,
                                  # kappa = 0.9, gamma = 0.35, t0=5,M_adapt=1,
                                   #M_diag = 1/c(1,2,1,rep(0.5,30),200),
                                   M_diag = 1/c(1,2,1,rep(0.5,30),1000),
                                   logit_r_hat=natl.dir$logit.est,logit_r_hat_var=natl.dir$var.est,
                                 data_list=data_list, n_iter=10000,warmup=50,chains=4)
Sys.time() - start.time

iter.id <- c(1001:5000,11001:15000,21001:25000,31001:35000)
plot(exp(NUTS_5.2a$reg_params[iter.id,1]),type='l')
plot(expit(NUTS_5.2a$reg_params[iter.id,2]),type='l')
plot(exp(NUTS_5.2a$reg_params[iter.id,3]),type='l')
plot(NUTS_5.2a$reg_params[iter.id,4],type='l')
plot(NUTS_5.2a$reg_params[iter.id,5],type='l')
plot(NUTS_5.2a$reg_params[iter.id,6],type='l')
plot(NUTS_5.2a$reg_params[iter.id,7],type='l')
plot(NUTS_5.2a$reg_params[iter.id,21],type='l')
plot(NUTS_5.2a$reg_params[iter.id,24],type='l')
plot(NUTS_5.2a$Yplus[iter.id],type='l')

plot(NUTS_5.2a$reg_params[,1],type='l')
plot(NUTS_5.2a$reg_params[,2],type='l')
plot(NUTS_5.2a$reg_params[,3],type='l')
plot(NUTS_5.2a$reg_params[,4],type='l')
plot(NUTS_5.2a$reg_params[,5],type='l')
plot(NUTS_5.2a$reg_params[,6],type='l')
plot(NUTS_5.2a$Yplus,type='l')

save(NUTS_5.2a,file="Handcoded MCMC Simulations/Alg 5.2a BYM2 NUTS 240223, MWI 2011-2015, 4 chains of 10k iterations posterior dist.rda")

each_chain <- length(NUTS_5.2a$Yplus)/chains
  plot(cumsum(NUTS_5.2a$Yplus[1:each_chain])/(1:each_chain),type='l',col=1,ylab='',
       ylim=c(min(NUTS_5.2a$Yplus),max(NUTS_5.2a$Yplus)),main = paste0('Running mean of Yplus_',k))
  for(c in 2:chains){
    lines(cumsum(NUTS_5.2a$Yplus[((c-1)*each_chain+1):(c*each_chain)])/(1:each_chain),col=c)
  }
  
  
plot(cumsum(NUTS_5.2a$alphaU[1:each_chain])/(1:each_chain),type='l',col=1,ylab='',
       ylim=c(min(NUTS_5.2a$alphaU),max(NUTS_5.2a$alphaU)),main = paste0('Running mean of Yplus_',k))
  for(c in 2:chains){
    lines(cumsum(NUTS_5.2a$alphaU[((c-1)*each_chain+1):(c*each_chain)])/(1:each_chain),col=c)
  }

plot(cumsum(NUTS_5.2a$alphaR[1:each_chain])/(1:each_chain),type='l',col=1,ylab='',
     ylim=c(min(NUTS_5.2a$alphaR),max(NUTS_5.2a$alphaR)),main = paste0('Running mean of Yplus_',k))
for(c in 2:chains){
  lines(cumsum(NUTS_5.2a$alphaR[((c-1)*each_chain+1):(c*each_chain)])/(1:each_chain),col=c)
}


#######################################################################
#######################################################################
# Alg 5.3a with NUTS ---------
tau_init=40
phi_init=0.7
d_init = 1
alphaU_init = -3.5
alphaR_init= -3.2 #initial values
alpha_tau = 0.01
U_tau = 1
alpha_phi=0.5
beta_phi=1
mu_U = -3.5
sd_U = 3
mu_R = -3.2
sd_R = 3
mu_logd = 0
sd_logd = 1
eps_adapt = 100
DA=F
#hyperpriors  


eps0 = 0.05
M_init = solve(Sigma)
logit_r_hat = dir.est$logit.est
logit_r_hat_var = dir.est$var.est
n_iter = 10
chains=4

postsamp_Alg5.3a_NUTS_MCMC <- function(tau_init=40, phi_init=0.7, d_init = 1, alphaU_init = -3.5, alphaR_init= -3.2, #initial values
                                       alpha_tau = 0.01,U_tau = 1,alpha_phi=0.5,beta_phi=1,mu_U = -3.5,sd_U = 3,mu_R = -3.2,sd_R = 3,mu_logd = 0,sd_logd = 1,#hyperpriors  
                                       eps0, M_init,eps_adapt = 100, DA=F,  #tuning params for regression param,d, and Yplus proposals
                                       logit_r_hat,logit_r_hat_var,data_list, 
                                       n_iter, chains=4, warmup1=0,  #how long to wait before updating M
                                       warmup2=n_iter){ #when to stop updating M
  
  data <- data_list$obs_dat
  data <- data[order(data$A2),]
  pop_dat <- data_list$pop_strata_dat
  #make data tables
  data_table <- as.data.table(data)
  pop_dt <- as.data.table(pop_dat)
  
  admin1_pop_dat <- pop_dat %>% group_by(A1) %>% summarise(N = sum(N))
  Yplus_init <- round(admin1_pop_dat$N*expit(logit_r_hat))
  W_inv_init <- solve(W_fn(phi_init,data_list$Q_scaled_inv))
  
  n_admin1 <- length(unique(data_list$pop_strata_dat$A1))
  n_admin2 <- nrow(data_list$Q_scaled_inv)
  
  cons_params <- list(pop_dt=pop_dt,Q_inv=data_list$Q_scaled_inv,
                     logit_r_hat = logit_r_hat, logit_r_hat_var=logit_r_hat_var,
                     mu_logd=mu_logd,sd_logd=sd_logd,mu_U=mu_U,mu_R=mu_R,sd_U=sd_U,sd_R=sd_R,
                     alpha_tau=alpha_tau,U_tau=U_tau,alpha_phi=alpha_phi,beta_phi=beta_phi)
  
  parallel_results <- mclapply(1:chains, function(j){
    
    b_init <- rnorm(n_admin2,0,1/sqrt(tau_init))
    reg_params_init <-  c(log(tau_init),logitlink(phi_init),log(d_init),alphaU_init,alphaR_init,b_init)
      
    Y_init <- sapply(1:n_admin1,function(i){
        
        subdat <- data %>% filter(A1==i)
        subpop_dat <- pop_dat %>% filter(A1==i)
        
        Nr_obs_init <- subdat$N*exp(subdat$U*reg_params_init[4] + (1-subdat$U)*reg_params_init[5] + b_init[subdat$A2])
        Nr_sum_init <- sum(subpop_dat$N*exp(subpop_dat$U*reg_params_init[4] + (1-subpop_dat$U)*reg_params_init[5] + b_init[subpop_dat$A2]))
        
        Y_init_t <- rmultinom(1,Yplus_init[i],prob=c(Nr_obs_init,Nr_sum_init-sum(Nr_obs_init)))
        Y_init_t <- Y_init_t[-length(Y_init_t)]
        Y_init_t[Y_init_t<subdat$Z] <- subdat$Z[Y_init_t<subdat$Z]
        
        return(Y_init_t)
      })
    data_table[, ':='(Y=unlist(Y_init))]
    Y_sum_obs_init <- (data_table[, sum(Y),by=A1][order(A1)])$V1
    
    ## initialize data frames 
    reg_params_postsamp_t <- rbind(reg_params_init, matrix(NA,n_iter-1,length(reg_params_init)))
    Yplus_postsamp_t <- rbind(Yplus_init, matrix(NA,n_iter-1,n_admin1))
    tree_depth_t <- rep(NA,n_iter)
    
    # set current states
    reg_params_current <- reg_params_init
    Y_current <- Y_init
    Yplus_current <- Yplus_init
    Y_sum_obs_current <- Y_sum_obs_init
    par_list <- list(eps=eps0, M=M_init, eps_adapt=eps_adapt)
    
    for(i in 2:n_iter){
      
      # draw regression parameters using NUTS -----
      cons_params$data_table <- data_table
      cons_params$Y_sum_obs <- Y_sum_obs_current
      cons_params$Yplus <- Yplus_current
      
      # if(i>warmup1 & i<warmup2){
      #   new_cov <- cov(reg_params_postsamp_t[1:(i-1),])
      #   if(is.positive.semidefinite(new_cov)){
      #     par_list$M <- solve(new_cov)
      #   }
      # }
      
      nuts <- NUTS_one_step(theta = reg_params_current, iter = i, f = log_target_5.3a_reg_params_4NUTS, grad_f=grad_log_target_5.3a_4NUTS, cons_params, 
                            par_list = par_list, max_treedepth = 10, verbose = F, DA=DA)
      
      reg_params_current <- nuts$theta
      d_current <- exp(reg_params_current[3])
      par_list <- nuts$pars
      tree_depth_t[i] <- nuts$tree_depth
      print(paste0(i, '--', nuts$tree_depth))
      
      Nr_sum_current <- sapply(1:n_admin1,function(i){
        subpop_dat <- pop_dat %>% filter(A1==i)
        Nr_sum_current_t <- sum(subpop_dat$N*exp(subpop_dat$U*reg_params_current[4] + (1-subpop_dat$U)*reg_params_current[5] + reg_params_current[5+subpop_dat$A2]))
        return(Nr_sum_current_t)
      })
      
      Nr_obs_current <- sapply(1:n_admin1,function(i){
        subdat <- data %>% filter(A1==i)
        Nr_obs_t <- subdat$N*exp(subdat$U*reg_params_current[4] + (1-subdat$U)*reg_params_current[5] + reg_params_current[5+subdat$A2])
        return(Nr_obs_t)
      })
      
      # draw proposal for Yplus ----
        for(area in 1:n_admin1){
          domain <- round(0.5*Yplus_current[area]):(2*Yplus_current[area])
          f <- log_target_cond_5.3_Yplus(domain,pop_dat[pop_dat$A1==area,],Y_current[[area]],Nr_sum_current[area],Nr_obs_current[[area]],
                                         d_current,logit_r_hat[area],logit_r_hat_var[area])
          
          domain <- domain[f>-Inf]
          f <- f[f>-Inf]
          del_f <- f-mean(f) #normalize
          #restrict domain to densities >0
          zero_dens <- which(exp(del_f)==0)
          while(length(zero_dens)>0){
            domain <- domain[-zero_dens]
            f <- f[-zero_dens]
            del_f <- f-mean(f)
            zero_dens <- which(exp(del_f)==0)
          }
          
          target_density <- (exp(del_f)/sum(exp(del_f)))
          target_F <- cumsum(target_density)
          
          #take draw from inverse CDF
          U <- runif(1)
          Yplus_current[area] <- (min(domain[target_F>U]))
        }
      
      
      # draw new Ys from condtional distribution (ITS) ----
      for(k in 1:nrow(data)){
        area_id <- data$A1[k]
        within_area_id <- which(which(data$A1==area_id)==k)
        
        # derive CDF
        domain <- data$Z[k]:round(0.5*data$N[k])
        f <- log_target_cond_5.3_y(domain,data[data$A1==area_id,],Yplus_current[area_id],pop_dat[pop_dat$A1==area_id,],
                                   Y_current[[area_id]],within_area_id,Nr_sum_current[area_id],Nr_obs_current[[area_id]],d_current)
        domain <- domain[f>-Inf]
        f <- f[f>-Inf]
        del_f <- f-mean(f) #normalize
        #restrict domain to densities >0
        zero_dens <- which(exp(del_f)==0)
        while(length(zero_dens)>0){
          domain <- domain[-zero_dens]
          f <- f[-zero_dens]
          del_f <- f-mean(f)
          zero_dens <- which(exp(del_f)==0)
        }
        
        target_density <- (exp(del_f)/sum(exp(del_f)))
        target_F <- cumsum(target_density)
        #take draw from inverse CDF
        U <- runif(1)
        Y_current[[area_id]][within_area_id] <- min(domain[target_F>U])
      }
      
      # record current state ----
      reg_params_postsamp_t[i,] <- reg_params_current
      Yplus_postsamp_t[i,] <- Yplus_current
      data_table[, ':='(Y=unlist(Y_current))]
      Y_sum_obs_current <- (data_table[, sum(Y),by=A1][order(A1)])$V1
      
    } 
    
    results <- list(reg_params_t=reg_params_postsamp_t, 
                    Yplus_t=Yplus_postsamp_t,
                    tree_depth_t = tree_depth_t)
    
    return(results)
  },mc.cores = chains)
  
  #combine chains
  b_postsamp <- reg_params_postsamp <- Yplus_postsamp <- NULL
  for(j in 1:chains){
    b_postsamp <- rbind(b_postsamp,parallel_results[[j]]$reg_params_t[,6:(5+n_admin2)])
    reg_params_postsamp <- rbind(reg_params_postsamp,parallel_results[[j]]$reg_params_t)
    Yplus_postsamp <- rbind(Yplus_postsamp,parallel_results[[j]]$Yplus_t)
  }
  logtau_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,1]}))
  logitphi_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,2]}))
  logd_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,3]}))
  alphaU_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,4]}))
  alphaR_postsamp <- unlist(lapply(1:chains,function(x){parallel_results[[x]]$reg_params_t[,5]}))
  
  eta_postsamp <- cbind(alphaU_postsamp + b_postsamp,alphaR_postsamp + b_postsamp)
  r_postsamp <- exp(eta_postsamp)
  
  tree_depth <- lapply(1:chains,function(x){parallel_results[[x]]$tree_depth_t})
  
  return(list(alphaU = alphaU_postsamp, alphaR = alphaR_postsamp,
              b=b_postsamp, eta=eta_postsamp, r=r_postsamp, 
              logtau=logtau_postsamp, 
              logitphi=logitphi_postsamp,
              reg_params=reg_params_postsamp,
              logd=logd_postsamp, 
              Yplus=Yplus_postsamp,
              tree_depth = tree_depth))
}


