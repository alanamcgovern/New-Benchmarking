# required packages
suppressMessages({
 library(data.table) # loaded in cluster
  library(LaplacesDemon) # loaded in cluster
  library(writexl) # loaded in cluster
  library(stats)
})


leapfrog_step = function(theta, r, eps, M_inv, f, grad_f, cons_params){
  r_tilde <- r + 0.5 * eps * grad_f(theta,cons_params)
  theta_tilde <- theta + eps * M_inv%*%r_tilde
  r_tilde <- r_tilde + 0.5 * eps * grad_f(theta_tilde,cons_params)
  list(theta = theta_tilde, r = r_tilde)
}

find_reasonable_epsilon = function(theta, f, grad_f, cons_params, M, eps = .1, verbose = TRUE){
  r <- as.vector(rmvnorm(1, rep(0,length(theta)), M))
  proposed <- leapfrog_step(theta, r, eps, M_inv, f, grad_f, cons_params)
  log_ratio <- f(params = proposed$theta,cons_params = cons_params) - 0.5*sum(proposed$r**2 / M) - 
    f(params = theta,cons_params = cons_params) - 0.5*sum(r**2 / M)
  alpha <- as.numeric(ifelse(exp(log_ratio) > 0.5, 1, -1))
  if(is.nan(alpha)) alpha <- -1
  count <- 1
  while(is.nan(log_ratio) || alpha * log_ratio > (-alpha)*log(2)){
    eps <- 2**alpha * eps
    proposed <- leapfrog_step(theta, r, eps, M_inv, f, grad_f, cons_params)
    log_ratio <- f(params = proposed$theta,cons_params = cons_params) -  0.5 * as.numeric(t(proposed$r)%*%M_inv%*%proposed$r) - 
      f(params = theta,cons_params = cons_params) -  0.5 * as.numeric(t(r)%*%M_inv%*%r)
    count <- count + 1
    if(count > 100) {
      stop("Could not find reasonable epsilon in 100 iterations!")
    }
  }
  if(verbose) message("Reasonable epsilon = ", eps, " found after ", count, " steps")
  eps
}

check_NUTS_scaled = function(s, theta_plus, theta_minus, r_plus, r_minus, M){
  if(is.na(s)) return(0)
  #condition1 <- sum(M%*%(theta_plus - theta_minus)*r_minus) >= 0
  # condition2 <- sum(M%*%(theta_plus - theta_minus)*r_plus) >= 0
  condition1 <- crossprod(M%*%(theta_plus - theta_minus), r_minus) >= 0
  condition2 <- crossprod(M%*%(theta_plus - theta_minus), r_plus) >= 0
  s && condition1 && condition2
}

build_tree = function(theta, r, u, v, j, eps, theta0, r0, f, grad_f, cons_params, M, M_inv, Delta_max = 1000){
  if(j == 0){
    proposed <- leapfrog_step(theta, r, v*eps, M_inv, f, grad_f, cons_params)
    theta <- proposed$theta
    r <- proposed$r
    log_prob <- f(params = theta,cons_params = cons_params) - 0.5 * as.numeric(t(r)%*%M_inv%*%r)
    log_prob0 <- f(params = theta0,cons_params = cons_params) - 0.5*as.numeric(t(r0)%*%M_inv%*%r0)
    n <- (log(u) <= log_prob)
    s <- (log(u) < Delta_max + log_prob)
    alpha <- min(1, exp(log_prob - log_prob0))
    if(is.nan(alpha)) stop()
    if(is.na(s) || is.nan(s)){
      s <- 0
    }
    if(is.na(n) || is.nan(n)){
      n <- 0
    }
    return(list(theta_minus=theta, theta_plus=theta, theta=theta, r_minus=r,
                r_plus=r, s=s, n=n, alpha=alpha, n_alpha=1))
  } else{
    obj0 <- build_tree(theta, r, u, v, j-1, eps, theta0, r0, f, grad_f, cons_params, M, M_inv)
    theta_minus <- obj0$theta_minus
    r_minus <- obj0$r_minus
    theta_plus <- obj0$theta_plus
    r_plus <- obj0$r_plus
    theta <- obj0$theta
    if(obj0$s == 1){
      if(v == -1){
        obj1 <- build_tree(obj0$theta_minus, obj0$r_minus, u, v, j-1, eps, theta0, r0, f, grad_f, cons_params, M, M_inv)
        theta_minus <- obj1$theta_minus
        r_minus <- obj1$r_minus
      } else{
        obj1 <- build_tree(obj0$theta_plus, obj0$r_plus, u, v, j-1, eps, theta0, r0, f, grad_f, cons_params, M, M_inv)
        theta_plus <- obj1$theta_plus
        r_plus <- obj1$r_plus
      }
      n <- obj0$n + obj1$n
      if(n != 0){
        prob <- obj1$n / n
        if(runif(1) < prob){
          theta <- obj1$theta
        }
      }
      s <- check_NUTS_scaled(obj1$s, theta_plus, theta_minus, r_plus, r_minus, M)
      alpha <- obj0$alpha + obj1$alpha
      n_alpha <- obj0$n_alpha + obj1$n_alpha
      
    } else{
      n <- obj0$n
      s <- obj0$s
      alpha <- obj0$alpha
      n_alpha <- obj0$n_alpha
    }
    if(is.na(s) || is.nan(s)){
      s <- 0
    }
    if(is.na(n) || is.nan(n)){
      n <- 0
    }
    return(list(theta_minus=theta_minus, theta_plus=theta_plus, theta=theta,
                r_minus=r_minus, r_plus=r_plus, s=s, n=n, alpha=alpha, n_alpha=n_alpha))
  }
}

NUTS_one_step <- function(theta, iter, f, grad_f, cons_params, par_list, max_treedepth = 10, verbose = TRUE,
                          DA=F, 
                          kappa = 0.75, t0 = 10, # stabilizes initial iterations
                          gamma = 0.05, #controls shrinkage towards mu
                          delta = 0.5 ){ # target acceptance prob
  
  if(DA){
    eps_adapt <- par_list$eps_adapt
  }
  
  M <- par_list$M
  M_inv <- solve(M)
  
  if(!DA){
    eps <- par_list$eps
  }else{
    if(iter == 2){
      eps <- find_reasonable_epsilon(theta, f, grad_f, cons_params, M, eps = eps, verbose = verbose)
      mu <- log(10*eps)
      H <- 0
      eps_bar <- 1
    } else{
      eps <- par_list$eps
      eps_bar <- par_list$eps_bar
      H <- par_list$H
      mu <- par_list$mu
    }
  }
  
  chol_M <- t(chol(M))
  r0 <- chol_M%*%rnorm(length(theta),0,1)
  
  suppressWarnings({
    u <- runif(1, 0, exp(f(params = theta,cons_params = cons_params) - 0.5 * t(r0)%*%M_inv%*%r0))
  })
  
  if(is.nan(u)){
   # warning("NUTS: sampled slice u is NaN")
    u <- runif(1, 0, 1e5)
  }
  theta_0 <- theta
  theta_minus <- theta
  theta_plus <- theta
  r_minus <- r0
  r_plus <- r0
  j=0
  n=1
  s=1
  
  if(DA){
    if(iter > eps_adapt){
      eps <- runif(1, 0.9*eps_bar, 1.1*eps_bar)
    }
  }
  
  while(s == 1){
    # choose direction {-1, 1}
    direction <- sample(c(-1, 1), 1)
    if(direction == -1){
      temp <- build_tree(theta_minus, r_minus, u, direction, j, eps, theta_0, r0, f, grad_f, cons_params, M, M_inv)
      theta_minus <- temp$theta_minus
      r_minus <- temp$r_minus
    } else{
      temp <- build_tree(theta_plus, r_plus, u, direction, j, eps, theta_0, r0, f, grad_f, cons_params, M, M_inv)
      theta_plus <- temp$theta_plus
      r_plus <- temp$r_plus
    }
    if(is.nan(temp$s)) temp$s <- 0
    if(temp$s == 1){
      if(runif(1) < temp$n / n){
        theta <- temp$theta
      }
    }
    n <- n + temp$n
    s <- check_NUTS_scaled(temp$s, theta_plus, theta_minus, r_plus, r_minus, M)
    j <- j + 1
    if(j > max_treedepth){
      warning("NUTS: Reached max tree depth")
      break
    }
  }
  
  if(DA){
    if(iter <= eps_adapt){
      H <- (1 - 1/(iter + t0))*H + 1/(iter + t0) * (delta - temp$alpha / temp$n_alpha)
      log_eps <- mu - sqrt(iter)/gamma * H
      eps_bar <- exp(iter**(-kappa) * log_eps + (1 - iter**(-kappa)) * log(eps_bar))
      eps <- exp(log_eps)
    } else{
      eps <- eps_bar
    }
  }
  
  
  if(!DA){
    return(list(theta = theta, tree_depth=j-1,
                pars = list(eps = eps, M = M)))
  }else{
    return(list(theta = theta, tree_depth=j-1,
                pars = list(eps = eps, eps_bar = eps_bar, H = H, mu = mu, eps_adapt = eps_adapt, M = M)))
  }
  
}

W_fn <- function(phi,Q_inv){
  return(diag(1-phi,nrow(Q_inv)) + phi*Q_inv)
}

log_target_cond_5.3_Yplus <- function(Yplus,pop_dat,Y,Nr_sum,Nr_obs,d,logit_r_hat,logit_r_hat_var){
  term1 <- -1/(2*logit_r_hat_var)*(LaplacesDemon::logit(Yplus/sum(pop_dat$N)) - logit_r_hat)^2
  term2 <- Yplus*log(d/(1+d)) 
  term3 <- lgamma(Yplus - sum(Y) + (1/d)*(Nr_sum-sum(Nr_obs))) - lgamma(Yplus - sum(Y) + 1)
  
  return(term1 + term2 + term3)
}

log_target_cond_5.3_y <- function(y,data,Yplus,Y,i,Nr_sum,Nr_obs,d){
  term1 <- lgamma(Yplus - sum(Y[-i]) - y + 1/d*(Nr_sum - sum(Nr_obs))) - lgamma(Yplus - sum(Y[-i]) -y + 1)
  term2 <- lgamma(data$N[i]-y+1) - lgamma(y-data$Z[i]+1) - lgamma(data$N[i]-y-data$n[i]+data$Z[i] +1)
  term3 <- lgamma(y + 1/d*Nr_obs[i])
  
  return(term1 + term2 + term3)
}

log_target_5.4_reg_params <- function(params,cons_params){
  tau <- exp(params[1])
  phi <- 1/(1+exp(-params[2]))
  W <- W_fn(phi,cons_params$Q_inv)
  W_inv <- solve(W)
  d <- exp(params[3])
  alphaU <- params[4]
  alphaR <- params[5]
  beta <- c(0,params[6:(4+cons_params$n_admin1)])
  b <- params[(5+cons_params$n_admin1):(length(params))]
  
  # add Nr_obs to data
  cons_params$data_table[, ':='(Nr_obs=N*exp(U*alphaU + (1-U)*alphaR + beta[A1] + b[A2]))]
  
  # get grouped sums faster w data.table package
  N_sum <- (cons_params$pop_dt[, sum(N), by=A1][order(A1)])$V1
  Nr_sum <- (cons_params$pop_dt[, sum(N*exp(U*alphaU + (1-U)*alphaR + beta[A1] + b[A2])), by=A1][order(A1)])$V1
  Nr_sum_obs <- (cons_params$data_table[, sum(Nr_obs),by=A1][order(A1)])$V1
  
  out <- sum(cons_params$Yplus)*log(d/(1+d)) - sum(log(1+d)*Nr_sum/d) +
    sum(lgamma(cons_params$Yplus - cons_params$Y_sum_obs + (1/d)*(Nr_sum-Nr_sum_obs))- lgamma((1/d)*(Nr_sum-Nr_sum_obs))) +
    sum(lgamma(cons_params$data_table$Y+1/d*cons_params$data_table$Nr_obs) - lgamma(1/d*cons_params$data_table$Nr_obs)) +
    -0.5*tau*t(b)%*%W_inv%*%b -
    (alphaU-cons_params$mu_U)^2/(2*cons_params$sd_U^2) - (alphaR-cons_params$mu_R)^2/(2*cons_params$sd_R^2) -
    sum((beta - cons_params$mu_beta)^2/(2*cons_params$sd_beta^2)) + 
    log(d) - cons_params$lambda_d*d +
    (length(b)-1)/2*log(tau) + log(cons_params$alpha_tau)/cons_params$U_tau*tau^(-0.5) - 0.5*logdet(W) + cons_params$alpha_phi*(log(phi)) + cons_params$beta_phi*log(1-phi)
  
  return(out)
}

grad_log_target_5.4 <- function(params,cons_params){
  
  tau <- exp(params[1])
  phi <- 1/(1+exp(-params[2]))
  W <- W_fn(phi,cons_params$Q_inv)
  W_inv <- solve(W)
  d <- exp(params[3])
  alphaU <- params[4]
  alphaR <- params[5]
  beta <- c(0,params[6:(4+cons_params$n_admin1)])
  b <- params[(5+cons_params$n_admin1):(length(params))]
  
  # add Nr_obs to data
  cons_params$data_table[, ':='(Nr_obs_urban=U*N*exp(U*alphaU + beta[A1]+ b[A2]))]
  cons_params$data_table[, ':='(Nr_obs_rural=(1-U)*N*exp((1-U)*alphaR + beta[A1] + b[A2]))]
  Nr_sum_obs_urban <- (cons_params$data_table[, sum(Nr_obs_urban),by=A1][order(A1)])$V1
  Nr_sum_obs_rural <- (cons_params$data_table[, sum(Nr_obs_rural),by=A1][order(A1)])$V1
  Nr_sum_obs <- (cons_params$data_table[, sum(Nr_obs_urban+Nr_obs_rural),by=A1][order(A1)])$V1
  
  # get grouped sums faster w data.table package
  N_sum <- (cons_params$pop_dt[, sum(N), by=A1][order(A1)])$V1
  
  Nr_sum_urban <- (cons_params$pop_dt[U==1, sum(N*exp(U*alphaU + (1-U)*alphaR + beta[A1] + b[A2])), by=A1][order(A1)])$V1
  Nr_sum_rural <- (cons_params$pop_dt[U==0, sum(N*exp(U*alphaU + (1-U)*alphaR + beta[A1] + b[A2])), by=A1][order(A1)])$V1
  Nr_sum <- (cons_params$pop_dt[, sum(N*exp(U*alphaU + (1-U)*alphaR + beta[A1] + b[A2])), by=A1][order(A1)])$V1
  
  Nr_sum_area <- cons_params$pop_dt[,sum(N*exp(U*alphaU + (1-U)*alphaR + beta[A1] + b[A2])), by=A2][order(A2)]$V1
  Nr_sum_obs_area <- rbind(cons_params$data_table[, sum(Nr_obs_urban+Nr_obs_rural),by=A2][order(A2)],
                           list(cons_params$missing_adm2_id,rep(0,length(cons_params$missing_adm2_id))))[order(A2)]$V1 # add in 0s for unobserved areas
  
  
  # digamma terms to be included in gradients
  digamma_1 <- digamma(cons_params$Yplus - cons_params$Y_sum_obs + 1/d*(Nr_sum - Nr_sum_obs)) - digamma(1/d*(Nr_sum - Nr_sum_obs))
  digamma_2 <- digamma(cons_params$data_table$Y + 1/d*(cons_params$data_table$Nr_obs_urban+cons_params$data_table$Nr_obs_rural)) - 
    digamma(1/d*(cons_params$data_table$Nr_obs_urban+cons_params$data_table$Nr_obs_rural))
  cons_params$data_table$digamma_2 <- digamma_2
  
  
  grad_logtau <- 0.5*(length(b)-1) - 0.5*log(cons_params$alpha_tau)/cons_params$U_tau*tau^(-0.5) - 0.5*tau*t(b)%*%W_inv%*%b
  
  grad_logitphi <- -0.5*phi*(1-phi)*tr(W_inv%*%(cons_params$Q_inv - diag(1,length(b)))) + 
    0.5*tau*phi*(1-phi)*t(b)%*%W_inv%*%(cons_params$Q_inv - diag(1,length(b)))%*%W_inv%*%b +
    cons_params$alpha_phi*(1-phi) - cons_params$beta_phi*phi
  
  grad_logd <- sum(cons_params$Yplus)*(1/(1+d)) - sum(Nr_sum)*(1/(1+d)-log(1+d)/d) - (1/d)*sum((Nr_sum - Nr_sum_obs)*digamma_1) +
    (1/d)*sum((cons_params$data_table$Nr_obs_urban + cons_params$data_table$Nr_obs_rural)*(digamma((cons_params$data_table$Nr_obs_urban + cons_params$data_table$Nr_obs_rural)/d) - 
                                                                                             digamma(cons_params$data_table$Y + (cons_params$data_table$Nr_obs_urban + cons_params$data_table$Nr_obs_rural)/d))) +
    1 - cons_params$lambda_d*d
  
  grad_alphaU <- sum(-log(1+d)/d*Nr_sum_urban + 1/d*(Nr_sum_urban - Nr_sum_obs_urban)*digamma_1) +
    1/d*sum(cons_params$data_table$Nr_obs_urban*digamma_2) - 
    (1/cons_params$sd_U^2)*(alphaU - cons_params$mu_U)
  
  grad_alphaR <- sum(-log(1+d)/d*Nr_sum_rural + 1/d*(Nr_sum_rural - Nr_sum_obs_rural)*digamma_1) +
    1/d*sum(cons_params$data_table$Nr_obs_rural*digamma_2) - 
    (1/cons_params$sd_R^2)*(alphaR - cons_params$mu_R)
  
  grad_beta <- -log(1+d)/d*Nr_sum + (1/d)*(Nr_sum - Nr_sum_obs)*digamma_1 +
    1/d*(cons_params$data_table[, sum(digamma_2*(Nr_obs_urban+Nr_obs_rural)),by=A1][order(A1)])$V1 - 
    (1/cons_params$sd_beta^2)*(beta - cons_params$mu_beta)
  grad_beta <- grad_beta[2:cons_params$n_admin1]
  
  grad_b <- -log(1+d)/d*Nr_sum_area +
    1/d*(Nr_sum_area - Nr_sum_obs_area)*digamma_1[unique(cons_params$pop_dt[order(A2),c('A1','A2')])$A1] +
    1/d*(rbind(cons_params$data_table[, sum(digamma_2*(Nr_obs_urban+Nr_obs_rural)),by=A2],
               list(cons_params$missing_adm2_id,rep(0,length(cons_params$missing_adm2_id))))[order(A2)])$V1 - 
    tau*W_inv%*%b
  
  return(c(grad_logtau,grad_logitphi,grad_logd,grad_alphaU,grad_alphaR,grad_beta,grad_b))
}
# 
# alpha_tau = 0.01
# U_tau = 1
# alpha_phi=0.5
# beta_phi=0.5
# mu_U = -3.5
# sd_U = 3
# mu_R = -3.5
# sd_R = 3
# lambda_d = 1.5
# mu_beta = 0
# sd_beta=1#hyperpriors
# eps_adapt = 100
# DA=F

postsamp_Alg5.4_NUTS_MCMC_onechain <- function(alpha_tau = 0.01, U_tau = 1, alpha_phi=0.5, beta_phi=0.5, mu_U = -3.5, sd_U = 3, mu_R = -3.5, sd_R = 3, lambda_d = 1, 
                                               mu_beta = 0, sd_beta=1,#hyperpriors  
                                               eps0, M_init, eps_adapt = 100, DA=F,  #tuning params for regression param,d, and Yplus proposals
                                               logit_r_hat, logit_r_hat_var, data_list, 
                                               n_iter,
                                               add=F,current_state = NULL){ # adding more iterations to existing chain
  
  data <- data_list$obs_dat
  data <- data[order(data$A1),] ## this matters for when new Y's added to data table
  pop_dat <- data_list$pop_strata_dat
  
  #make data tables
  data_table <- as.data.table(data)
  pop_dt <- as.data.table(pop_dat)
  
  n_admin1 <- length(unique(pop_dat$A1))
  n_admin2 <- nrow(data_list$Q_scaled_inv)
  # are there any admin2 areas with no data?
  missing_adm2_id <- which(!(1:n_admin2 %in% unique(data$A2)))
  
  cons_params <- list(pop_dt=pop_dt,Q_inv=data_list$Q_scaled_inv, n_admin1 = n_admin1, n_admin2 = n_admin2, missing_adm2_id = missing_adm2_id,
                      lambda_d=lambda_d,mu_U=mu_U,mu_R=mu_R,sd_U=sd_U,sd_R=sd_R,mu_beta=mu_beta,sd_beta=sd_beta,
                      alpha_tau=alpha_tau,U_tau=U_tau,alpha_phi=alpha_phi,beta_phi=beta_phi)
  
  if(add){
    Yplus_init <- unlist(c(current_state[1,(6+n_admin1+3*n_admin2):(5+2*n_admin1+3*n_admin2)]))
    tau_init <- exp(current_state$logtau)
    phi_init <- 1/(1+exp(-1*current_state$logitphi))
    d_init <- exp(current_state$logd)
    alphaU_init <- current_state$alphaU
    alphaR_init <- current_state$alphaR
    beta_init <- unlist(c(current_state[1,6:(5+n_admin1)]))
    b_init <- unlist(c(current_state[1,(6+n_admin1):(5+n_admin1+n_admin2)]))
    
  }else{
    Yplus_init <- round(pop_dt[,sum(N),by=A1][order(A1)]$V1*(1/(1+exp(-logit_r_hat))))
    tau_init <- exp(rnorm(1,4,0.5))
    phi_init <- 1/(1+exp(-rnorm(1,0,0.25)))
    d_init <- exp(rnorm(1,-1,0.5))
    alphaU_init <- rnorm(1,-3.5,0.25)
    alphaR_init <- rnorm(1,-3.5,0.25)
    beta_init <- c(0,rnorm(n_admin1-1,0,0.1))
    b_init <- rnorm(n_admin2,0,1/sqrt(tau_init))
  }
  
  reg_params_init <-  c(log(tau_init),LaplacesDemon::logit(phi_init),log(d_init),alphaU_init,alphaR_init,beta_init[2:n_admin1],b_init)
  
  Y_init <- sapply(1:n_admin1,function(i){
    
    subdat <- data_table[A1==i,]
    subpop_dat <- pop_dt[A1==i,]

    Nr_obs_init <- subdat[,N*exp(U*reg_params_init[4] + (1-U)*reg_params_init[5] + beta_init[A1] + b_init[A2])]
    Nr_sum_init <- subpop_dat[,sum(N*exp(U*reg_params_init[4] + (1-U)*reg_params_init[5] + beta_init[A1] + b_init[A2]))]

    Y_init_t <- rmultinom(1,Yplus_init[i],prob=c(Nr_obs_init,Nr_sum_init-sum(Nr_obs_init)))
    Y_init_t <- Y_init_t[-length(Y_init_t)]
    Y_init_t[Y_init_t<subdat$Z] <- subdat$Z[Y_init_t<subdat$Z]
    
    return(Y_init_t)
  })
  data_table[, ':='(Y=unlist(Y_init))]
  Y_sum_obs_init <- (data_table[, sum(Y),by=A1][order(A1)])$V1
  
  ## initialize data frames 
  reg_params_postsamp <- rbind(reg_params_init, matrix(NA,n_iter-1,length(reg_params_init)))
  Yplus_postsamp <- rbind(Yplus_init, matrix(NA,n_iter-1,n_admin1))
  Y_postsamp <- sapply(1:n_admin1,function(i){rbind(Y_init[[i]],matrix(NA,n_iter-1,length(Y_init[[i]])))})
  tree_depth <- rep(NA,n_iter)
  tries <- matrix(NA,n_iter,n_admin1)
  colnames(tries) <- paste0('tries',1:n_admin1)
  
  # set current states
  reg_params_current <- reg_params_init
  Y_current <- Y_init
  Yplus_current <- Yplus_init
  Y_sum_obs_current <- Y_sum_obs_init
  par_list <- list(eps=eps0, M=M_init, eps_adapt=eps_adapt)
  
  data_table$rho <- data_table$wt/data_table[,by=A1,sum(wt*n)]$V1[data_table$A1]
  
  for(i in 2:n_iter){
    # check if got stuck during warmup
    if(i==5){
      if(nrow(unique.array(reg_params_postsamp))<4){
        # force an error
        return(NULL)
      }
    }
    # cat(i,'\n')
    # draw regression parameters using NUTS -----
    cons_params$data_table <- data_table
    cons_params$Y_sum_obs <- Y_sum_obs_current
    cons_params$Yplus <- Yplus_current

    nuts <- NUTS_one_step(theta = reg_params_current, iter = i, f = log_target_5.4_reg_params, grad_f=grad_log_target_5.4, cons_params,
                          par_list = par_list, max_treedepth = 10, verbose = F, DA=DA)

    reg_params_current <- nuts$theta
    par_list <- nuts$pars
    tree_depth[i] <- nuts$tree_depth
    #print(paste0(i, '--', nuts$tree_depth))

    d_current <- exp(reg_params_current[3])
    beta_current <- c(0,reg_params_current[6:(4+n_admin1)])

    Nr_sum_current <- (pop_dt[,sum(N*exp(U*reg_params_current[4] + (1-U)*reg_params_current[5] + beta_current[A1] + reg_params_current[4+n_admin1+A2])),by=A1][order(A1)])$V1

    Nr_obs_current <- sapply(1:n_admin1,function(i){
      Nr_obs_t <- data_table[A1==i,N*exp(U*reg_params_current[4] + (1-U)*reg_params_current[5] + beta_current[A1] + reg_params_current[4+n_admin1+A2])]
      return(Nr_obs_t)
    })

    # draw new Ys and Yplus from condtional distribution (ITS) ----
    if(i>50 & add==F){
    for(area_id in 1:n_admin1){
      Nr_sum <- Nr_sum_current[area_id]
      Nr_obs <- Nr_obs_current[[area_id]]
      r_hat <- 1/(1+exp(-logit_r_hat[area_id]))
      N <- data_table[A1==area_id,]$N
      n <- data_table[A1==area_id,]$n
      Z <- data_table[A1==area_id,]$Z
      rho <- data_table[A1==area_id,]$rho
      
      try <- T
      try_count <- 0
      while(try){
        Yplus <- Yplus_current[area_id]
        Y <- Y_current[[area_id]]
        rho_z_samples <- sapply(1:length(N),function(j){rho[j]*rhyper(2e3,Y[j],N[j]-Y[j],n[j])})
        
        domain <- round(0.5*Yplus):(2*Yplus)
        f <- log_target_cond_5.3_Yplus(domain,pop_dt[A1==area_id,],Y,Nr_sum,Nr_obs,
                                       d_current,logit_r_hat[area_id],logit_r_hat_var[area_id])
        
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
        if(is.na((min(domain[target_F>U]))))
          return(NULL)
        Yplus <- (min(domain[target_F>U]))
        
        for(within_area_id in 1:length(N)){
          
          N_k <- N[within_area_id]
          n_k <- n[within_area_id]
          Z_k <- Z[within_area_id]
          rho_k <- rho[within_area_id]
          
          domain <- Z_k:max(10,2*Y[within_area_id])
          
          term1 <- lgamma(N_k-domain + 1) + lgamma(domain + 1) -
            lgamma(N_k-domain-n_k+Z_k+1) - lgamma(domain-Z_k+1)
          
          term2 <- lgamma(domain + 1/d_current*Nr_obs[within_area_id]) - lgamma(domain + 1) +
            lgamma(Yplus - sum(Y[-within_area_id]) - domain + 1/d_current*(Nr_sum - sum(Nr_obs))) -
            lgamma(Yplus - sum(Y[-within_area_id]) - domain +1)
          
          term3 <- sapply(domain,function(y){
            rho_z_samples[,within_area_id] <- rho_k*rhyper(2e3,y,N_k-y,n_k)
            # if(min(rowSums(rho_z_samples))<0)
            #   stop('check')
            kde <- density(rowSums(rho_z_samples),from=r_hat,to=r_hat,n=1)
            return(min(-1*log(kde$y),-2))
          })
          
          if(sum(term3==0)==length(term3)){
            try_count <- try_count + 1
            if(try_count==50){
              stop(paste('We are stuck at area ',area_id,' at iteration ', i, '\n'))
            }
           # cat('Trying again for area', area_id,' (stopped at cluster',within_area_id,')\n')
            try <- T
            break
          }
          
          f <- term1 + term2 + term3
          f[f==-Inf] <- NA
          
          del_f <- f-mean(f,na.rm=T) #normalize
          
          target_density <- (exp(del_f)/sum(exp(del_f),na.rm=T))
          target_density[is.na(target_density)] <- 0
          target_F <- cumsum(target_density)
          
          # take inverse draw from CDF
          U <- runif(1)
          if(is.na(min(domain[target_F>U]))|min(domain[target_F>U])==Inf)
            stop(paste0('Y_',area_id,'_',within_area_id,' is about to be set to NA at iteration ', i))
          Y[within_area_id] <- min(domain[target_F>U])
          rho_z_samples[,within_area_id] <- rho_k*rhyper(2e3,min(domain[target_F>U]),N_k-min(domain[target_F>U]),n_k)
          
          # if we got to the end of the loop
          if(within_area_id==length(N))
            try <- F
        }
        
      }
      Y_current[[area_id]] <- Y
      Yplus_current[area_id] <- Yplus
      Y_postsamp[[area_id]][i,] <- Y
      tries[i,area_id] <- try_count
    }
      
    }

    # record current state ---
    reg_params_postsamp[i,] <- reg_params_current
    Yplus_postsamp[i,] <- Yplus_current
    data_table[, ':='(Y=unlist(Y_current))]
    Y_sum_obs_current <- unlist(lapply(Y_current,sum))
  
  }
  
  #combine chains
  beta_postsamp <- cbind(0,reg_params_postsamp[,6:(4+n_admin1)])
  b_postsamp <- reg_params_postsamp[,(5+n_admin1):(4+n_admin1+n_admin2)]
  logtau_postsamp <- reg_params_postsamp[,1]
  logitphi_postsamp <- reg_params_postsamp[,2]
  logd_postsamp <- reg_params_postsamp[,3]
  alphaU_postsamp <- reg_params_postsamp[,4]
  alphaR_postsamp <- reg_params_postsamp[,5]
  
  #admin.key <- pop_dat %>% dplyr::select(A1,A2) %>% unique()
  admin.key <- unique.array(pop_dat[,c('A1','A2')])
  eta_postsamp <- matrix(NA,n_iter,2*n_admin2)
  for(i in 1:nrow(admin.key)){
    eta_postsamp[,i] <- alphaU_postsamp + beta_postsamp[,admin.key$A1[i]] + b_postsamp[,admin.key$A2[i]]
    eta_postsamp[,n_admin2 + i] <- alphaR_postsamp + beta_postsamp[,admin.key$A1[i]] + b_postsamp[,admin.key$A2[i]]
  }
  r_postsamp <- exp(eta_postsamp)
  
  out <- data.frame(cbind(logtau_postsamp,logitphi_postsamp,logd_postsamp,alphaU_postsamp,alphaR_postsamp,beta_postsamp,b_postsamp,
                          r_postsamp,Yplus_postsamp,tree_depth,tries))
  colnames(out)[1:5] <- c('logtau','logitphi','logd','alphaU','alphaR')
  colnames(out)[6:(n_admin1+5)] <- paste0('beta',1:n_admin1)
  colnames(out)[(n_admin1 + 6):(n_admin2 + n_admin1+5)] <- paste0('b',1:n_admin2)
  colnames(out)[(n_admin2 + n_admin1+6):(2*n_admin2 + n_admin1+5)] <- paste0('rU',1:n_admin2)
  colnames(out)[(2*n_admin2 + n_admin1+6):(3*n_admin2 + n_admin1+5)] <- paste0('rR',1:n_admin2)
  colnames(out)[(3*n_admin2 + n_admin1+6):(3*n_admin2 + 2*n_admin1+5)] <- paste0('Yplus',1:n_admin1)
  
  return(out)
}



