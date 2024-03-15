#' No-U-Turn sampler --- https://github.com/kasparmartens/NUTS/blob/master/R
#'
#' @param theta Initial value for the parameters
#' @param f log-likelihood function (up to a constant)
#' @param grad_f the gradient of the log-likelihood function
#' @param n_iter Number of MCMC iterations
#' @param M_diag Diagonal elements of the mass matrix in HMC. Defaults to ones.
#' @param M_adapt Parameter M_adapt in algorithm 6 in the NUTS paper
#' @param delta Target acceptance ratio, defaults to 0.5
#' @param max_treedepth Maximum depth of the binary trees constructed by NUTS
#' @param eps Starting guess for epsilon
#' @return Matrix with the trace of sampled parameters. Each mcmc iteration in rows and parameters in columns.
#' @export


NUTS <- function(theta, f, grad_f, n_iter, M = NULL, M_adapt = 50, delta = 0.5, max_treedepth = 10, eps = 1, verbose = TRUE){
  theta_trace <- matrix(0, n_iter, length(theta))
  par_list <- list(M_adapt = M_adapt)
  for(iter in 1:n_iter){
    nuts <- NUTS_one_step(theta, iter, cons_params, par_list, delta = delta, max_treedepth = max_treedepth, eps = eps, verbose = verbose)
    theta <- nuts$theta
    par_list <- nuts$pars
    theta_trace[iter, ] <- theta
    print(iter)
  }
  theta_trace
}

NUTS_one_step <- function(theta, iter, f, grad_f, cons_params, par_list, max_treedepth = 10, verbose = TRUE,
                          DA=T, 
                          kappa = 0.75, t0 = 10, # stabilizes initial iterations
                          gamma = 0.05, #controls shrinkage towards mu
                          delta = 0.5 ){ # target acceptance prob
  
  if(DA){
    eps_adapt <- par_list$eps_adapt
  }
  
  M <- par_list$M
  M_inv <- solve(M)
  M_diag <- diag(M)
  
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
  
  r0 <- as.vector(rmvnorm(1, rep(0,length(theta)), M)) #optimize later
  u <- runif(1, 0, exp(f(params = theta,cons_params = cons_params) - 0.5 * t(r0)%*%M_inv%*%r0))
  
  if(is.nan(u)){
    warning("NUTS: sampled slice u is NaN")
    u <- runif(1, 0, 1e5)
  }
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
      temp <- build_tree(theta_minus, r_minus, u, direction, j, eps, theta, r0, f, grad_f, cons_params, M, M_inv)
      theta_minus <- temp$theta_minus
      r_minus <- temp$r_minus
    } else{
      temp <- build_tree(theta_plus, r_plus, u, direction, j, eps, theta, r0, f, grad_f, cons_params, M, M_inv)
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

