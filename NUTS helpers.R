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

#is this even right?? I think it should be dot product, not cross
check_NUTS = function(s, theta_plus, theta_minus, r_plus, r_minus){
  if(is.na(s)) return(0)
  condition1 <- crossprod(theta_plus - theta_minus, r_minus) >= 0
  condition2 <- crossprod(theta_plus - theta_minus, r_plus) >= 0
  s && condition1 && condition2
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

