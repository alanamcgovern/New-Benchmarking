//
// This Stan program defines an MCMC algorithm for estimating P(r|Y,Y+) where Y~Pois, Y+~Multinom and log(rj)= a + bj, b~BYM2(tau,phi) with PC priors on tau and phi
//

functions {
  real pc_tau_lpdf(real tau, real U, real alpha){
    real lpdf;
    real theta;
    theta = -log(alpha)/U;
    lpdf = log(0.5*theta) - 1.5*log(tau) - theta*tau^(-0.5);
    return lpdf;
  }
  
  real kld_phi(real phi, matrix Q_inv){
    real out;
    out = phi*trace(Q_inv) - phi*rows(Q_inv) - log_determinant(diag_matrix(rep_vector(1-phi,rows(Q_inv))) + phi*Q_inv);
    return out;
  }
  
  real pc_phi_lpdf(real phi, real lambda, matrix Q_inv, vector eigen){
    real out;
    real kld;
    
    kld = kld_phi(phi,Q_inv);
    out = log(lambda) - 0.5*log(8*kld) - lambda*sqrt(2*kld) + log(trace(Q_inv)-rows(Q_inv)-sum((eigen-1)./(1+phi*(eigen-1))));
    return out;
  }
}

data {
  int<lower=0> lenA;
  int<lower=0> Yplus;
  int Y[lenA];
  int N[lenA];
  matrix[lenA,lenA] Q_scaled_inv;
  //eigenvalues of Q_scaled_inv
  vector[lenA] eigen;
  //lambda chosen based on U and alpha for PC prior on phi
  real<lower=0> lambda;
}

// The parameters accepted by the model. 
parameters {
  real a;
  vector[lenA] b;
  real<lower=0> tau;
  real<lower=0,upper=1> phi;
}

transformed parameters {
  real r[lenA];
  matrix[lenA,lenA] Sigma_b;
  
  vector[lenA] Nxr;
  // sum(N_j*r_j)
  real Nxrsum;
  // N_i*r_i/sum(N_j*r_j)
  vector[lenA] rsumw;
  
  for(i in 1:lenA){
    r[i] = exp(a + b[i]);
    Nxr[i] = N[i]*r[i];
  }
  Nxrsum=sum(Nxr);
  //probability for multinomial
  for(i in 1:lenA){
    rsumw[i] = Nxr[i]/Nxrsum;
  }
  
  //variance for BYM2 transformation
  Sigma_b= ((1-phi)*diag_matrix(rep_vector(1,lenA)) + phi*Q_scaled_inv)*(1/tau);
  
}

// The model to be estimated.
model {
  target += pc_tau_lpdf(tau | 1, 0.01);
  target += pc_phi_lpdf(phi | lambda, Q_scaled_inv, eigen);
  a ~ normal(-6,10);
  b ~ multi_normal(rep_vector(0,lenA),Sigma_b);
  Yplus ~ poisson(Nxrsum);
  Y ~ multinomial(rsumw);
}

