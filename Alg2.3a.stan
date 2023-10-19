//
// This Stan program defines an MCMC algorithm for estimating P(b|Y) where Y~Pois and b~N(0,1/tau)
//

data {
  int<lower=0> lenA;
  int Y[lenA];
  int N[lenA];
}

// The parameters accepted by the model. 
parameters {
  real a;
  vector[lenA] b;
  real<lower=0> tau;
}

transformed parameters {
  real r[lenA];
  matrix[lenA,lenA] Sigma_b;
  
  for(i in 1:lenA){
    r[i] = exp(a + b[i]);
  }
  //variance for BYM2 transformation
  Sigma_b= diag_matrix(rep_vector(1/tau,lenA));
}

// The model to be estimated.
model {
  tau ~ gamma(0.01,0.01); //same as INLA
  a ~ normal(-6,10);
  b ~ multi_normal(rep_vector(0,lenA),Sigma_b);
  for(i in 1:lenA){
    Y[i] ~ poisson(N[i]*r[i]);
  }
}

