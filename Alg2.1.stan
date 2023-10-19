//
// This Stan program defines an MCMC algorithm for estimating P(b|Y) where Y~Pois and b's have normal prior
//

data {
  int<lower=0> lenA;
  int Y[lenA];
  int N[lenA];
}

// The parameters accepted by the model. 
parameters {
  real b[lenA];
}

transformed parameters {
  real r[lenA];
  for(i in 1:lenA){
   r[i] = exp(b[i]);
  }
}


// The model to be estimated.
model {
  for(i in 1:lenA){
    b[i] ~ normal(-2,2);
    Y[i] ~ poisson(N[i]*r[i]);
  }
}


