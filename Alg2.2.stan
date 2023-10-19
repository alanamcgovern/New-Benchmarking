//
// This Stan program defines an MCMC algorithm for estimating P(b|Y) where Y~Pois and b's have normal prior
//

data {
  int<lower=0> lenA;
  int<lower=0> Yplus;
  int<lower=0> Y[lenA];
  vector[lenA] N;
}

// The parameters accepted by the model. 
parameters {
  real b[lenA];
}

transformed parameters {
  vector[lenA] r;
  // N_i*r_i
  vector[lenA] Nxr;
  // sum(N_j*r_j)
  real Nxrsum;
  // N_i*r_i/sum(N_j*r_j)
  vector[lenA] rsumw;
  
  for(i in 1:lenA){
    r[i] = exp(b[i]);
    Nxr[i] = N[i]*r[i];
  }
  Nxrsum=sum(Nxr);
  for(i in 1:lenA){
    rsumw[i] = Nxr[i]/Nxrsum;
  }
}


// The model to be estimated.
model {
  for(i in 1:lenA){
    b[i] ~ normal(-2,10);
  }
  Yplus ~ poisson(Nxrsum);
  Y ~ multinomial(rsumw);
}
