//
// This Stan program defines an MCMC algorithm for estimating P(r|Y,Y+) 
// where Y~Pois, Y+~Multinom and log(rj)= a + bj, b~N(0,1/tau)
//

data {
  int<lower=0> lenA;
  int<lower=0> Yplus;
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
  Sigma_b= diag_matrix(rep_vector(1/tau,lenA));
  
}

// The model to be estimated.
model {
  tau ~ gamma(0.01,0.01);
  a ~ normal(-6,10);
  b ~ multi_normal(rep_vector(0,lenA),Sigma_b);
  Yplus ~ poisson(Nxrsum);
  Y ~ multinomial(rsumw);
}

