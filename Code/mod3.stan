//
// This Stan program fits cluster level negative binomial model with U/R intercept and BYM2 spatial effect
//

functions {
  
  real pc_tau_lpdf(real tau){
    real lpdf;
    real theta;
    theta = -log(0.01)/1;
    lpdf = log(theta) - 1.5*log(tau) - theta*tau^(-0.5);
    return lpdf;
  }
  
}

data {
  //int<lower=0> lenA1;
  int<lower=0> lenA2;//number of admin2 areas
  int<lower=0> lenC;//number of clusters
 // array[lenC] int admin1_id;
  array[lenC] int urban_id; //1=urban, 0=rural
  array[lenC] int admin2_id;//which admin2 area does each cluster belong to
  array[lenC] int Y;
  array[lenC] int N;
  matrix[lenA2,lenA2] Q_scaled_inv;
}

// The parameters accepted by the model. 
parameters {
  vector[2] alpha;
  vector[lenA2] b;
  real<lower=0> tau;
  real<lower=0, upper=1> phi;
  real<lower=0> d;
}

transformed parameters {
  vector[lenC] Nr;
  matrix[lenA2,lenA2] Sigma_b;
  
  //variance for BYM2 transformation
 Sigma_b= (diag_matrix(rep_vector(1-phi,lenA2)) + phi*Q_scaled_inv)*(1/tau);

 for(i in 1:lenC){
   Nr[i] = N[i]*exp(alpha[2-urban_id[i]] + b[admin2_id[i]]);
 }

}

// The model to be estimated.
model {
  target += pc_tau_lpdf(tau);
  phi ~ beta(0.5,0.5);
  d ~ exponential(1);
  alpha ~ multi_normal(rep_vector(-3.5,2),diag_matrix(rep_vector(9,2)));
  b ~ multi_normal(rep_vector(0,lenA2),Sigma_b);
  for(i in 1:lenC){
    Y[i] ~ neg_binomial(Nr[i]/d,1/d);
  }
  
}

