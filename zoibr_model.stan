// generated with brms 2.16.3
functions {
  
/* zero-one-inflated beta log-PDF of a single response 
   * Args: 
   *   y: response value 
   *   mu: mean parameter of the beta part
   *   phi: precision parameter of the beta part
   *   zipp: zero-inflation probability parameter
   *   coi: conditional one-inflation probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zoib2_lpdf(real y, real mu, real phi,
                                    real zipp, real coi) {
     row_vector[2] shape = [mu * phi, (1 - mu) * phi]; 
     if (y == 0) { 
       return bernoulli_lpmf(0 | zipp); 
     } else if (y == 1) {
       return bernoulli_lpmf(1 | zipp) + bernoulli_lpmf(1 | coi);
     } else { 
       return bernoulli_lpmf(1 | zipp) + bernoulli_lpmf(0 | coi) + beta_lpdf(y | shape[1], shape[2]);
     } 
   }

}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_zipp;  // number of population-level effects
  matrix[N, K_zipp] X_zipp;  // population-level design matrix
  int<lower=1> K_coi;  // number of population-level effects
  matrix[N, K_coi] X_coi;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K] b;  // population-level effects
  real Intercept_phi;  // temporary intercept for centered predictors
  vector[K_zipp] b_zipp;  // population-level effects
  vector[K_coi] b_coi;  // population-level effects
}
transformed parameters {
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = X * b;
    // initialize linear predictor term
    vector[N] phi = Intercept_phi + rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] zipp = X_zipp * b_zipp;
    // initialize linear predictor term
    vector[N] coi = X_coi * b_coi;
    for (n in 1:N) {
      // apply the inverse link function
      phi[n] = exp(phi[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      zipp[n] = inv_logit(zipp[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      coi[n] = inv_logit(coi[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      mu[n] = inv_logit(mu[n]);
    }
    for (n in 1:N) {
      target += zoib2_lpdf(Y[n] | mu[n], phi[n], zipp[n], coi[n]);
    }
  }
  // priors including constants
  target += student_t_lpdf(Intercept_phi | 3, 0, 2.5);
}
generated quantities {
  // actual population-level intercept
  real b_phi_Intercept = Intercept_phi;
}
