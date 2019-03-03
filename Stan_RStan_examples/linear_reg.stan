data {
  // Define the data that must be passed to Stan 

  int N; // number of observations
  int k; // number of covariates including intercept
  matrix[N, k] X; //covariate matrix including intercept
  vector[N] y; //outcome vector
}
parameters {
  // Define the parameters to estimate, as well as any restrictions on the parameter values   
  // (standard deviations can't be negative...)

  vector[k] beta;         // the regression coefficients
  real<lower = 0> sigma;  // the residual standard deviation (restricted to be non-negative)
}
model {
  // This is where we write out the probability model

  // Define the priors  
  beta ~ normal(0, 5);     // using same prior for all betas here (not necessary)
  // sigma ~ cauchy(0, 2.5);  // Example uses a Cauchy for sigma ?!
  sigma ~ gamma(2,0.001);  // suggested gamma prior for sigma
  // The likelihood
  y ~ normal(X*beta, sigma);
}

   

  