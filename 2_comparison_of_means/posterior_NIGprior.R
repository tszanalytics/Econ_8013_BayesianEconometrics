# Normal-Gamma prior for x NOrmmal with unknown mean and variance

##p(mu|x) = Student-t(mu1, b1/(lam1*a1),2*a1)


## input: x = data (one variable)
## output: mu_draws = M draws from the t-posterior for mu = mean of x

posterior_NIGprior = function(x,mu0=0.0, lam0=0.00001,a0 = 0.0001, b0=0.0001,M=100000) {
  ## prior parameter values here:
  #mu0 = 0.0       # prior mean
  #lam0 = 0.00001  # prior precision = 1/var
  #a0 = 0.0001     # prior "degrees of freedom"
  #b0 = 0.0001     # prior scale parameter for Gamma
  #M = 100000      # number of draws from posterior

  # compute posterior parameter values
  n = length(x)
  mx = mean(x)
  s2 = var(x)
  mu1 = (lam0*mu0 + n*mx)/(lam0 + n)
  lam1 = lam0 + n
  a1 = a0 + n/2
  b1 = b0 + 0.5*(n*s2 + (lam0*n*(mx - mu0)^2)/(lam0 + n))
  tsd = sqrt(b1/(lam1*a1))
  tdf = 2*a1
  mu_draws = mu1 + tsd*rt(M, df=tdf)
return(mu_draws)
}



