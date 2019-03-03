rm(list=ls())

## using RStan example

## Need following 4 lines before use:
Sys.setenv(USE_CXX14 = 1)
library("rstan") 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


### generate data
set.seed(42)
N = 200
k = 5
sigma = 1
X = cbind(rep(1,N),matrix(rnorm(N*k), N, k))
betas = rnorm(k+1)

# normal errors
y = X%*%betas + sigma*rnorm(N)
y = as.vector(y)
# Student-t errors
v = 5
yt = X%*%betas + sigma*rt(N,df=v)
yt = as.vector(yt)

# examine data
par(mfrow=c(2,2))
plot(X[,2],y,lwd=2,col=2)
plot(X[,3],y,lwd=2,col=2)
plot(X[,4],y,lwd=2,col=2)
plot(X[,5],y,lwd=2,col=2)

r = lm(y~X[,2:(k+1)])
summary(r)

hist(y,breaks=20,col=4)


# Cauchy prior for sigma?
s_draws = rcauchy(10000,scale=2.5)
plot(s_draws)
hist(s_draws,breaks=100)
plot(density(s_draws))

s2_draws = s_draws^2
plot(s2_draws)
hist(s2_draws,breaks=100)
plot(density(s2_draws))

tau_gamdraws = rgamma(10000,shape = 0.001, scale=0.001)
plot(tau_gamdraws)
hist(tau_gamdraws,breaks=100)
plot(density(tau_gamdraws))

s2_gamdraws = 1/tau_gamdraws
plot(s2_gamdraws)
hist(s2_gamdraws,breaks=100)
plot(density(s2_gamdraws))


## use Stan to estimate linear regression
# Specify the data list that we will pass to Stan. This gives Stan everything declared in the data{} block. 
data_reg = list(X = X, N = N, y = y, k = (k+1))   # k covariates + intercept

# Call Stan. You'll need to give it a file (.stan file), 
# or a fitted Stan object (fit)
# You should also pass Stan a data list, number of cores to estimate on, 
# the number of Markov chains to run (4 by default)
# and number of iterations (2000 by default). 
# We use multiple chains to make sure that the posterior distribution that we converge on 
# is stable, and not affected by starting values. 

# The first time you run the models, they will take some time to compile before sampling. 
# On subsequent runs, it will only re-compile if you change the model code. 

model_fit <- stan(file='linear_reg.stan', data = data_reg, cores = 4, chains = 4, iter = 2000)

## examine results
print(model_fit)
plot(model_fit)
#pairs(model_fit, pars = c("beta", "sigma", "lp__"))
pairs(model_fit, pars = c("beta[2]","beta[3]","beta[4]", "sigma"))

chains <- extract(model_fit, permuted = TRUE) # return a list of arrays 

# posterior for mu
mu <- chains$beta[,2] 
plot(density(mu))
abline(v=mean(mu))
mean(mu)
sd(mu)

# posterior for sigma
s <- chains$sigma 
hist(s,col=3,breaks=100,freq=F)
lines(density(s),col=4)
abline(v=mean(s))
mean(s)
sd(s)

 
#pairs(model_fit, pars = c("beta"))

b1 = chains$beta[,1] 
b2 = chains$beta[,2]

mycol <- rgb(0, 0, 255, max = 255, alpha = 25, names = "blue50")
plot(b1,b2,pch=16,col=mycol)
abline(v=mean(b1),col=2)
abline(h=mean(b2),col=2)


# heteroskedastic (t-errors) data
data_reg = list(X = X, N = N, y = yt, k = (k+1))   # k covariates + intercept


### Student-t errors model
model_fit2 <- stan(file='linear_reg_t_errors.stan', data = data_reg, cores = 4, chains = 4, iter = 2000)

## examine results
print(model_fit2)
plot(model_fit2)
pairs(model_fit2, pars = c("beta[2]","beta[3]","beta[4]", "nu"))
# pairs(model_fit, pars = c("beta", "sigma", "lp__"))

chains2 <- extract(model_fit2, permuted = TRUE) # return a list of arrays 

# posterior for mu
nu_draws <- chains2$nu 
par(mfrow=c(1,1))
hist(nu_draws,breaks=50,freq=F)
lines(density(nu_draws),col=4,lwd=2)
abline(v=mean(nu_draws),lwd=2)
mean(nu_draws)
sd(nu_draws)
quantile(nu_draws,probs=(c(0.005,0.995)))
quantile(nu_draws,probs=(c(0.05,0.95)))

# test vs. Normality (nu = 30 or nu >= 30)
# nu = 30
source(file="postoddsmc.R")
postoddsmc(nu_draws,null=30.0)

# prob >= 30
probgt30 = sum(nu_draws >= 30.0)/length(nu_draws)
probgt30  # "Bayesian p-value" = area in tail

# odds against >=30
problt30 = 1.0 - probgt30
oddsgt30 = problt30/probgt30
oddsgt30
