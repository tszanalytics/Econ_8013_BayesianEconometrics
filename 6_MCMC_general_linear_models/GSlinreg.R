# GSlinreg.R
#
# MCMC Gibbs for simple linear regression example
# Simulating the data
set.seed(1244)
n=100; alpha=0; beta=2; sig2=0.5; true=c(alpha,beta,sig2)
x=rnorm(n)
y=rnorm(n,alpha+beta*x,sqrt(sig2))

# Prior hyperparameters
alpha0=0; tau2a=10; beta0=0; tau2b=10; nu0=3; s02=1; nu0s02=nu0*s02

# Setting up starting values
alpha=0; beta=0; sig2=1

# Gibbs sampler
M = 1000
draws = matrix(0,M,3)
for (i in 1:M){
var = 1/(1/tau2a+n/sig2)
mean = var*(sum(y-beta*x)/sig2+alpha0/tau2a)
alpha = rnorm(1,mean,sqrt(var))
var = 1/(1/tau2b+sum(x^2)/sig2)
mean = var*(sum((y-alpha)*x)/sig2+beta0/tau2b)
beta = rnorm(1,mean,sqrt(var))
sig2 = 1/rgamma(1,(nu0+n)/2,(nu0s02+sum((y-alpha-beta*x)^2)/2))
draws[i,] = c(alpha,beta,sig2)
}
# Markov chains + marginal posterior
names = c("alpha","beta","sig2")
ind = 101:M
par(mfrow=c(3,3))
for (i in 1:3){
ts.plot(draws[,i],xlab="iterations",ylab="",main=names[i])
abline(v=ind[1],col=4)
abline(h=true[i],col=2,lwd=2)
acf(draws[ind,i],main="")
hist(draws[ind,i],prob=T,main="",xlab="")
abline(v=true[i],col=2,lwd=2)
}
