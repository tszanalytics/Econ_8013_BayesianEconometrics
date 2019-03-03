# binomialposteriors.R
#
# Evaluates and plots posteriors for binomial experiment with 
# three different priors: normal, beta and uniform

# set vectors to store values
R = 1001
th = rep(0,R)
normbin = rep(0,R)
betaprior = rep(0,R)
betapost = rep(0,R)
unifpost = rep(0,R)
normprior = rep(0,R)

# parameters for beta prior
a = 2
b = 2

# n = number of obs, y = number of 'successes' in n trials
n = 20
y = 4

# prior parameters for normal prior
m0 = 0.9
s20 = 100

# univariate numerical integration to determine normalizing constant for normbin spec
int1 <- function(z) (exp(-((z-m0)^2)/(2*s20)))*(z^y)*(1-z)^(n-y)
cth = integrate(int1, lower = 0, upper = 1.0)
# c is the normalizing constant for use with normbin
c = cth$value
c

# Evaluate each posterior from 0 to 1, and the beta and normal priors
# normbin = posterior with normal prior
# betapost = posterior with beta prior
# unifpost = posterior with uniform prior
# betaprior = beta prior
# normprior = normal prior
for (i in 1:R)
{
th[i] = i*0.001 - 0.001
normbin[i] = (1/c)*(exp(-((th[i]-m0)^2)/(2*s20)))*(th[i]^y)*(1-th[i])^(n-y)
betaprior[i] = (gamma(a+b)/(gamma(a)*gamma(b)))*(th[i]^(a-1))*(1-th[i])^(b-1)
betapost[i] = (gamma(n+a+b)/(gamma(y+a)*gamma(n-y+b)))*(th[i]^(y+a-1))*(1-th[i])^(n-y+b-1)
unifpost[i] = (gamma(n+2)/(gamma(y+1)*gamma(n-y+1)))*(th[i]^y)*(1-th[i])^(n-y)
normprior[i] = (1/sqrt(2*3.142*s20))*exp(-((th[i]-m0)^2)/(2*s20))

}

# plot of everything on one graph
yy = cbind(normbin,betaprior,betapost,unifpost,normprior)
matplot(th,yy,type='l', col= 1:5,lwd=2)
legend("topright",legend=c("Norm-post","Beta-prior", "Beta-post", "Unif-post", "Norm-prior"),col=c(1:5),lty=c(1:5),lwd=2,bty="n",cex=1.1)
# NB: use lty=c(1,1) in above (matplot and legend) for all solid lines


par(mfrow=c(1,2))

# plot of priors
# uniform prior = 1 for all 0 < theta < 1
unif = rep(1,R)
yy = cbind(unif,betaprior,normprior)
matplot(th,yy,type='l', col= 1:3,lwd=2)
legend("topright",legend=c("Uniform prior","Beta prior", "Normal prior"),col=c(1:3),lty=c(1:3),lwd=2,bty="n",cex=1.1)
# NB: use lty=c(1,1) in above (matplot and legend) for all solid lines


# plot of posteriors
yy = cbind(unifpost,betapost,normbin)
matplot(th,yy,type='l', col= 1:3,lwd=2)
legend("topright",legend=c("Uniform post","Beta post", "Normal post"),col=c(1:3),lty=c(1:3),lwd=2,bty="n",cex=1.1)
# NB: use lty=c(1,1) in above (matplot and legend) for all solid lines




