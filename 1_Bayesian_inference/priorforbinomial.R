#
# Priors and posteriors with binomial likelihood
# 
# NB no need to use loops, just define the sequence (as in normdensity.R)

# R = number of evaluations of distns
R = 1001

# uniform prior = 1 for all 0 < theta < 1
uniform = rep(1,R)

# univariate numerical integration to determine normalizing constant
# specific case works (setting parameters in function)
# int1 <- function(z) (exp(-((z-0.5)^2)/(2*100)))*(z^5)*(1-z)^5
# cth = integrate(int1, lower = 0, upper = 1.0)
# c = cth$value
# c

# prior parameters for truncated normal (and normalizing constant evaluation)
m0 = 0.9
s20 = 1

int1 <- function(z) (exp(-((z-m0)^2)/(2*s20)))
cth = integrate(int1, lower = 0, upper = 1.0)
c = cth$value
c


# evaluate and plot normal prior
th = rep(0,R)
normpr = rep(0,R)
for (i in 1:R)
{
th[i] = i*0.001 - 0.001
normpr[i] = (1/c)*(exp(-((th[i]-m0)^2)/(2*s20)))

}
plot(th,normpr,type='l',lwd=2)
title("Truncated normal prior centered at 0.9, var = 1.0")

# evaluate and plot beta-binomial specification
th = rep(0,R)
normbin = rep(0,R)
betaprior = rep(0,R)
betapost = rep(0,R)
unifpost = rep(0,R)

# parameters for beta prior and likelihood
a = 2
b = 1
n = 5
y = 1


# univariate numerical integration to determine normalizing constant for normbin
int1 <- function(z) (exp(-((z-m0)^2)/(2*s20)))*(z^y)*(1-z)^(n-y)
cth = integrate(int1, lower = 0, upper = 1.0)
c = cth$value
c

for (i in 1:R)
{
th[i] = i*0.001 - 0.001
normbin[i] = (1/c)*(exp(-((th[i]-m0)^2)/(2*s20)))*(th[i]^y)*(1-th[i])^(n-y)
betaprior[i] = (gamma(a+b)/(gamma(a)*gamma(b)))*(th[i]^(a-1))*(1-th[i])^(b-1)
betapost[i] = (gamma(n+a+b)/(gamma(y+a)*gamma(n-y+b)))*(th[i]^(y+a-1))*(1-th[i])^(n-y+b-1)
unifpost[i] = (gamma(n+2)/(gamma(y+1)*gamma(n-y+1)))*(th[i]^y)*(1-th[i])^(n-y)
}
yy = cbind(normpr,normbin,betaprior,betapost,unifpost)
matplot(th,yy,type='l', col= 1:5,lwd=2)
legend("topright",legend=c("Normal prior","Normal posterior","Beta prior", "Beta posterior", "Uniform posterior"),col=c(1:5),lty=c(1:5),lwd=2,bty="n",cex=1.1)
title("Priors and posteriors with 5 observations")


