# example of single equation regression iid MCMC draws
# see ch. 2 Ross-Allenby-McCulloch (2005)
library(bayesm)

# rng seed
set.seed(66)

# n = obs, R = iterations of MCMC
n=200
R=10000

# Simulate data
X=cbind(rep(1,n),runif(n),runif(n)); beta=c(1,2,3); sigsq=.25
y=X%*%beta+rnorm(n,sd=sqrt(sigsq))

out=runireg(Data=list(y=y,X=X),Mcmc=list(R=R))

cat("Summary of beta/sigma-sq draws",fill=TRUE)
summary(out$betadraw,tvalues=beta)
summary(out$sigmasqdraw,tvalues=sigsq)

if(0){
## plotting examples
plot(out$betadraw)
}

hist(out$betadraw[,2],col="blue")
plot(out$betadraw[,2],type='l',col="blue",lwd=2)
acf(out$betadraw[,2],xlab="Lag",ylab="ACF",main="")
pacf(out$betadraw[,2],xlab="Lag",ylab="ACF",main="")

plot(density(out$betadraw[,2]),col="blue")