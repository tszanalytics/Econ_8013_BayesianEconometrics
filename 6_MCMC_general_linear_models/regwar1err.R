# regwar1err.R
# regression with ar(1) errors
# simulates data and calls acreg to perform MCMC

# Simulated data
set.seed(12357)
n      = 50
sigma2 = 5
p      = 2
beta   = c(2,4)
rho    = 0.95
sigma  = sqrt(sigma2)
x      = matrix(rnorm(p*n),n,p)
Om      = matrix(0,n,n)
for (i in 1:n)
  Om[i,] = rho^abs(((1-i):(n-i)))
Om = Om/(1-rho^2)  
y     = 1 + x%*%beta+sigma*t(chol(Om))%*%rnorm(n)

# MCMC for AR(1) errors model
R <- 3000
mc1 <- acreg(y,x,R) # using cov. Omega construction and imposing |rho|<1
# mc1 <- acregalt(y,x,R) # avoiding Omega constr. and choice of |rho|<1 or not
str(mc1)

bs <- mc1$bs
s2s <- mc1$s2s
rs <- mc1$rs

par(mfrow=c(2,3))
hist(bs[,1],col=4)
abline(v=mean(bs[,1]),col=2)
hist(bs[,2],col=4)
abline(v=mean(bs[,2]),col=2)
hist(bs[,3],col=4)
abline(v=mean(bs[,3]),col=2)
hist(rs,col=4)
abline(v=mean(rs),col=2)
hist(s2s,col=4)
abline(v=mean(s2s),col=2)

mean(bs[,2])
sd(bs[,2])

mean(rs)
sd(rs)

k <- ncol(bs)

par(mfrow=c(k,3))
for (i in 1:k){
  ts.plot(bs[,i],ylab="",xlab="iteration",main=paste("beta",i,sep=""))
  abline(h=beta[i],col=2)
  acf(bs[,i],ylab="",main="")
  hist(bs[,i],prob=T,xlab="",ylab="",main="")
  abline(v=beta[i],col=2)
}

par(mfrow=c(2,3))
ts.plot(s2s,ylab="",xlab="iteration",main=expression(sigma^2))
abline(h=sigma2,col=2)
acf(s2s,ylab="",main="")
hist(s2s,prob=T,xlab="",ylab="",main="")
abline(v=sigma2,col=2)
ts.plot(rs,ylab="",xlab="iteration",main=expression(rho))
abline(h=rho,col=2)
acf(rs,ylab="",main="")
hist(rs,prob=T,xlab="",ylab="",main="")
abline(v=rho,col=2)



