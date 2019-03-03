###########################################################################################
#
#  General linear model with Student t errors
#
###########################################################################################
#
# Author : Hedibert Freitas Lopes
#          Booth School of Business
#          University of Chicago
#          5807 South Woodlawn Avenue
#          Chicago, Illinois, 60637
#          Email : hlopes@ChicagoBooth.edu
#
###########################################################################################
rm(list=ls())

lfull=function(nu){(n*nu/2)*log(nu/2)-n*lgamma(nu/2)-eta*nu}

# Simulated data
# --------------
set.seed(12434)
n      = 100
nu     = 10
sigma2 = 1.0
p      = 3
beta   = c(0.5,1.0,2.0)
sigma  = sqrt(sigma2)
x      = matrix(rnorm(p*n),n,p)
omega  = 1/rgamma(n,nu/2,nu/2)
O      = sigma2*diag(omega)
error  = rnorm(n,0,sqrt(sigma2*omega))
y      = x%*%beta+error

# pdf(file="t-error1.pdf",width=10,height=8)
par(mfrow=c(1,2))
ts.plot(error,xlab="Observation",ylab=expression(epsilon[i]))
xxx = seq(min(error),max(error),length=40)
hist(error,prob=TRUE,breaks=xxx,xlab="",ylab="",main=expression(epsilon[i]))
xxx = seq(min(error),max(error),length=1000)
lines(xxx,dt(xxx/sigma,nu)/sqrt(sigma),col=4)
# dev.off()

# prior
# -----
b0      = rep(0,p)
V0      = diag(10,p)
iV0     = solve(V0)
iV0b0   = iV0%*%b0
nu0     = 5
s02     = 1.0
nu0.lam = 25
nu0n    = nu0+n
nu0s02  = nu0*s02

# Initial values
# --------------
ols   = lm(y~x-1)
b     = ols$coef
s2    = mean(ols$res^2)
iO    = diag(1/omega)
xiOx  = t(x)%*%iO%*%x
xiOy  = t(x)%*%iO%*%y
nu.l  = 10

# MCMC scheme
# -----------
set.seed(23164)
dd    = 0.25
burn  = 10000
M     = 1000
LAG   = 100
niter = burn+LAG*M
s2s   = rep(0,niter)
nuls  = rep(0,niter)
bs    = matrix(0,niter,3)
lams  = matrix(0,niter,n)
for (i in 1:niter){
  print(i)
  par2  = nu0s02 + t(y-x%*%b)%*%iO%*%(y-x%*%b)
  s2    = 1/rgamma(1,nu0n/2,par2/2)
  V1    = solve(iV0+xiOx/s2)
  b1    = V1%*%(iV0b0+xiOy/s2)
  b     = b1+t(chol(V1))%*%rnorm(p)
  lam   = rgamma(n,(nu.l+1)/2,(nu.l+(y-x%*%b)^2/s2)/2)
  eta   = 1/nu0.lam + 0.5*sum(log(1/lam)+lam)
  nu.ll = rnorm(1,nu.l,dd)
  accep = min(0,lfull(nu.ll)-lfull(nu.l))
  if (log(runif(1))<accep){
    nu.l = nu.ll
  }
  iO       = diag(lam)
  xOx      = t(x)%*%iO%*%x
  xiOy     = t(x)%*%iO%*%y
  s2s[i]   = s2
  bs[i,]   = t(b)
  lams[i,] = lam
  nuls[i]  = nu.l  
}
bs     = bs[seq(burn+1,niter,by=LAG),]
lams   = lams[seq(burn+1,niter,by=LAG),]
s2s    = s2s[seq(burn+1,niter,by=LAG)]
nuls   = nuls[seq(burn+1,niter,by=LAG)]
omegas = 1/lams

# pdf(file="t-error2.pdf",width=10,height=8)
par(mfrow=c(3,3))
for (i in 1:p){
  ts.plot(bs[,i],ylab="",xlab="iteration",main=paste("beta",i,sep=""))
  abline(h=beta[i],col=2)
  acf(bs[,i],ylab="",main="")
  hist(bs[,i],prob=T,xlab="",ylab="",main="")
  abline(v=beta[i],col=2)
}
# dev.off()


# pdf(file="t-error3.pdf",width=10,height=8)
par(mfrow=c(2,3))
ts.plot(s2s,ylab="",xlab="iteration",main=expression(sigma^2))
abline(h=sigma2,col=2)
acf(s2s,ylab="",main="")
hist(s2s,prob=T,xlab="",ylab="",main="")
abline(v=sigma2,col=2)
ts.plot(nuls,ylab="",xlab="iteration",main=expression(nu))
abline(h=nu,col=2)
acf(nuls,ylab="",main="")
hist(nuls,prob=T,xlab="",ylab="",main="")
abline(v=nu,col=2)
# dev.off()



