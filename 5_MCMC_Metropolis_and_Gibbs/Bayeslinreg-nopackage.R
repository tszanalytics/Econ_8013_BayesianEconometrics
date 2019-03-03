#  Linear regression code without using a package
#
rm(list=ls())

# DGP
set.seed(13579)
n      = 1000
beta   = c(1,0.5)
k = length(beta)
x      = cbind(1,rnorm(n))
s = 1.0
y     = x%*%beta + sqrt(s)*rnorm(n)
error = y-x%*%beta

# pdf(file="dgp-linreg.pdf",width=10,height=8)
plot(x[,2],y,xlab="x",ylab="y")
# dev.off()


# prior
b0        = rep(0,k)
V0        = diag(100,k)
iV0       = solve(V0)
iV0b0     = iV0%*%b0
b0iV0b0   = t(b0)%*%iV0b0
d0 = 0.01333
nu0 = 5


# Initial values
ols   = lm(y~x-1)
b     = ols$coef
e     = y-x%*%b
s2ols = t(e)%*%e / (n-k)
s2 = s2ols
nu = nu0 + n

# MCMC 
sd    = 0.1
burn  = 1000
M     = 5000
niter = burn+M
bs    = matrix(0,niter,k)
taus    = rep(0,niter)

for (iter in 1:niter){
    print(iter)

# updating formulas, see Greenberg (2013), p.48 
    iV1 = iV0+t(x)%*%x
    V1 = solve(iV1)
    b1 = V1%*%(iV0b0+t(x)%*%y)
    b  = b1 + t(chol(as.numeric(s2)*V1))%*%rnorm(k)  + rnorm(k)
    d1 = d0 + t(y)%*%y + b0iV0b0 - t(b1)%*%iV1%*%b1
    tau = rgamma(1,shape = nu/2,rate = d1/2)
    s2 = 1/tau
    bs[iter,]   = t(b)
    taus[iter]   = tau
}
bs = bs[seq(burn+1,niter),]
taus = taus[seq(burn+1,niter)]


# pdf(file="hetero-error2.pdf",width=10,height=8)

# regression coefficients
par(mfrow=c(2*k,3))
for (i in 1:k){
  ts.plot(bs[,i],ylab="",xlab="iteration",main=paste("beta",i,sep=""))
  abline(h=beta[i],col=2)
  acf(bs[,i],ylab="",main="")
  hist(bs[,i],prob=T,breaks=100,xlab="",ylab="",main="")
  abline(v=beta[i],col=2)
}

# precision
  ts.plot(taus,ylab="",xlab="iteration",main=paste("tau",i,sep=""))
  abline(h=1/(s^2),col=2)
  acf(taus,ylab="",main="")
  hist(taus,prob=T,breaks=100,xlab="",ylab="",main="")
  abline(v=1/(s^2),col=2)

# variance
  s2s = 1/taus
  ts.plot(s2s,ylab="",xlab="iteration",main=paste("sigma sq.",i,sep=""))
  abline(h=s^2,col=2)
  acf(s2s,ylab="",main="")
  hist(s2s,prob=T,breaks=100,xlab="",ylab="",main="")
  abline(v=s^2,col=2)

# dev.off()

# means and SDs
mean(bs[,1]); sd(bs[,1])
mean(bs[,2]); sd(bs[,2])
mean(s2s); sd(s2s)



