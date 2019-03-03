#
# heterosked1.R
# heteroskedasticity
#
set.seed(23457)
n <- 100
beta <- c(10,1,1)
alpha <- c(-3,0.3)

alpha[2]
beta[3]

x2 <- runif(n,min=10,max=20)
x3 <- runif(n,min=5,max=10)
cbind(x2,x3)

# heterskedastic variance as a function of x2
sig2 <- exp(alpha[1] + alpha[2]*x2)
sig = sqrt(sig2)

y <- beta[1] + beta[3]*x2 + beta[3]*x3 + rnorm(n,mean=0,sd=sig)

# Suppose we wish to model the heterosked. as exp(a1 + a2*x2)
# (i.e. somehow we know the right variable and the exact form!)


# Generalized least squares
r1 <- lm(y~x2+x3)
ehat <- r1$residuals

q <- log(ehat^2)
z <- cbind(rep(1,n),x2)
ahat <- solve(t(z)%*%z)%*%t(z)%*%q
lm(q~x2)

sh2 <- exp(z%*%ahat)

cbind(sig2,sh2)

V <- diag(n)
diag(V) <- sh2
Vi = solve(V)


X <- cbind(rep(1,n),x2,x3)
k = ncol(X)
bhols <- solve(t(X)%*%X)%*%t(X)%*%y
bhgls <- solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi%*%y
sh2gls <- t(y - X%*%bhgls)%*%Vi%*%(y - X%*%bhgls)/(n-k)
sh2g = as.numeric(sh2gls)
covbgls <- sh2g*solve(t(X)%*%Vi%*%X)
sebgls <-sqrt(diag(covbgls))


# Since we generated the data ourselves, we can compare with the 
# covariance matrix based on the actual values of the parameters

# "Infeasible GLS"  Vk = known cov. matrix

Vk <- diag(n)
diag(Vk) <- sig2
Vki = solve(Vk)
bkgls <- solve(t(X)%*%Vki%*%X)%*%t(X)%*%Vki%*%y
covbgls <- solve(t(X)%*%Vki%*%X)
sebgls <-sqrt(diag(covbgls))

# BP test for heterosked.
sh2ols <- t(y - X%*%bhols)%*%(y - X%*%bhols)/n
checksh2 <- t(ehat)%*%ehat/n

bpy <- ehat^2/sh2ols
rbp <-lm(bpy~x2)
ebp <- rbp$residuals
mbpy = mean(bpy)
tss = sum((bpy-mbpy)^2)
rss = t(ebp)%*%ebp
bplm <- 0.5*(tss - rss)
lmpval <- 1 - pchisq(bplm, df=1)


# nonlinear ML estimation - not correct (for ahml most likely), fix!!?

R = 100
bhml = rep(0,R*k)
dim(bhml) <- c(R,k)
p = length(ahat)
ahml = rep(0,R*p)
dim(ahml) <- c(R,p)
bhml[1,] <- t(bhgls)
ahml[1,] <- t(ahat)

for (i in 2:R) {
s2 <- exp(z%*%ahml[i-1,])
diag(V) <- s2
Vi = solve(V)

bhi <- bhml[i-1,] + solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi%*%(y - X%*%bhml[i-1,])
bhml[i,] <- bhi

ahpart = exp(-z%*%ahml[i-1,])*(y - X%*%bhi)*(y - X%*%bhi) - 1
Zy <- c(z[,1]%*%ahpart, z[,2]%*%ahpart)
ZZi <- solve(t(z)%*%z)
ahi <- ahml[i-1,] + ZZi%*%Zy
ahml[i,] <- ahi
}


# Bayesian inference
p = length(bhols)
q = length(ahat)

# prior
b0        = rep(0,p)
V0        = diag(100,p)
a0        = rep(0,p)
C0        = diag(100,p)
iV0       = solve(V0)
iV0b0     = iV0%*%b0
iC0       = solve(C0)
iC0a0     = iC0%*%a0

# Initial values
x = X
ols   = lm(y~x-1)
b     = ols$coef
e     = y-x%*%b
a <- ahat

# MCMC
sd    = 0.1
burn  = 1000
R     = 200000
niter = burn+R
bs    = matrix(0,niter,p)
as    = matrix(0,niter,q)
for (iter in 1:niter){
A  = exp(z%*%a/2)
  y1 = y/A
  x1 = x/matrix(A,n,p)  
  V1 = solve(iV0+t(x1)%*%x1)
  b1 = V1%*%(iV0b0+t(x1)%*%y1)
  b  = b1+t(chol(V1))%*%rnorm(p)
  e  = y-x%*%b
  a1 = rnorm(q,a,sd)
  la = sum(dnorm(e,0,exp(z%*%a1/2),log=TRUE))-sum(dnorm(e,0,exp(z%*%a/2),log=TRUE))
  if (log(runif(1))<la){
    a = a1
  }
  bs[iter,]   = t(b)
  as[iter,]   = t(a)
}
bs = bs[(burn+1):R,]
as = as[(burn+1):R,]


bmean <- apply(bs, 2, mean)
bsd <- apply(bs, 2, sd)
amean <- apply(as, 2, mean)
asd <- apply(as, 2, sd)

cbind(bmean,bsd)
cbind(amean,asd)

par(mfrow=c(p,3))
for (i in 1:p){
  ts.plot(bs[,i],ylab="",xlab="iteration",main=paste("beta",i,sep=""))
  abline(h=beta[i],col=2)
  acf(bs[,i],ylab="",main="")
  hist(bs[,i],prob=T,xlab="",ylab="",main="")
  abline(v=beta[i],col=2)
}


for (i in 1:q){
  ts.plot(as[,i],ylab="",xlab="iteration",main=paste("alpha",i,sep=""))
  abline(h=alpha[i],col=2)
  acf(as[,i],ylab="",main="")
  hist(as[,i],prob=T,xlab="",ylab="",main="")
  abline(v=alpha[i],col=2)
}





