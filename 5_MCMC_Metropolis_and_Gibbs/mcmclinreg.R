# mcmclinreg.R
# 
# DIY MCMC for the linear regression model example

# rng seed
set.seed(12357)

# n = obs, R = iterations of MCMC
n=200
R=2000

# data from Koop chapter 3 and 4
df = read.table("HPRICE.txt", header = TRUE)
str(df)
summary(df)
# try using this to reproduce the results in chapter 4
hprice <- df$price
lsize <- df$lot_siz
bed <- df$bed
bath <- df$bath
stor <- df$stor
y <- hprice
n <- length(y)
X <- cbind(rep(1,n),lsize,bed,bath,stor)
k = ncol(X)

# Simulate data
# X=cbind(rep(1,n),runif(n),runif(n)); beta=c(1,2,3); sigsq=0.5
# y=X%*%beta+rnorm(n,sd=sqrt(sigsq))
# k = ncol(X)

# Informative prior Bayesian analysis

k = ncol(X)
n = nrow(X)

# Prior hyperparameters
# variance v0 picked as 3*sqrt(v0) giving reasonable "99%" interval

# *** INCORRECT INFORMATIVE PRIOR TO DEMSONSTRATE THE EFFECT
# b0      = c(1,2,3) # correct parameter values
# b0      = c(0,0,0)
b0        = rep(0,k) # default uninformative prior
# Vb0    = c(0.01,0.01,0.01)*diag(1,k)  # very informative prior
# default diffuse prior:
Vb0   = diag(10^10,k)

nu0 = 3
s02 = 1
iVb0 = solve(Vb0)

nu1      = nu0+n
nu0s02    = nu0*s02

# Initial values
ols   = lm(y~X-1)
bols     = ols$coef
s2ols    = mean(ols$res^2)


V1 = solve(iVb0+t(X)%*%X)
b1 = V1%*%(iVb0%*%b0+ t(X)%*%y)
b  = b1+t(chol(V1))%*%rnorm(k)


R <- 2000
burn  = 0
LAG   = 1
niter = burn+LAG*R
s2s   = rep(0,niter)
bs    = matrix(0,niter,k)
for (iter in 1:niter){
  print(iter)
  e  = y-X%*%b
  s2 = 1/rgamma(1,nu1/2,(nu0s02+t(e)%*%e)/2)  
  V1 = solve(iVb0+t(X)%*%X/s2)
  b1 = V1%*%(iVb0%*%b0+t(X)%*%y/s2)
  b  = b1+t(chol(V1))%*%rnorm(k)
  s2s[iter]   = s2
  bs[iter,]   = t(b)
}
s2s    = s2s[seq(burn+1,niter,by=LAG)]
bs     = bs[seq(burn+1,niter,by=LAG),1:k]
s2 <- s2s
b <-bs

par(mfrow=c(k+1,3))
for (i in 1:k){
  ts.plot(bs[,i],ylab="",xlab="iteration",main=paste("beta",i,sep=""))
  abline(h=beta[i],col=2)
  acf(bs[,i],ylab="",main="")
  hist(bs[,i],prob=T,xlab="",ylab="",main="")
  abline(v=beta[i],col=2)
}

  ts.plot(s2,ylab="",xlab="iteration",main=paste("sig2",sep=""))
  abline(h=sigsq,col=2)
  acf(s2,ylab="",main="")
  hist(s2,prob=T,xlab="",ylab="",main="")
  abline(v=sigsq,col=2)




par(mfrow=c(2,2))
for (i in 1:k) { plot(density(b[,i]),col="blue", main = cbind("Beta",i))
 }
plot(density(s2),col="blue", main = "Sigma sq."  )
title("Sigma sq.")

bmean = colMeans(b)
vb = diag(var(b))
seb = sqrt(vb)
ci = cbind(rep(0,k),rep(0,k),rep(0,k),rep(0,k),rep(0,k),rep(0,k))
for (i in 1:k) {
ci[i,] = quantile(b[,i], probs = c(0.5, 2.5, 5, 95, 97.5, 99.5)/100) }
cat("probability intervals for parameters")
rbind(c(0.5, 2.5, 5, 95, 97.5, 99.5),ci)
l95 = ci[,2]
u95 = ci[,5]


ypred = X%*%bmean
resid = y - ypred
apr = cbind(y,ypred,resid)
par(mfrow=c(1,1))
matplot(apr,col=1:k,type='l',main="Actual, Predicted, Residual")
rss = t(resid)%*%resid
tss = t(y - mean(y))%*%(y - mean(y))
R2 = 1 - rss/tss
aic = log(rss/n) + 2*(k+1)/n
bic = log(rss/n) + log(n)*(k+1)/n
df = n-k
ll = c("R2","aic","bic","rss","tss")
Coefficient.Estimates = cbind(bmean,seb,l95,u95)
Regression.Statistics = rbind(ll,round(cbind(R2,aic,bic,rss,tss),4))

s2m <- mean(s2)
ses2 <- sd(s2)
s2nt = quantile(s2, probs = c(2.5,97.5)/100) 
cat("0.95 probability interval for sigma sq.")
s2nt
Sigmasq.Estimates = cbind(s2m,ses2,s2nt[1],s2nt[2])


Coefficient.Estimates
Sigmasq.Estimates
Regression.Statistics





