# triathlete example
# READ IN DATA 
df <- read.table("tri1.dat", header = F)
str(df)
summary(df)
# time = total time taken
# swim, bike, run = distances of each leg of the event
# yrmonth = yr and month of event, e.g 905 = May 1990
obs = df$V1
time = df$V2
swim = df$V3
bike = df$V4
run = df$V5
yrmonth = df$V6

cbind(time,swim,bike,run)

par(mfrow=c(2,2))
plot(swim,time)
plot(bike,time)
plot(run,time)

# OLS results (noninformative prior)
utri = lm(time~swim+bike+run)
summary(utri)

# Coding error in data originally - obs 4 recorded as 83.67 instead
# of 93.67 - which led to the following.
# OLS provides horrible results. The intercept would represent the 
# transition time between legs of the event (getting out of the water
# and on to the bike, and getting off the bike and out to the run)

# Informative prior Bayesian analysis

y = time
X = cbind(1,swim,bike,run)
k = ncol(X)
n = nrow(X)
In = diag(1,n)

cor(X)

# Prior hyperparameters
# variance v0 picked as 3*sqrt(v0) giving reasonable "95%" interval
# mean transition time = 2, var = 0.3
# mean swim time for 1 mile = 30 mins, var = 1
# mean bike time for 1 mile = 2.7 mins, var = 0.015
# mean run time for 1 mile = 8 mins, var = 0.2

# order of coeffs is intercept swim bike run
b0      = c(2,28,2.7,7.5)
Vb0    = c(0.3,1,0.015,0.15)*diag(1,k)
# more diffuse prior:
# Vb0   = diag(10^100,k)
#lambda = 0.1333
lambda = 1
nu = 3

ss0 = 1
iVb0 = solve(Vb0)*1.0

# plot priors
R = 1000
theta = rep(0,R)
nm = rep(0,R)
normal2 = rep(0,R)

# evaluate each distn from mean - 4*se to mean + 4*se
# so start at mean - 4*se + 8*se*inc, where inc goes from 0 to 1)
# then evaluate the distribution at those values
i = seq(1,1000)
par(mfrow=c(2,2))

btm = b0[1]
btv = Vb0[1,1]
bt = btm -4*sqrt(btv) + 8*sqrt(btv)*i/1000
pt = dnorm(bt,mean = btm, sd = sqrt(btv))
plot(bt,pt,type = 'l',lwd=2,col=4)
legend("topleft",legend=c("Prior for transition"),col=4,lty=1,lwd=2,bty="n",cex=1.1)

# priors for regression coefficients
bsm = b0[2]
bsv = Vb0[2,2]
bs = bsm -4*sqrt(bsv) + 8*sqrt(bsv)*i/1000
ps = dnorm(bs,mean = bsm, sd = sqrt(bsv))
plot(bs,ps,type = 'l',lwd=2,col=4)
legend("topleft",legend=c("Prior for swim"),col=4,lty=1,lwd=2,bty="n",cex=1.1)

bbm = b0[3]
bbv = Vb0[3,3]
bb = bbm -4*sqrt(bbv) + 8*sqrt(bbv)*i/1000
pb = dnorm(bb,mean = bbm, sd = sqrt(bbv))
plot(bb,pb,type = 'l',lwd=2,col=4)
legend("topleft",legend=c("Prior for bike"),col=4,lty=1,lwd=2,bty="n",cex=1.1)

brm = b0[4]
brv = Vb0[4,4]
br = brm -4*sqrt(brv) + 8*sqrt(brv)*i/1000
pr = dnorm(br,mean = brm, sd = sqrt(brv))
plot(br,pr,type = 'l',lwd=2,col=4)
legend("topleft",legend=c("Prior for run"),col=4,lty=1,lwd=2,bty="n",cex=1.1)

# try different values of lambda and nu
# nu = 3
# lambda = 1
# prior for regression error variance, sigma squared
a = nu/2
b = nu*lambda/2
sig2 = a - 4*sqrt(b) + 8*sqrt(b)*i/1000
sig = sig2  - min(sig2) + 0.00001
psi = dgamma(sig,shape = a, scale = b)
plot(sig,psi,type = 'l',lwd=2,col=4)
legend("topleft",legend=c("Prior for variance"),col=4,lty=1,lwd=2,bty="n",cex=1.1)

# Likelihood/OLS values
xxiols  = solve(t(X)%*%X)
bols = xxiols%*%(t(X)%*%y)
sseols = t(y - X%*%bols)%*%(y - X%*%bols)
s2ols = sseols/(n-k)
sebols = sqrt(s2ols*diag(xxiols))

# Marginal posterior distributions of sigma2 and beta
d0   = (nu+n)/2

# posterior means and SDs
XXI    = solve(t(X)%*%X+solve(Vb0))
b1   = XXI%*%(t(X)%*%y+solve(Vb0)%*%b0)
d1  = (nu*lambda+t(y)%*%y+t(b0)%*%solve(Vb0)%*%b0-t(b1)%*%solve(XXI)%*%b1)/2
s21   = d1/(d0+1)
sds2  = s21/sqrt(d0+2)
nu1lambda1=nu*lambda+t(y-X%*%b0)%*%(In-X%*%solve(solve(Vb0)+t(X)%*%X)%*%t(X))%*%(y-X%*%b0)
sdb=sqrt((nu1lambda1/(nu+n))*diag(XXI))

ym = mean(y)
tss = t(y-ym)%*%(y-ym)
yp = X%*%b1
res = y - yp
rss = t(res)%*%res
rsq = 1-rss/tss
aic = log(rss/n) + 2*(k+1)/n
bic = log(rss/n) + log(n)*(k+1)/n
df = n-k



# get following to work at some point
# mdat <- matrix(matout, nrow = 3, ncol=2, byrow=FALSE,
#               dimnames = list(c("constant", "beta", "sigmasq"),
#                               c("mean", "s.e.")))
# print(mdat)


matout = round(rbind(cbind(b1,sdb),c(s21,sds2)),5)
cat(" posterior Mean (col 1) & Std Dev (col 2) for beta (rows 1 to k) & sigma (row k+1)",fill=TRUE)
print(matout)

cat(" R-squared, AIC, BIC")
print(cbind(rsq,aic,bic))

ress = cbind(utri$residuals, res)
ypred = cbind(y,utri$fitted.values,yp)
par(mfrow=c(2,2))
matplot(ress,type='l',lwd=2,col=1:2)
par(mfrow=c(1,1))
matplot(ypred,type='l',lwd=2,col=1:3,lty=1)
legend("topleft",legend=c("Actual","OLS","Bayes"),col=1:3,lty=1,lwd=2,bty="n",cex=1.1)

# compare uninformative/OLS and informative prior coefficients
cbind(bols,sebols,b1,sdb)

# likelihood for regression coefficients


# transition time coefficient
bsm = b0[1]
bsv = Vb0[1,1]
bslm = bols[1]
bslv = sebols[1]
bspm = b1[1]
bspv = sdb[1]
bsp = bspm -4*bspv + 8*bspv*i/1000
bsl = bslm -4*bslv + 8*bslv*i/1000
ps = dnorm(bsl,mean = bsm, sd = sqrt(bsv))
psl = dnorm(bsl,mean = bslm, sd = bslv)
psp = dnorm(bsl,mean = bspm, sd = bspv)
psy = cbind(ps,psl,psp)
matplot(bsl,psy,type = 'l',lwd=2,col=1:3, lty=1, main = "Transition time coefficient")
legend("topleft",legend=c("Prior", "Likelihood", "Posterior"),col=1:3,lty=1,lwd=2,bty="n",cex=1.1)


par(mfrow=c(2,2))
# transition time coefficient
bsm = b0[1]
bsv = Vb0[1,1]
bslm = bols[1]
bslv = sebols[1]
bspm = b1[1]
bspv = sdb[1]
bsp = bspm -4*bspv + 8*bspv*i/1000
bsl = bslm -4*bslv + 8*bslv*i/1000
ps = dnorm(bsp,mean = bsm, sd = sqrt(bsv))
psl = dnorm(bsp,mean = bslm, sd = bslv)
psp = dnorm(bsp,mean = bspm, sd = bspv)
psy = cbind(ps,psl,psp)
matplot(bsp,psy,type = 'l',lwd=2,col=1:3, lty=1, main = "Transition time coefficient")
legend("topleft",legend=c("Prior", "Likelihood", "Posterior"),col=1:3,lty=1,lwd=2,bty="n",cex=1.1)

# swim time coefficient
bsm = b0[2]
bsv = Vb0[2,2]
bslm = bols[2]
bslv = sebols[2]
bspm = b1[2]
bspv = sdb[2]
bsp = bspm -4*bspv + 8*bspv*i/1000
bsl = bslm -4*bslv + 8*bslv*i/1000
ps = dnorm(bsp,mean = bsm, sd = sqrt(bsv))
psl = dnorm(bsp,mean = bslm, sd = bslv)
psp = dnorm(bsp,mean = bspm, sd = bspv)
psy = cbind(ps,psl,psp)
matplot(bsp,psy,type = 'l',lwd=2,col=1:3, lty=1, main = "Swim coefficient")
legend("topleft",legend=c("Prior", "Likelihood", "Posterior"),col=1:3,lty=1,lwd=2,bty="n",cex=1.1)

# bike time coefficient
bsm = b0[3]
bsv = Vb0[3,3]
bslm = bols[3]
bslv = sebols[3]
bspm = b1[3]
bspv = sdb[3]
bsp = bspm -4*bspv + 8*bspv*i/1000
bsl = bslm -4*bslv + 8*bslv*i/1000
ps = dnorm(bsp,mean = bsm, sd = sqrt(bsv))
psl = dnorm(bsp,mean = bslm, sd = bslv)
psp = dnorm(bsp,mean = bspm, sd = bspv)
psy = cbind(ps,psl,psp)
matplot(bsp,psy,type = 'l',lwd=2,col=1:3, lty=1, main = "Bike coefficient")
legend("topleft",legend=c("Prior", "Likelihood", "Posterior"),col=1:3,lty=1,lwd=2,bty="n",cex=1.1)


# run time coefficient
bsm = b0[4]
bsv = Vb0[4,4]
bslm = bols[4]
bslv = sebols[4]
bspm = b1[4]
bspv = sdb[4]
bsp = bspm -4*bspv + 8*bspv*i/1000
bsl = bslm -4*bslv + 8*bslv*i/1000
ps = dnorm(bsp,mean = bsm, sd = sqrt(bsv))
psl = dnorm(bsp,mean = bslm, sd = bslv)
psp = dnorm(bsp,mean = bspm, sd = bspv)
psy = cbind(ps,psl,psp)
matplot(bsp,psy,type = 'l',lwd=2,col=1:3, lty=1, main = "Run coefficient")
legend("topleft",legend=c("Prior", "Likelihood", "Posterior"),col=1:3,lty=1,lwd=2,bty="n",cex=1.1)





# plot(bsl,psl,type = 'l',lwd=2,col=4, main = "Swim coefficient")
# legend("topleft",legend=c("Likelihood for swim"),col=4,lty=1,lwd=2,bty="n",cex=1.1)


# posteriors for regression coefficients
# bspm = b1[2]
# bspv = sdb[2]
# bsp = bsm -4*bspv + 8*bspv*i/1000
# psp = dnorm(bsp,mean = bspm, sd = bspv)
psy = cbind(ps,psl,psp)
bx = cbind(bsp,bsp,bsp)
matplot(bx,psy,type = 'l',lwd=2,col=4, main = "Swim coefficient")
legend("topleft",legend=c("Prior", "Likelihood", "Posterior"),col=4,lty=1,lwd=2,bty="n",cex=1.1)

bbm = b0[3]
bbv = Vb0[3,3]
bb = bbm -4*sqrt(bbv) + 8*sqrt(bbv)*i/1000
pb = dnorm(bb,mean = bbm, sd = sqrt(bbv))
plot(bb,pb,type = 'l',lwd=2,col=4)
legend("topleft",legend=c("Prior for bike"),col=4,lty=1,lwd=2,bty="n",cex=1.1)

brm = b0[4]
brv = Vb0[4,4]
br = brm -4*sqrt(brv) + 8*sqrt(brv)*i/1000
pr = dnorm(br,mean = brm, sd = sqrt(brv))
plot(br,pr,type = 'l',lwd=2,col=4)
legend("topleft",legend=c("Prior for run"),col=4,lty=1,lwd=2,bty="n",cex=1.1)

# try different values of lambda and nu
# nu = 3
# lambda = 1
# prior for regression error variance, sigma squared
a = nu/2
b = nu*lambda/2
sig2 = a - 4*sqrt(b) + 8*sqrt(b)*i/1000
sig = sig2  - min(sig2) + 0.00001
psi = dgamma(sig,shape = a, scale = b)
plot(sig,psi,type = 'l',lwd=2,col=4)
legend("topleft",legend=c("Prior for variance"),col=4,lty=1,lwd=2,bty="n",cex=1.1)



# Now use MCMC!
# load package bayesm
library(bayesm)

# rng seed
set.seed(12357)

# n = obs, R = iterations of MCMC
n=10
R=5000

# Simulate data
#X=cbind(rep(1,n),runif(n),runif(n)); beta=c(1,2,3); sigsq=0.5
#y=X%*%beta+rnorm(n,sd=sqrt(sigsq))
#k = ncol(X)
# Prior = list(b0,Vb0, nu, lambda)
# generate MCMC draws using bayesm routine 'runireg'

o=runireg(Data=list(y=y,X=X),Prior=list(betabar=b0,A=iVb0,nu=nu,ssq=ss0),Mcmc=list(R=R))

b = o$betadraw
s2 = o$sigmasqdraw

summary(b)
summary(s2)
sb = summary(b)


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

Coefficient.Estimates
Regression.Statistics

yy = cbind(y,ypred)
matplot(yy,type='l')

#RMSEs
yprior = X%*%b0
# yp = posterior predictions from analytical
# ypred = posterior predictions from MCMC
ypols = X%*%bols
yy = cbind(y,yprior,ypols,yp,ypred)
matplot(yy,type='l',col=1:5,lty=1,lwd=2)

rmse = rep(0,4)
for (i in 1:4) {
rmse[i] = sqrt(sum((y-yy[,i+1])^2)/n)
}
rmse





# log functional form?

lswim = log(swim)
lbike = log(bike)
lrun = log(run)

# OLS results (noninformative prior)
utri = lm(time~lswim+lbike+lrun)
summary(utri)


# prior in log form? Genereate a sample from nonlog prior
bsm = b0[3]
bsv = Vb0[3,3]
bslm = bols[3]
bslv = sebols[3]
bspm = b1[3]
bspv = sdb[3]
bsp = bspm -4*bspv + 8*bspv*i/1000
bsl = bslm -4*bslv + 8*bslv*i/1000
bssim = rnorm(1000,mean = bsm, sd = sqrt(bsv))
lbssim = log(bssim)
hist(lbssim)
hist(bssim)

mean(bssim)
sd(bssim)
mean(lbssim)
sd(lbssim)


