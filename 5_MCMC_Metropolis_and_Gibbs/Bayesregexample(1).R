# example of single equation regression iid MCMC draws
# see ch. 2 Rossi-Allenby-McCulloch (2005)

# load package bayesm
library(bayesm)

# rng seed
set.seed(12357)

# n = obs, R = iterations of MCMC
n=200
R=50000

# Simulate data
X=cbind(rep(1,n),runif(n),runif(n)); beta=c(1,2,3); sigsq=0.5
y=X%*%beta+rnorm(n,sd=sqrt(sigsq))
k = ncol(X)

# generate MCMC draws using bayesm routine 'runireg'
o=runireg(Data=list(y=y,X=X),Mcmc=list(R=R))

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

bmeans = colMeans(b)
vb = diag(var(b))
seb = sqrt(vb)
ci = cbind(rep(0,k),rep(0,k),rep(0,k),rep(0,k),rep(0,k),rep(0,k))
for (i in 1:k) {
ci[i,] = quantile(b[,i], probs = c(0.5, 2.5, 5, 95, 97.5, 99.5)/100) }
cat("probability intervals for parameters")
rbind(c(0.5, 2.5, 5, 95, 97.5, 99.5),ci)
l95 = ci[,2]
u95 = ci[,5]


ypred = X%*%bmeans
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
Coefficient.Estimates = cbind(bmeans,seb,l95,u95)
Regression.Statistics = rbind(ll,round(cbind(R2,aic,bic,rss,tss),4))

Coefficient.Estimates
Regression.Statistics

coeffs = Coefficient.Estimates
regstats = Regression.Statistics

regout = list(coeffs,regstats)
regout



# uniformative aka using OLS routine
r5 <-lm(y~X-1)

tstats <- summary(r5)$coeff[,3]
tstats

v <- summary(r5)$df[2]
v

# precise vs. precise odds
odds <- (1 + (tstats^2)/v)^(-(v+1)/2)
odds
oddsagainstzero <- 1/odds

cbind(summary(r5)$coeff,oddsagainstzero)

# Odds using MCMC sample

ts <- bmeans/seb
v <- n - k
odds <- (1 + (ts^2)/v)^(-(v+1)/2)
odds
oddsagainstzero <- 1/odds
cbind(coeffs,oddsagainstzero)



# to test coeff number 2, b2 = 2.0 for example
ts <- abs(bmeans[2]-2.0)/seb[2]
v <- n - k
odds <- (1 + (ts^2)/v)^(-(v+1)/2)
odds
oddsagainsttwo <- 1/odds
oddsagainsttwo

