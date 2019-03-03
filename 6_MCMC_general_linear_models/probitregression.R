###################################################################
#
# Example 6.2.1 from Press (2003)
#
###################################################################
#
# Occurrence or nonoccurrence of infection following birth by 
# Ceasarian section.
#
# Response variable
#
#      y=1 if caesarian birth resulted in an infection
#      y=0 if caesarian birth DID NOT result in an infection
#
# Explanatory variables
#
#      x1=1 if caesarian was nonplanned
#      x2=1 if risk factors were present at the time of birth
#      x3=1 if antibiotics were given as a prophylaxis
#
###################################################################
# Auxiliary functions and densities
# ---------------------------------
like=function(beta){
  prob = pnorm(x%*%beta)
  return(prod((prob^y)*((1-prob)^(1-y))))
}
like1=function(x,beta){
  prob = pnorm(x%*%beta)
  return(prod((prob^y)*((1-prob)^(1-y))))
}
prior  = function(beta){prod(dnorm(beta,0,2.236068))}
rprior = function(n){matrix(rnorm(4*n,0,2.236068),n,4)}
post   = function(beta){like(beta)*prior(beta)}
post1  = function(x,beta){like1(x,beta)*prior(beta)}
q025   = function(x){quantile(x,0.025)}
q975   = function(x){quantile(x,0.975)}

# Dataset and exploratory analysis
# --------------------------------
n     = c(98,18,2,26,58,9,40)
count = c(11,1,0,23,28,0,8)
y = NULL
for (i in 1:7) y = c(y,rep(0,n[i]-count[i]),rep(1,count[i]))
x = matrix(c(rep(c(1,1,1),n[1]),rep(c(0,1,1),n[2]),rep(c(0,0,1),n[3]),
    rep(c(1,1,0),n[4]),rep(c(0,1,0),n[5]),rep(c(1,0,0),n[6]),
    rep(c(0,0,0),n[7])),ncol=3,byrow=T)
x = cbind(1,x)

par(mfrow=c(3,1))
plot(y,xlab="birth",ylab="",main="x1=1 if caesarian was nonplanned")
lines(x[,2],col=2)
plot(y,xlab="birth",ylab="",main="x2=1 if risk factors were present at the time of birth")
lines(x[,3],col=2)
plot(y,xlab="birth",ylab="",main="x3=1 if antibiotics were given as a prophylaxis")
lines(x[,4],col=2)

# Maximum likelihood estimates
# ----------------------------
bhat = c(-1.093022,0.607643,1.197543,-1.904739)
Vhat = matrix(0,4,4)
Vhat[1,1] =  0.040745
Vhat[2,2] =  0.073101
Vhat[3,3] =  0.062292
Vhat[4,4] =  0.080788
Vhat[1,2] = -0.007038
Vhat[1,3] = -0.039399
Vhat[1,4] =  0.004829
Vhat[2,3] = -0.006940
Vhat[2,4] = -0.050162
Vhat[3,4] = -0.016803
for (i in 2:4)
  for (j in 1:(i-1))
    Vhat[i,j] = Vhat[j,i]
cVhat = t(chol(Vhat))

# MCMC scheme (Random walk Metropolis)
# ------------------------------------
set.seed(1234)
beta  = bhat
tau   = 1
betas = NULL
burn  = 1000
M     = 1000
LAG   = 20
niter = burn+LAG*M
for (iter in 1:niter){
  betadraw = beta + tau*cVhat%*%rnorm(4)
  accept = min(1,post(betadraw)/post(beta))
  if (runif(1)<accept){
     beta = betadraw
  }
  betas = rbind(betas,t(beta))
}
betas = betas[seq(burn+1,burn+LAG*M,by=LAG),]

# Table 6.2, page 126, Press (2003)
# ---------------------------------
m=round(cbind(apply(betas,2,mean),
sqrt(apply(betas,2,var)),
apply(betas,2,q025),
apply(betas,2,q975)),3)
------------------------------------
beta     Mean   s.d.   Lower   Upper
------------------------------------
beta0  -1.080  0.213  -1.514  -0.685 
beta1   0.603  0.238   0.134   1.069 
beta2   1.177  0.256   0.708   1.689 
beta3  -1.897  0.262  -2.403  -1.400 
------------------------------------

# Figure 6.1, page 128, Press (2003)
# ----------------------------------
par(mfrow=c(2,4))
breaks = seq(-2,0,length=25)
hist(betas[,1],main=expression(beta[0]),ylab="",xlab="",
     breaks=breaks,prob=T,ylim=c(0,2))
breaks = seq(-0.5,2,length=25)
hist(betas[,2],main=expression(beta[1]),ylab="",xlab="",
     breaks=breaks,prob=T,ylim=c(0,2))
breaks = seq(0,2.5,length=25)
hist(betas[,3],main=expression(beta[2]),ylab="",xlab="",
     breaks=breaks,prob=T,ylim=c(0,2))
breaks = seq(-3,-1,length=25)
hist(betas[,4],main=expression(beta[3]),ylab="",xlab="",
     breaks=breaks,prob=T,ylim=c(0,2))
acf(betas[,1],xlab="lag",ylab="",main="",lag.max=20)
acf(betas[,2],xlab="lag",ylab="",main="",lag.max=20)
acf(betas[,3],xlab="lag",ylab="",main="",lag.max=20)
acf(betas[,4],xlab="lag",ylab="",main="",lag.max=20)


# Computing p(y|M) by using the MC estimator
# ------------------------------------------
ind     = matrix(0,7,4)
ind[1,] = c(1,2,3,4)
ind[2,] = c(1,2,3,0)
ind[3,] = c(1,2,4,0)
ind[4,] = c(1,3,4,0)
ind[5,] = c(1,2,0,0)
ind[6,] = c(1,3,0,0)
ind[7,] = c(1,4,0,0)
betap   = rprior(10000)
py.mc   = NULL
for (j in 1:7){
  ll = NULL
  for (i in 1:10000)
    ll = c(ll,like1(x[,ind[j,]],betap[i,ind[j,]]))
  py.mc = c(py.mc,mean(ll))
}
pm = py.mc/sum(py.mc)
cbind(ind,round(100*pm,2))
----------------------------------------
Model  X1   X2   X3   X4  100Pr(M|data)%
----------------------------------------
M1      1    2    3    4       73.01
M2      1    2    3    0        0.00
M3      1    2    4    0        0.01
M4      1    3    4    0       26.98
M5      1    2    0    0        0.00
M6      1    3    0    0        0.00
M7      1    4    0    0        0.00
----------------------------------------

# Posterior inference for model M4
# --------------------------------
set.seed(1234)
beta  = bhat[ind[4,]]
tau   = 1
betas1= NULL
burn  = 1000
M     = 1000
LAG   = 20
niter = burn+LAG*M
for (iter in 1:niter){
  cVhat = t(chol(Vhat[ind[4,],ind[4,]]))
  betadraw = beta + tau*cVhat%*%rnorm(3)
  accept = min(1,post1(x[,ind[4,]],betadraw)/post1(x[,ind[4,]],beta))
  if (runif(1)<accept){
     beta = betadraw
  }
  betas1 = rbind(betas1,t(beta))
}
betas1 = betas1[seq(burn+1,burn+LAG*M,by=LAG),]

m=round(cbind(apply(betas1,2,mean),
sqrt(apply(betas1,2,var)),
apply(betas1,2,q025),
apply(betas1,2,q975)),3)
------------------------------------
beta     Mean   s.d.   Lower   Upper
------------------------------------
beta0  -0.977  0.211  -1.406  -0.583
beta2   1.235  0.249   0.756   1.730
beta3  -1.526  0.217  -1.956  -1.111
------------------------------------

par(mfrow=c(2,2))
plot(density(betas1[,1]),xlab="",ylab="",main=expression(beta[0]))
lines(density(betas[,1]),lty=2)
plot(density(betas[,2]),xlab="",ylab="",main=expression(beta[1]),lty=2)
segments(0,0,0,1.5,lty=1)
plot(density(betas1[,2]),xlab="",ylab="",main=expression(beta[2]))
lines(density(betas[,3]),lty=2)
plot(density(betas1[,3]),xlab="",ylab="",main=expression(beta[3]))
lines(density(betas[,4]),lty=2)

# Bayesian model averaging (Models M1 and M4)
# -------------------------------------------
xx = matrix(c(1,1,1,0,1,1,0,0,1,1,1,0,0,1,0,1,0,0,0,0,0),7,3,byrow=T)
xx = cbind(1,xx)
prob  = pnorm(betas%*%t(xx))
prob1 = pnorm(betas1%*%t(xx[,ind[4,]]))
prob2 = 0.73*prob+0.27*prob1
par(mfrow=c(3,2))
for (i in c(1,2,7,6,5,4)){
  f = density(prob[,i])
  f1 = density(prob1[,i])
  f2 = density(prob2[,i])
  U = max(f$y,f1$y,f2$y)
  plot(f2,xlab="",ylab="",main="",xlim=c(0,1),ylim=c(0,U),axes=F)
  axis(1);axis(2)
  lines(f2,lty=1,lwd=1.5)
  lines(f,lty=2,lwd=1.5)
  lines(f1,lty=3,lwd=1.5)
  title(paste("x1=",xx[i,2], " - x2=",xx[i,3]," - x3=",xx[i,4],sep=""))
}

par(mfrow=c(1,1))
f2 = density(prob2[,1])
plot(f2,xlab="",ylab="",main="",xlim=c(0,1),ylim=c(0,25),axes=F)
axis(1);axis(2)
j=0
for (i in c(1,2,7,6,5,4)){
  j = j + 1
  f2 = density(prob2[,i])
  lines(f2,lty=1,col=j,lwd=1.5)
}
legend(0.4,20,legend=c(
"x1=1  x2=1  x3=1",
"x1=0  x2=1  x3=1",
"x1=0  x2=0  x3=0",
"x1=1  x2=0  x3=0",
"x1=0  x2=1  x3=0",
"x1=1  x2=1  x3=0"),col=1:6,lwd=rep(1.5,6),bty="n")
title("x1=1 - caesarian nonplanned 
       x2=1 - risk factors present birth
       x3=1 - antibiotics given")















