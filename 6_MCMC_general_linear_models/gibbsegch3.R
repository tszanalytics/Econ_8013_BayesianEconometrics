# Gibbs sampler for conjugate normal linreg

# set obs, iterations, seed
set.seed(66)
R=2000
n=100

# simulate data
X=cbind(rep(1,n),runif(n)); beta=c(1,2); sigsq=.25
y=X%*%beta+rnorm(n,sd=sqrt(sigsq)

# MCMC
Data1=list(y=y,X=X); Mcmc1=list(R=R) 

out=runiregGibbs(Data=Data1,Mcmc=Mcmc1)

# summarize output
cat("Summary of beta and Sigma draws",fill=TRUE)
summary(out$betadraw,tvalues=beta)
summary(out$sigmasqdraw,tvalues=sigsq)

# plots
if(0){
## plotting examples
plot(out$betadraw)
}
plot(out$sigmasqdraw)

# mean and sd of beta draws
sd(out$betadraw)
meanb = colSums(out$betadraw)/R
meanb

# one column at a time version
mean(out$betadraw[1:R,1])
sqrt(var(out$betadraw[1:R,1]))
sqrt(var(out$betadraw[1:R,2]))


