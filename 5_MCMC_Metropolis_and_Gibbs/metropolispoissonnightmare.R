# metropolispoissonnightmare.R

# simulating from a probability density using a random walk metropolis
# Two jumping proposal distributions considered, U[-s,s] and N(0,s^2).
# s is the tuning parameter in either case.

# Target distribution is a poisson likelihood x prior
# If prior is a gamma, then posterior is gamma
# - see Hahn (2014), p. 52-55

mu <- 1
sig2 <- 1  # if these are unknown, we need conditionals for these in the 
           # the MCMC chain. 

## *** try different values of s ***
# For uniform: s= 0.5 leads to too much autocorr in chain
#                 1.5 is a bit better, but still too much autocorr
#                 3.0 looks better

# For normal:  s= 0.5 leads to too much autocorr in chain
#                 1.5 is a bit better, but still too much autocorr
#                 2.0 looks better
#                 3.0 seems acceptance rate a bit low
#                10.0 shows what you get with a low acceptance rate

s <- 0.08 # tuning parameter - changes variance of proposal moves

set.seed(2340)
n.iter <- 1 + 500000
burnin <- 1000
thetadraw <- acc <- rep(0,n.iter)

sumx <- 10
n <- 100
a <- 30

# starting value theta(0)
thetadraw[1] <- 0.1

accprob <- runif(n.iter,0,1)  # acceptance probabilities, a(i)  

# uniform proposal move
propmove <- runif(n.iter,-s,s)

# normal proposal move
# propmove <- rnorm(n.iter,mean=0,sd=s)

for (i in 2:n.iter) {
proptheta <- thetadraw[i-1] + propmove[i]
# if (proptheta < 0) proptheta = thetadraw[i-1] + runif(1,0,s)
if (proptheta < 0) proptheta = abs(proptheta)

logacceptnum <- -n*proptheta + sumx*log(proptheta) + proptheta*log(a) # - lgamma(proptheta[i]+1)
logacceptden <- -n*thetadraw[i-1] + sumx*log(thetadraw[i-1]) + thetadraw[i-1]*log(a) # - lgamma(thetadraw[i-1]+1)

# logacceptnum <- -((proptheta - mu)^2) / (2*sig2)
# logacceptden <- -((thetadraw[i-1] - mu)^2) / (2*sig2)
densratio <- exp(logacceptnum - logacceptden)
  if (accprob[i] < densratio) {
thetadraw[i] <- proptheta # accept proposal
acc[i] <- 1 }  # keep count of acceptances
else 
thetadraw[i] <- thetadraw[i-1] # reject proposal
}

# acceptance rate
sum(acc)/n.iter

par(mfrow=c(3,1))
plot(thetadraw,type='l')
acf(thetadraw,col=2,lwd=4)
hist(thetadraw[burnin:n.iter],freq=F,breaks=100,col=4)
lines(density(thetadraw[burnin:n.iter]),lwd=2)

mean(thetadraw[burnin:n.iter])
sd(thetadraw[burnin:n.iter])


burnin <- 1000
par(mfrow=c(3,1))
plot(thetadraw[burnin:n.iter],type='l')
acf(thetadraw[burnin:n.iter],col=2,lwd=4)
hist(thetadraw[burnin:n.iter],freq=F,breaks=100,col=4)
lines(density(thetadraw[burnin:n.iter]),lwd=2)

# Take a look at beginning of MCMC chain
start <- 1
end <- 100
par(mfrow=c(3,1))
plot(thetadraw[start:end],type='l')
acf(thetadraw[start:end],col=2,lwd=4)
hist(thetadraw[start:end],freq=F,breaks=100,col=4)
lines(density(thetadraw[start:end]),lwd=2)

# plot prior*likelihood function
n <- 100
sumx <- 10
a <- 30
proptheta <- 0 + 0.5*(1:1000/1000)
logacceptnum <- rep(0,length(proptheta))
for (i in 1:1000) {
logacceptnum[i] <- -n*proptheta[i] + sumx*log(proptheta[i]) + proptheta[i]*log(a) # - lgamma(proptheta[i]+1)
}
plot(proptheta,exp(logacceptnum),type='l')

