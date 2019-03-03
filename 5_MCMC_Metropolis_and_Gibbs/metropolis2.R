# metropolis2.R

# metropolis.R obtained draws from a Normal. Here we obtain draws
# from a t density
# *** plots and variance calculation show it has fatter tails when
# *** v is small (say 5)

# simulating from a probability density using a random walk metropolis
# Two jumping proposal distributions considered, U[-s,s] and N(0,s^2).
# s is the tuning parameter in either case.

# Target distribution is a Student-t distribution
mu <- 1
sig2 <- 1  # if these are unknown, we need conditionals for these in the 
v <- 5    # the MCMC chain. 

## *** try different values of s ***
# For uniform: s= 0.5 leads to too much autocorr in chain
#                 1.5 is a bit better, but still too much autocorr
#                 3.0 looks better

# For normal:  s= 0.5 leads to too much autocorr in chain
#                 1.5 is a bit better, but still too much autocorr
#                 2.0 looks good
#                 3.0 seems acceptance rate a bit low
#                10.0 shows what you get with a low acceptance rate

s <- 2.0  # tuning parameter - changes variance of proposal moves

set.seed(2340)
n.iter <- 1 + 10000
burnin <- 200
thetadraw <- acc <- rep(0,n.iter)

# starting value theta(0)
thetadraw[1] <- 2.0

accprob <- runif(n.iter,0,1)  # acceptance probabilities, a(i)  

# uniform proposal move
propmove <- runif(n.iter,-s,s)
# normal proposal move
# propmove <- rnorm(n.iter,mean=0,sd=s)

for (i in 2:n.iter) {
proptheta <- thetadraw[i-1] + propmove[i]
propt <- (proptheta - mu)/sqrt(sig2)
currentt <- (thetadraw[i-1] - mu)/sqrt(sig2)
logacceptnum <- log((1 + (propt^2)/v)^(-(v+1)/2))
logacceptden <- log((1 + (currentt^2)/v)^(-(v+1)/2))
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

# generate draws from true distribution for comparison
thetas <- rnorm(10000,mean=mu,sd=sqrt(sig2))
lines(density(thetas),col=2,lwd=2)

mean(thetadraw[burnin:n.iter])
sd(thetadraw[burnin:n.iter])

par(mfrow=c(1,1))
hist(thetadraw[burnin:n.iter],freq=F,breaks=100,col=4)
lines(density(thetadraw[burnin:n.iter]),lwd=2)
lines(density(thetas),col=2,lwd=2)



# see metropolis.R for drawing from a Normal density 



