# Examples from Albert (2009) Bayesian Computation with R
# 

# dicrete example from the book
p = seq(0.05, 0.95, by = 0.1)
prior = c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0)
prior = prior/sum(prior)
plot(p, prior, type = "h", lwd=3, ylab="Prior Probability")

#
# six sided die problem
# x = no. of dots observed
x = seq(1, 6)
priorx = c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6)
plot(x, priorx, type = "h", lwd=3, col=4,ylab="Prior Probability")

# suppose we roll the die 10 times and get the sequence
# 6325314314
# the posterior distribution is now
r = c(6,3,2,5,3,1,4,3,1,4)

ff = tabulate(r,6)   # also works without the 6 defining number of bins
plot(x,ff,type = "h", lwd=3, col=4,ylab="data frequencies")


# normalized likelihood
likn = ff/sum(ff)

# posterior
post = ff*priorx/sum(ff*priorx)

plot(x,post,type = "h", lwd=3, col=4,ylab="Posterior probabilities")

# If we want to assign a prior that indicates we are fairly confident that 1/6
# is actually the right value for each probability (we have good reason to 
# think the die is actually fair), then we assign a hyperprior to each value
# of x
# Assign a beta distribution to each value of x and update this
# or create a pseudo-sample to represent the informative prior

# priorinform is equivalent to 100 observations from a fair die
priorinform = c(10,10,10,10,10,10) # each value observed 10 times

# priorvinf is equivalent to 1000 observations from a fair die
priorvinf = c(100,100,100,100,100,100) # each value observed 100 times


postinf = (ff+priorinform)/sum(ff+priorinform)

postvinf = (ff+priorvinf)/sum(ff+priorvinf)


pp = cbind(post,postinf,postvinf)
pp

par(mfrow=c(1,1))
matplot(x,pp,type = "h", lwd=c(15,9,3), lty=1, col=c(3,4,2),ylab="Posterior probabilities")
legend("topright",c("uniform","Informative (10)","Very inf. (100)"),lty=1,lwd=c(15,9,3), col=c(3,4,2))



par(mfrow=c(2,2))
matplot(x,pp,type = "h", lwd=c(4,2,1), lty=1, col=c(3,4,2),ylab="Posterior probabilities")
plot(x,post,type = "h", lwd=3, col=4,main="Posterior uninformative prior")
plot(x,postinf,type = "h", lwd=3, col=4,main="Posterior informative prior")
plot(x,postvinf,type = "h", lwd=3, col=4,main="Posterior very informative prior")



# VARIANCE OF ESTIMATE OF PROPORTION = p(1-p)/n
varunif = 0.1*(1-0.1)

varinf = postinf[1]*(1-postinf[1])/10
varvinf = postvinf[1]*(1-postvinf[1])/100

sdinf = sqrt(varinf)
sdvinf = sqrt(varvinf)
sdu = sqrt(varunif)
cbind(sdu,sdinf,sdvinf)

par(mfrow=c(2,3))
# posteriors for x = i (normal approx.)
for (i in 1:6) {

sdinf = sqrt(postinf[i]*(1-postinf[i])/10)
sdvinf = sqrt(postvinf[i]*(1-postvinf[i])/100)
curve(dnorm(x,mean=postvinf[i],sd=sdvinf),col=1,lty=1,lwd=2)
curve(dnorm(x,mean=postinf[i],sd=sdinf),add=T,col=2,lty=1,lwd=2)
curve(dnorm(x,mean=post[i],sd=sdu),add=T,col=3,lty=1,lwd=2)
}

# some other stuff below - playing around with Beta distributions

# mean of beta = a/(a+b), var = m(1-m)/(a+b+1) - 
# WARNING: THIS DOES NOT WORK FOR PROPORTION (as above)
# NEED TO USE BINOMIAL DISTN INSTEAD - mean goes to zero for beta!

par(mfrow=c(1,1))
# s = no. of successes in f trials
s = 2; f = 10

unifm = 1/2
p10m = (1+s)/(1+s + 10+f)
p100m = (1+s)/(1+s + 100+f)



# posterior with a very informative prior (100 pseudo obs.)
a = 1; b = 100
curve(dbeta(x,a+s,b+f), lty=1,lwd=3,col=2,from=0, to=1,xlab="p",ylab="Density")
# posterior with an informative prior (10 pseudo obs.)
a = 1; b = 10
curve(dbeta(x,a+s,b+f), add=TRUE, lty=1,lwd=3,col=1)
# posterior with flat prior observing 2 x=1 out of 10 trials gives, for x = 1
a = 1; b = 1  # give uniform prior with a beta density
curve(dbeta(x,a+s,b+f),add=TRUE, lty=1,lwd=3,col=3)
# flat prior
curve(dbeta(x,a,b),add=TRUE,lty=1,lwd=3,col=4)


# Simulating draws
draws = runif(1000)
x = draws <=0.1
sum(x)
# mean of binomial = np, var = np(1-p)
n = 1000
p = sum(x)/n

mdraw = n*p
vardraw = n*p*(1-p)

# VARIANCE OF ESTIMATE OF PROPORTION = p(1-p)/n
varp = p*(1-p)/n

# Example from chapter 2 of Albert (2009)
# prior, likelihood and posterior BETA densities
a = 3.26
b = 7.19
s = 11
f = 16
curve(dbeta(x,a+s,b+f), from=0, to=1,xlab="p",ylab="Density",lty=1,lwd=4)
curve(dbeta(x,s+1,f+1),add=TRUE,lty=2,lwd=4)
curve(dbeta(x,a,b),add=TRUE,lty=3,lwd=4,col=4)
legend(.7,4,c("Prior","Likelihood","Posterior"),lty=c(3,2,1),lwd=c(3,3,3))

# 95% posterior confidence interval
qbeta(c(0.05, 0.95), a + s, b + f)

# simulate a sample of 1000 obs from the posterior
ps = rbeta(100000, a + s, b + f)
hist(ps,xlab="p",breaks=100,main="")

# probability that x (a proportion) is greater than 0.5
sum(ps >= 0.5)/100000

# 95% posterior confidence interval from simulation
quantile(ps, c(0.05, 0.95))
