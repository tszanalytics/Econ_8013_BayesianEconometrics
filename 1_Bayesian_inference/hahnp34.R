# hahn, p.34 example

# since we have just propn to likelihood, we can drop the 5 from Hahn's formula

set.seed(2346)
n.iter = 10000
theta = runif(n.iter)
p = runif(n.iter)
ind = ifelse(p < (theta^4)*(1-theta),1,0)
mean(ind)

th = 1:1000/1000
pth = (th^4)*(1-th)
plot(theta,p,pch='.',col="red")
lines(th,pth,type='l',ylim=c(0,1),col="blue")

# NB: (1/6)*(1/5)= 0.03333333



# problem with the above is that the acceptance rate is very low
# so multiplying by some value actually increases the accuracy.
set.seed(6981)
n.iter = 10000  # 100000 max or it will give a solid block of red!
theta = runif(n.iter)
p = runif(n.iter)
ind = ifelse(p < 10*(theta^4)*(1-theta),1,0)
mean(ind)
pth = 10*(th^4)*(1-th)
plot(theta,p,pch='.',col="red")
lines(th,pth,type='l',ylim=c(0,1),col="blue")




# an acceptance rate of between about 0.25 and 0.4 works best 
# in practice.

# p. 35: A more sophisticated way to get the same result is
# based on the mean value theorem
# = evaluating the area under the curve by partitioning the line and 
#  computing the area of each rectangle under the curve.
#  This is much more efficient (far fewer iterations required)

set.seed(42)
n.iter = 1000
theta = runif(n.iter)
ptheta = (theta^4)*(1-theta)
mean(ptheta)

# we don't really need the randomness here - we can partition into
# equal spaces, but it can be easier to code the random draws.

# plot of posterior using the random values of theta
plot(theta,ptheta,col="blue")  # would have to order them to plot as line

# doing the same with a sequence (very easy in this case)
theta = 1:1000/1000
ptheta = (theta^4)*(1-theta)
mean(ptheta)
lines(theta,ptheta,type='l',col="magenta")


### p.36 The powerful approach!

# Using the law of large numbers for numerical accuracy of approximations.
# Formula (3.3), p. 36 is remarkably powerful.  We can use ANY function of y

# All we need to be able to do is to obtain random (independent) draws 
# from the distribution for y, then we can use (3.3) to obtain the 
# posterior distribution and moments, etc., for any f(y).


# Example from p. 37
# Suppose we have two variables from two different distributions 
# They don't have to even be the same type of distribution - see further below)
# In the example, both distributions are binomial, just with different parameters
set.seed(8831)
n.iter = 10000
pi.w = rbeta(n.iter,15,24)  # draws from women's distn.
pi.m = rbeta(n.iter,4,20)  # draws from men's distn.

plot(density(pi.w),col=1,xlim=c(0,0.7),ylim=c(0,6))
lines(density(pi.m),col=2)

# summary stats
mean(pi.w); sd(pi.w)
mean(pi.m); sd(pi.m)

# Now we can look at the posterior distn. for the difference in pi.w - pi.m
# THINK about how you might do this using a frequentist approach (good luck!).

difpi = pi.w - pi.m

hist(difpi,breaks=100,freq=F,col=4)
lines(density(difpi),col="green",lwd=2)
abline(v=0,col=2) # where is zero (no difference null hypothesis)?

# summary stats
mean(difpi); sd(difpi)
summary(difpi)

# "confidence" aka "credible" aka probability intervals
q95 = quantile(difpi,probs = c(0.025,0.975))  # 95%
q99 = quantile(difpi,probs = c(0.005,0.995))  # 99%
q95; q99

# further out in the tail (as with the 99% interval), the more draws
# needed for accuracy

# lines for the 95% and 99% intervals
abline(v=q95,col="orange")
abline(v=q99,col="magenta")
abline(v=mean(difpi),col=1)

# note the skewness to the left tail so the intervals are not symmetric 
# about the mean

mean(difpi) - q95[1]
mean(difpi) - q95[2]
mean(difpi) - q99[1]
mean(difpi) - q99[2]

# This simulation based Bayesian approach automatically takes this into
# account.
# Standard frequentist intervals would be (incorrectly) symmetric about
# the mean.

# Now how to get the predictive distribution using this posterior
# what is the following line doing?  see p.39 and eqn. (2.7)
sales.w = rbinom(n.iter,500,pi.w)  # assuming 500 female customers

# plot the predictive distribution
hist(sales.w,breaks=100,col="green",freq=F)
lines(density(sales.w),col="blue",lwd=2)

mean(sales.w); sd(sales.w)
summary(sales.w)
q95s = quantile(sales.w,probs = c(0.025,0.975))  # 95%
q95s

# we obtain distributions of any function of the random variable
# by applying the function TO THE RANDOM DRAWS from the distribution

# trivial example, (p.39) if revenue per item = $5.20 then:
rev.per.item = 5.2

sales.rev.w = rev.per.item * sales.w

hist(sales.rev.w,breaks=100,col="green",freq=F)
lines(density(sales.rev.w),col="blue",lwd=2)

# To make it more interesting, revenue could be some nonlinear 
# function of total amount sold, etc.


# Now draw from two normals and compare means
# [can draw from two completely different distributions (say beta and normal)]

m = 10000  # number of draw, i.e. numerical evaluations for accuracy

samplemean1 = 0.6
samplesd1 = 1.0
n1 = 40

samplemean2 = 1.2
samplesd2 = 0.7
n2 = 35

# see formula (3.23), p. 60 Hahn - completing the square
# for the likelihood
s1draws = rnorm(m, mean=samplemean1,sd=samplesd1/sqrt(n1))
s2draws = rnorm(m, mean=samplemean2,sd=samplesd2/sqrt(n2))

plot(density(s1draws),col=3,xlim=c(0.0,2.0),ylim=c(0.0,6.0))
lines(density(s2draws),col=4)
# plot(density(s2draws),col=4)

# posterior distribution of difference in means
diffdraws = s2draws - s1draws
plot(density(diffdraws),col=3)
mean(diffdraws)
sd(diffdraws)
qs = quantile(diffdraws,probs=c(0.005,0.025,0.5,0.975,0.995))
qs
abline(v=qs,col=4)

# comparison of means from beta and normal
# one sample from N(0.9.0,0.7), 35 obs
# second sample from Beta(17,20) (so n = 37)
s1draws = rnorm(m, mean=0.9,sd=0.7/sqrt(35)) # draws from normal distn.
plot(density(s1draws),col=3,xlim=c(0.0,2.0),ylim=c(0.0,6.0))
s2draws = rbeta(m,17,20)  # draws from beta distn
lines(density(s2draws),col=4)

# posterior means computed in exactly same way as previously
diffdraws = s1draws - s2draws
plot(density(diffdraws),col=3)
mean(diffdraws)
sd(diffdraws)
qs = quantile(diffdraws,probs=c(0.005,0.025,0.5,0.975,0.995))
qs
abline(v=qs,col=4)


