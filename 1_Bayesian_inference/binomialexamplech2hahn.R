# binomialexamplech2hahn.R

# Plot, prior, likelihood and posterior
# Find normalizing constant analytically and numerically

# Suppose likelihood is binomial and prior is beta (or uniform if a=b=1)
# posterior is then beta:

# Define prior parameters and data values
a = 1
b = 1
y = 20
n = 50

# create a range of values for theta 
# (over all reasonable possible values for theta)
theta <- 0:1000/1000

# Plot the likelihood over possible values of y for given theta:
thetal <- 0.3
ys <- 0:n
py <- (factorial(n)/(factorial(n-ys)*factorial(ys)))*(thetal^ys)*(1-thetal)^(n-ys)

plot(ys,py,col=4,pch=20)

# analytical evaluation of normalizing constant,c, for posterior
lc <- lgamma(a+b+n) - lgamma(a+y) - lgamma(b+n-y)
c <- exp(lc)
c

ptheta <- c*(theta^(a+y-1)) * (1-theta)^(b+n-y-1)

plot(theta,ptheta,type='l',col=4,lwd=2, main="posterior = blue, prior = red")
abline(h=1,col=2)
# since the prior = 1 for all theta, the posterior is proportional to 
# the likelihood

# numerical evaluation of normalizing constant, c
# = integral of numerator w.r.t. theta
int1 <- function(th) (th^(a+y-1)) * (1-th)^(b+n-y-1)
cth = integrate(int1, lower = 0, upper = 1.0)
# c is the normalizing constant = 1/integral 
c = 1/cth$value
c

# Use simulation to evaluate mean, mode and SD
# Draw 1,000,000 values of theta from the posterior
shape1 <- a+y
shape2 <- b+n-y
thdraws <- rbeta(1000000, shape1, shape2)
lines(density(thdraws),col=3)

# mean of posterior
mntheta <- mean(thdraws)
mntheta

# mode of posterior
postheta <- density(thdraws)
indp <- which.max(postheta$y)
modeth <- postheta$x[indp]
modeth

# median of posterior and other percentiles
intervals <- quantile(thdraws, probs = c(0.01, 0.025, 0.5, 0.975, 0.99))
medianth <- intervals[3]
medianth
intervals

# standard deviation and variance
var(thdraws)
sd(thdraws)

# Analytical calculations to check:
# see p. 22-23 Hahn
modetheta <- (a+y-1)/(a+b+n-2)
meantheta <- (a+y)/(a+b+n)
vartheta <- ((a+y)*(b+n-y)) / (((a+b+n)^2) * (a+b+n+1))

modetheta; meantheta; vartheta; sqrt(vartheta)


########## Hahn, p. 34 example
set.seed(234)
iters <- 100000  # try different number of draws
thetas <- runif(iters)
pthetas <- runif(iters)
indic <- ifelse(pthetas < 5*thetas^4 * (1-thetas),1,0)
mean(indic)

# We are using the law of large numbers here
# Consider the mean as iters increases
iters <- 1000
means <- rep(0,iters)
for (i in 1:iters) {
means[i] <- mean(indic[1:i])
}
plot(means,type='l')



