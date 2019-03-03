# normmixMHbolstad.R

# MH example from Bolstad (2010)


# look at function normMixMH for coding of RW and ind. MH algorithms

library(Bolstad2)
## Set up the normal mixture
theta0 <- c(0,1) # density 1, mean and S.D.
theta1 <- c(3,2) # density 2
p <- 0.8  # mixture proportion

# syntax
# normMixMH(theta0, theta1, p,  candidate, steps = 1000, type = 'ind',
#            randomSeed = NULL, startValue = NULL)



## Sample from an independent N(0,3^2) candidate density
# candidate <- c(0, 3) # try different S.Ds to get different rejection rates, etc.

candidate <- c(0, 50)  # ********  try SD = 50, 9 and 0.5 - 9 is about right! ******

MCMCsampleInd <- normMixMH(theta0, theta1, p, candidate, steps = 5000, startValue = 5,randomSeed=12357)

y <- MCMCsampleInd
m <- length(y)

rej <- 0
for (i in 2:m) {
d <- ifelse(y[i] == y[i-1], 1, 0)
rej <- rej + d
}
rej
rrate <- rej/m
rrate

par(mfrow = c(1, 2), pty = "s")
hist(y,breaks=50,col=3)
acf(y,lwd=3)

# calculate values to plot true density
 mu0 <- theta0[1]
    sigma0 <- theta0[2]
    mu1 <- theta1[1]
    sigma1 <- theta1[2]
    mu <- candidate[1]
    sigma <- candidate[2]


 theta <- seq(from = min(mu0 - 3 * sigma0, mu1 - 3 * sigma1), 
        to = max(mu0 + 3 * sigma0, mu1 + 3 * sigma1), by = 0.001)
    fx <- p * dnorm(theta, mu0, sigma0) + (1 - p) * dnorm(theta, 
        mu1, sigma1)

# plot MCMC sample and true densities
par(mfrow = c(1, 2), pty = "s")
x <-MCMCsampleInd
hist(x, prob=T,breaks=100,col=5) # light blue/turquoise
lines(theta, fx,col=4) # blue
plot(density(x),col=2) # red
summary(x)
lines(theta, fx,col=4) # blue

# first 50 draws
x[1:50]


# RW MH algorithm
## If we wish to use the alternative random walk N(0, 0.5^2)
## candidate density
# candidate <- c(0, 0.5) # do same as above - try different values

candidate <- c(0, 3) # try 50, 9, 0.5 and 2 (2 or 3 looks best)

MCMCsampleRW <- normMixMH(theta0, theta1, p, candidate, type = 'rw',steps = 5000, startValue = 0,randomSeed=12357)



y <- MCMCsampleRW
m <- length(y)
par(mfrow = c(1, 2), pty = "s")
hist(y,breaks=50,col=3)
acf(y,lwd=3)
rejrate(y)
summary(y)


# plot true density
# ptheta <- p*exp(-0.5*theta^2) + (1-p)*exp(((theta-theta1[1])^2)/(2*theta1[2]^2))

 mu0 <- theta0[1]
    sigma0 <- theta0[2]
    mu1 <- theta1[1]
    sigma1 <- theta1[2]
    mu <- candidate[1]
    sigma <- candidate[2]


 theta <- seq(from = min(mu0 - 3 * sigma0, mu1 - 3 * sigma1), 
        to = max(mu0 + 3 * sigma0, mu1 + 3 * sigma1), by = 0.001)
    fx <- p * dnorm(theta, mu0, sigma0) + (1 - p) * dnorm(theta, 
        mu1, sigma1)


targetSample <- MCMCsampleRW

oldPar <- par(mfrow = c(1, 2), pty = "s")
    h <- hist(targetSample, plot = FALSE)
    ymax <- max(c(h$density, fx)) * 1.05
    hist(targetSample, prob = TRUE, breaks=80, col = "light blue", xlim = range(theta), 
        ylim = c(0, ymax), main = "Sample from target density", 
        xlab = "x", ylab = "Density")
    lines(theta, fx)
    box()
    plot(targetSample, type = "l", main = "", ylab = "Target Sample")
    par(oldPar)
    invisible(targetSample)


plot(theta,fx,type='l')



