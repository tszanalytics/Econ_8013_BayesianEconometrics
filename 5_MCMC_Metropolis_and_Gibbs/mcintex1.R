
niter <- 10000
x <- runif(niter)
y <- runif(niter)*5

plot(x,y,col=8)

i <- 0:niter/niter
i <- i[1:niter]
z <- dbeta(i,shape1=8,shape2=16)

lines(i,z,col=2)

count <- ifelse(y<=z,1,0)
area <- sum(count)/niter
areaofbox = 1*5
area * areaofbox

length(z)


# Gibbs sampling!

niter <- 10000
set.seed(123)
y <- rnorm(200,mean=10,sd=17)

ymean <- mean(y)
s2 <- var(y)
n <- length(y)
v <- n-1
# storage for my random draws
tau <- rep(0,niter)
mu <- rep(0,niter)

# initial value
tau[1] <- 4
mu[1] <- rnorm(1,mean=ymean,sd=(1/sqrt(n*tau[1])))

for (i in 2:niter) {

mu[i] <- rnorm(1,mean=ymean,sd=(1/sqrt(n*tau[i-1])))
tau[i] <- rgamma(1,n/2, 0.5*(v*s2 + n*(mu[i] - ymean)^2))

}

mean(mu); sd(mu)
plot(mu,type='l',col=3)

hist(mu,freq=F)
lines(density(mu),col=4)

sigma2 <- 1/tau

mean(sigma2); sd(sigma2)
plot(sigma2,type='l',col=3)

hist(sigma2,freq=F)
lines(density(sigma2),col=4)

sqrt(mean(sigma2))

# Metropolis for a standard normal

niter <- 10000

mu <- rep(0,niter)
mu[1] <- 10

accprob <- runif(niter,0,1)  # acceptance probability



for (i in 2:niter) {
prop <- mu[i-1] + rnorm(1,mean=0,sd=1)
propratio <- dnorm(prop,0,1)/dnorm(mu[i-1],0,1)

mu[i] <- ifelse(propratio > accprob[i],prop,mu[i-1])

}

plot(mu,type='l')
mean(mu)
sd(mu)
# approx 0.95 interval
mean(mu)+2*sd(mu)
mean(mu)-2*sd(mu)

quantile(mu,prob=c(0.025,0.975))


#### Now put some other likelihood x prior in instead of normal

# Metropolis for a some function
# substitute in a function in place of dnorm below (in loop)

niter <- 10000

mu <- rep(0,niter)
mu[1] <- 10

accprob <- runif(niter,0,1)  # acceptance probability



for (i in 2:niter) {
prop <- mu[i-1] + rnorm(1,mean=0,sd=1)
propratio <- dnorm(prop,0,1)/dnorm(mu[i-1],0,1)

mu[i] <- ifelse(propratio > accprob[i],prop,mu[i-1])

}

plot(mu,type='l')

mean(mu)
sd(mu)
# approx 0.95 interval
mean(mu)+2*sd(mu)
mean(mu)-2*sd(mu)

quantile(mu,prob=c(0.025,0.975))










