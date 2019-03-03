# BeroulliAtBat.R
# bernoulli trials inference, baseball example

# informative beta prior for top 50 player batting ave

# Draw 1,000,000 values of theta from the prior and plot
# a <- 6   ### try different values, see fig. 2.4, p.20 of Hahn
# b <- 18   ### 6 and 18 look fine to me!
# a <- 12
# b <- 28


shape1 <- 12
shape2 <- 28
thdraws0 <- rbeta(1000000, shape1, shape2)
hist(thdraws0,col=4,breaks=100,freq=F)
lines(density(thdraws0),col=3,lwd=2)

# Uniform with a = 0.1 and b = 0.5
thetauab <- runif(1000000,min=0.1,max=0.5)
hist(thetauab,col=3)

# mean and SD of prior
mntheta <- mean(thdraws0)
mntheta; sd(thdraws0)

# median of posterior and other percentiles
intervals <- quantile(thdraws0, probs = c(0.01, 0.025, 0.5, 0.975, 0.99))
medianth <- intervals[3]
medianth
intervals



# data for baseball player Charlie Blackmon, April 2016
df <-  read.csv("Blackmonv2.csv", header = TRUE)
str(df)

atbats <- df$AB # has total for month at end
hits <- df$H

n <- length(atbats)


cbind(hits,atbats)
sum(hits)
sum(atbats)

# First five games
atbats5 <- atbats[1:5]
hits5 <- hits[1:5]
atbats5; hits5


batave5 <- sum(hits5)/sum(atbats5)
batave5

bataveapr <- sum(hits)/sum(atbats)
bataveapr

# running batting average
batave <- rep(0,n)
for (t in 1:n) {
batave[t] <- sum(hits[1:t])/sum(atbats[1:t])
}

games <- 1:n
plot(games,batave,type='l',col=3, main="Running batting average over 21 games")

# Question 5
h <- sum(hits5)
ab <- sum(atbats5)
# posterior with uniform prior is beta
# post5u <- theta^h * (1 - theta)^(ab-h)

# uniform is beta with a=b=1
a <- 1/2
b <- 1/2

# Draw 1,000,000 values of theta from the posterior
shape1 <- a+h
shape2 <- b+ab-h
thdrawsj <- rbeta(1000000, shape1, shape2)
plot(density(thdrawsj),col=3)
abline(h=1,col=2)

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

# Question 6
a <- 1
b <- 1

# Draw 1,000,000 values of theta from the posterior
shape1 <- a+h
shape2 <- b+ab-h
thdraws <- rbeta(1000000, shape1, shape2)

plot(density(thdrawsj),col=4)
lines(density(thdraws),col=2)


lines(density(thdraws0),col=3)



# Question 8

# Draw 1,000,000 values of theta from the posterior
h <- sum(hits)
ab <- sum(atbats)

shape1 <- a+h
shape2 <- b+ab-h
thdrawstot <- rbeta(1000000, shape1, shape2)
plot(density(thdrawstot),col=2,lwd=2, "red=beta prior posterior, black=uniform prior post")
lines(density(thdraws),col=4)
lines(density(thdraws0),col=3)

mean(thdrawstot)
sd(thdrawstot)
quantile(thdrawstot, probs = c(0.01, 0.025, 0.5, 0.975,0.99))


# posterior from uniform prior

shape1 <- 1+h
shape2 <- 1+ab-h
thdrawstotu <- rbeta(1000000, shape1, shape2)
lines(density(thdrawstotu),col=1,lwd=2)



######################
n <- 10
y <- 2
n2 <- 100
y2 <- 20


theta <- (0:1000)/1000

# beta posterior for a uniform prior


post <- (gamma(n+2)/(gamma(y+1)*gamma(n-y+1)))*(theta^y)*(1-theta)^(n-y)
post2 <- (gamma(n2+2)/(gamma(y2+1)*gamma(n2-y2+1)))*(theta^y2)*(1-theta)^(n2-y2)
postsame <- dbeta(theta,y+1,n-y+1)

posts <- cbind(post,postsame,post2)

matplot(theta,posts,type='l',col=c("green","blue"))
abline(h=1,col=2) # prior

# prob ave > 0.2
pgreater02 <- 1 - pbeta(0.2,(y+1),(n-y+1))

# beta prior

a <- 3
b <- 9  # equivalent to observing 12 obs with 3 successes

betaprior <- (gamma(a+b)/(gamma(a)*gamma(b)))*(theta^(a-1))*(1-theta)^(b-1)
betaprior2 <- dbeta(theta,3,9)

plot(theta,betaprior2,type='l',col="black")
lines(theta,betaprior,type='l',col="orange")

# posterior using the beta prior (and binomial likelihood)
postbeta <- (gamma(n+a+b)/(gamma(y+a)*gamma(n-y+b)))*(theta^(y+a-1))*(1-theta)^(n-y+b-1)

postbeta2 <- dbeta(theta,y+a,n-y+b)

plot(theta,postbeta,type='l',col=2)
lines(theta,betaprior,type='l',col="orange")
lines(theta,post,type='l',col="yellow")
lines(theta,postbeta2,type='l',col="blue")


# frequentist analysis
# Naive frequentist analysis (treat theta as a mean)
n <- sum(atbats)
thetahat <- sum(hits)/sum(atbats)
seth <- sd(hits)/sqrt(n)

dfab <- n-1

critt5 <- abs(qt(0.05/2, dfab))
critt5

ci95u <- thetahat + seth*critt5
ci95l <- thetahat - seth*critt5

cbind(thetahat,seth,ci95l,ci95u)

# Recognizing that theta is a proportion
sethp <- sqrt(thetahat*(1-thetahat)/atbats)

ci95u <- thetahat + sethp*critt5
ci95l <- thetahat - sethp*critt5

cbind(thetahat,sethp,ci95l,ci95u)

# with Normal interval

ci95u <- thetahat + sethp*1.96
ci95l <- thetahat - sethp*1.96

cbind(thetahat,sethp,ci95l,ci95u)



# using a uniform prior
theta <- (0:1000)/1000
n <- atbats
y <- sum(hits)
postu <- (gamma(n+2)/(gamma(y+1)*gamma(n-y+1)))*(theta^y)*(1-theta)^(n-y)
plot(theta,postu,type='l',col=2,main="All the data")



# Posterior with quantiles using a uniform prior
postbunif <- dbeta(theta,sum(hits)+1,atbats-sum(hits)+1)
plot(theta,postbunif,type='l',col="blue")
abline(h=1,col=2)
psample <- rbeta(100000,y+1,n-y+1)
pumn <- mean(psample)
pusd <- sd(psample)
cbind(pumn,pusd)
quantile(psample, probs = c(0.01,0.05,0.95,0.99))

#  posterior with beta prior from qu. 1
a <- 2
b <- 3
pinfsample <- rbeta(100000,y+a,n-y+b)
pumn <- mean(pinfsample)
pusd <- sd(pinfsample)
cbind(pumn,pusd)
quantile(pinfsample, probs = c(0.01,0.05,0.95,0.99))




# using the first game's data only (assuming 4 at bats)
hit1 <- hit[1:4]
n1 <- 4
n <- n1
y <- sum(hit1)
postu1 <- (gamma(n+2)/(gamma(y+1)*gamma(n-y+1)))*(theta^y)*(1-theta)^(n-y)
plot(theta,postu,type='l',col=1,main="All the data")
lines(theta,postu1,type='l',col=2,main="1st game only")

# game 2
n2 <- 3
hit2 <- hit[5:7]
y <- sum(hit2)
postu2 <- (gamma(n+2)/(gamma(y+1)*gamma(n-y+1)))*(theta^y)*(1-theta)^(n-y)
plot(theta,postu,type='l',col=1,main="All the data")
lines(theta,postu1,type='l',col=2,main="1st game only")
lines(theta,postu2,type='l',col=2,main="2nd game only")


# qu. 1 "informative" beta prior
theta <- (0:1000)/1000
a <- 2
b <- 3
betaprior <- dbeta(theta,a,b)
plot(theta,betaprior,type='l',col="black")

s1 <- sum(hits[1:4])
n1 <- atbatspergame[1]

# qu. 2 using a uniform prior
postbunif <- dbeta(theta,s1+1,n1-s1+1)
plot(theta,postbunif,type='l',col="blue")
abline(h=1,col=2)
psample <- rbeta(100000,s1+1,n1-s1+1)
pumn <- mean(psample)
pusd <- sd(psample)
cbind(pumn,pusd)


# qu. 3 posterior with beta prior from qu. 1
postbeta1 <- dbeta(theta,s1+a,n1-s1+b)
plot(theta,postbeta1,type='l',col="black")
lines(theta,betaprior,type='l',col="red")

# best estimate of batting ave.
theorymean <- (s1+a)/(s1+a+n1-s1+b)
pa <- s1+a
pb <- n1-s1+b
sdnum <- pa*pb
sdden <- ((pa+pb)^2)*(pa+pb+1)
theorysd <- sqrt(sdnum/sdden)
cbind(theorymean,theorysd)

# qu. 5

# game 2
n2 <- atbatspergame[2]
s2 <- sum(hits[(n1+1):(n1+n2)])

postbeta2 <- dbeta(theta,s1+s2+a,(n1+n2)-(s1+s2)+b)
plot(theta,postbeta2,type='l',col="blue")
lines(theta,postbeta1,type='l',col="black")
lines(theta,betaprior,type='l',col="red")

psample <- rbeta(100000,s1+s2+a,(n1+n2)-(s1+s2)+b)
pmn <- mean(psample)
psd <- sd(psample)
cbind(pmn,psd)

# game 3
n3 <- atbatspergame[3]
s3 <- sum(hits[(n1+n2+1):(n1+n2+n3)])

postbeta3 <- dbeta(theta,s1+s2+s3+a,(n1+n2+n3)-(s1+s2+s3)+b)
plot(theta,postbeta3,type='l',col="green")
lines(theta,postbeta2,type='l',col="blue")
lines(theta,postbeta1,type='l',col="black")
lines(theta,betaprior,type='l',col="red")

psample <- rbeta(100000,s1+s2+s3+a,(n1+n2+n3)-(s1+s2+s3)+b)
pmn <- mean(psample)
psd <- sd(psample)
cbind(pmn,psd)


# game 4
n4 <- atbatspergame[4]
s4 <- sum(hits[(n1+n2+n3+1):(n1+n2+n3+n4)])

postbeta4 <- dbeta(theta,s1+s2+s3+s4+a,(n1+n2+n3+n4)-(s1+s2+s3+s4)+b)
plot(theta,postbeta4,type='l',col="purple")
lines(theta,postbeta3,type='l',col="green")
lines(theta,postbeta2,type='l',col="blue")
lines(theta,postbeta1,type='l',col="black")
lines(theta,betaprior,type='l',col="red")

psample <- rbeta(100000,s1+s2+s3++s4+a,(n1+n2+n3+n4)-(s1+s2+s3+s4)+b)
pmn <- mean(psample)
psd <- sd(psample)
cbind(pmn,psd)






#draw 100,000 obs from the posterior
postsample <- rbeta(100000,s1+a,n1-s1+b)
empiricalmean <- mean(postsample)
empsd <- sd(postsample)

cbind(theorymean,empiricalmean)
cbind(theorysd,empsd)


