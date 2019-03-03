# MetropRegHahnp85-metroptau.R
# *** includes metropolis step for tau also ***

# Linear regression using a random walk metropolis

# All we need is prior x likelihood

# see comments on following function in Hahn(2014), p.85
# prior*likelihood function:
logpost <<- function(b0,b1,tau) {
logpost1 <- ((n/2)-1)*log(tau)
resid <- y - b0 - b1*xdif
logpost2 <- -(tau/2) * sum(resid^2)
return(logpost1 + logpost2)
}

# data
y <- c(6,7,8,12,3,6,9,8)
x <- rep(c(1,2,3,4),2)
n <- length(y)
xdif <- x - mean(x)



# Metropolis algorithm
set.seed(12347)
n.iter <- 10000 
burnin <- 500
b0length <- 5; b1length <- 5  ### tuning parameters
taulength <- 2.0

tau <- b0 <- b1 <- rep(0,n.iter)

# initial values
tau[1] <- 0.5
b0[1] <- 0
b1[1] <- 0

b0accprob <- runif(n.iter,0,1)  # acceptance probs
b1accprob <- runif(n.iter,0,1)
tauaccprob <- runif(n.iter,0,1)
b0acc <- b1acc <- tauacc <- 0  # tuning check

# uniform RW moves
b0move <- runif(n.iter,0,b0length) - b0length/2
b1move <- runif(n.iter,0,b1length) - b1length/2
taumove <- runif(n.iter,0,taulength) - taulength/2

# MC loop
for (i in 2:n.iter) {
  if(i%%10000==0) cat(i, " iterations complete", "\n")

# sample tau using Gibbs - try Metropolis for this step also!
  prop <- tau[i-1] + taumove[i]
  if(prop<=0) prop = tau[i-1] + abs(taumove[i])  # use antithetic value of prop <= 0
  propratio <- exp(logpost(b0[i-1],b1[i-1],prop) - 
                   logpost(b0[i-1],b1[i-1],tau[i-1]))

  if(propratio > tauaccprob[i]) {
    tau[i] <- prop     # accept proposal
    tauacc <- tauacc + 1  } # count acceptances
  else {
    tau[i] <- tau[i-1]  } # reject proposal


# sample b0
  prop <- b0[i-1] + b0move[i]
  propratio <- exp(logpost(prop,b1[i-1],tau[i]) - 
                   logpost(b0[i-1],b1[i-1],tau[i]))

  if(propratio > b0accprob[i]) {
    b0[i] <- prop     # accept proposal
    b0acc <- b0acc + 1  } # count acceptances
  else {
    b0[i] <- b0[i-1]  } # reject proposal

# sample b1
  prop <- b1[i-1] + b1move[i]
  propratio <- exp(logpost(b0[i],prop,tau[i]) - 
                   logpost(b0[i],b1[i-1],tau[i]))

  if(propratio > b1accprob[i]) {
    b1[i] <- prop     # accept proposal
    b1acc <- b1acc + 1  } # count acceptances
  else {
    b1[i] <- b1[i-1]  } # reject proposal

}  # end of MCMC loop

# eliminate burnin 
b0 <- b0[burnin:n.iter]
b1 <- b1[burnin:n.iter]
tau <- tau[burnin:n.iter]

# posterior means and SDs
header <- cbind("param","mean","SD")
params <- rbind("b0","b1","tau")
results <- rbind(round(cbind(mean(b0),sd(b0)),digits=3),
             round(cbind(mean(b1),sd(b1)),digits=3),round(cbind(mean(tau),sd(tau)),digits=3))
results2 <- rbind(header,cbind(params,results))


# acceptance rates
b0acc/n.iter
b1acc/n.iter
tauacc/n.iter

noquote(results2)  # remove the quotation marks from the output


# plot b0
par(mfrow=c(3,3))
plot(b0,type='l')
acf(b0,col=2,lwd=4)
hist(b0,freq=F,breaks=100,col=4)
lines(density(b0),lwd=2)
abline(v=mean(b0),col=3,lwd=3) # posterior mean
abline(v=median(b0),col=2,lwd=3) # posterior median

# plot b1
#par(mfrow=c(3,1))
plot(b1,type='l')
acf(b1,col=2,lwd=4)
hist(b1,freq=F,breaks=100,col=4)
lines(density(b1),lwd=2)
abline(v=mean(b1),col=3,lwd=3) # posterior mean
abline(v=median(b1),col=2,lwd=3) # posterior median

# plot tau
#par(mfrow=c(3,1))
plot(tau,type='l')
acf(tau,col=2,lwd=4)
hist(tau,freq=F,breaks=100,col=4)
lines(density(tau),lwd=2)
abline(v=mean(tau),col=3,lwd=3) # posterior mean
abline(v=median(tau),col=2,lwd=3) # posterior median

############################


