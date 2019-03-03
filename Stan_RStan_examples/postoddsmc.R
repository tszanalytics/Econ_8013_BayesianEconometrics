# postoddsmc.R


postoddsmc <- function(mcsample,null=0.0){

# use: set postoddsmc(mcsample,null=nullvalue)
#      mcsample = sample from the posterior
#      nullvalue = value under the null hypothesis (to be tested)

# purpose: evaluate the posterior odds against a point null hypothesis
#          give a sample from the posterior density (from MC or MCMC)

# output = odds ratio vs. null of zero 
#          posterior "p-values" for left and right areas around null value
#          from an MC sample from the posterior density

# Jeff Mills and Hamed Namavari, 2016

testPoint = null       #### null hypothesis cut-off point

DenB1 <- density(mcsample,n=1024)  #### posterior density

# par(mfrow=c(1 ,1))
# plot(DenB1$x,DenB1$y, type="l",main=expression('H'[0]),ylab="density", xlab=expression(paste(beta[2])),col=4,lwd=1)
# abline(v=testPoint,col=2,lwd=2)

index <- which(DenB1$x-testPoint==min(abs(DenB1$x-testPoint)))
if (length(index) == 0) {
  index <- which(-DenB1$x+testPoint==min(abs(DenB1$x-testPoint)))
}

# an alternative to above?
# if (testPoint<min(DenB1$x) | testPoint>max(DenB1$x)) num = min(DenB1$y)

denom <-  max(DenB1$y)
num <- DenB1$y[index]
oddsratio <- denom/num


####  NEXT LINE TO HAVE oddsratio > 10000 reported ######
if(oddsratio > 10000) {oddsratio ="> 10000"}


# cat("Odds Against Null",fill=TRUE)
# oddsratio

list(oddsratio = oddsratio)
}


