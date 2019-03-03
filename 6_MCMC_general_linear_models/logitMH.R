# logitMH.R
# logit model

# generate some data
set.seed(1235)
x = cbind(rep(1,100),rnorm(100),runif(100))
z = x%*%c(1.0, 2.0, 3.0) + 0.5*rnorm(100)
summary(z)

p <- exp(z)/(1+exp(z))

c <- 0.95

y <- ifelse(p > c, 1, 0)

# syntax: 
# BayesLogistic(y, x, steps = 1000,
#                     priorMean = NULL, priorVar = NULL,
#                     mleMean = NULL, mleVar,
#                     startValue = NULL, randomSeed = NULL,
#                     plots = FALSE)


xx <- x[,2:3]
# startValue = c(1.0, 2.0, 3.0),
mhb <- BayesLogistic(y, xx, steps = 5000, randomSeed = 1235)

b <- mhb$beta
matplot(b,type='l')
plot(b[,1],type='l')
plot(b[,2],type='l')

par(mfrow = c(2, 2), pty = "s")
acf(b[,1])
acf(b[,2])
acf(b[,3])

numEff(b[,2])
mean(b[,2])
sd(b[,2])
plot(density(b[,2]),col=3)
hist(b[,2],breaks=100,col=4)



numEff(b[,1])
mean(b[,1])
sd(b[,1])
plot(density(b[,1]),col=3)
hist(b[,1],breaks=100,col=4)




R = 1000
Data1=list(y=y,X=xx,p=2); Mcmc1=list(R=R,keep=1)
out=rmnlIndepMetrop(Data=Data1,Mcmc=Mcmc1)

cat("Summary of beta draws",fill=TRUE)
summary(out$betadraw,tvalues=beta)

if(0){
## plotting examples
plot(out$betadraw)
}

#
df = read.table("chdage.dat", header = TRUE)
str(df)
summary(df)
age <-df$age
chd <-df$chd

mhb <- BayesLogistic(chd, age, steps = 5000, randomSeed = 1235)
b <-mhb$beta

plot(b[,1],type='l')
plot(b[,2],type='l')

par(mfrow = c(2, 2), pty = "s")
acf(b[,1])
acf(b[,2])
acf(b[,3])

numEff(b[,2])
mean(b[,2])
sd(b[,2])
plot(density(b[,2]),col=3)
hist(b[,2],breaks=100,col=4)

numEff(b[,1])
mean(b[,1])
sd(b[,1])
plot(density(b[,1]),col=3)
hist(b[,1],breaks=100,col=4)

nParameters <- 2
beta.df = data.frame(b)
    names(beta.df) = paste("b", 0:(ncol(beta.df) - 1), sep = "")
    describe(beta.df)
    Mean.beta = sapply(beta.df, mean)
    StdDev.beta = sapply(beta.df, sd)
    Z.beta = Mean.beta/StdDev.beta
    print(data.frame(Mean.beta, StdDev.beta, Z.beta))
         nRows = nParameters
        nCols = 2
        oldPar = par(mfrow = c(nRows, nCols))
        nms = names(beta.df)
        for (i in 1:nParameters) {
            plot(ts(beta.df[, i]), main = paste("Time series plot of", 
                nms[i]), ylab = nms[i])
            plot(acf(beta.df[, i], plot = FALSE), main = paste("Autocorrelation plot of", 
                nms[i])) }




# bayesm RW MH

library(bayesm)
R = 1000
y <- as.vector(chd)
X <- as.matrix(age)
lgtdata=NULL
lgtdata[[1]]=list(y=y,X=X)
Data=list(lgtdata=lgtdata)

out=rhierBinLogit(Data=list(lgtdata=lgtdata),Mcmc=list(R=R))


summary(out)

b <- out$betadraw
vb <- out$betadraw
d <- out$Deltadraw
summary(b)
summary(vb)
summary(d)

hist(b,breaks=100,col=3)
acf(b)
b <- as.numeric(out$betadraw)

# horrible acceptance rate!
plot(b,type='l',col=4)

# sbeta = 0.01 works much better
R <- 5000
out=rhierBinLogit(Data=list(lgtdata=lgtdata),Mcmc=list(sbeta=0.01,R=R))

b <- as.numeric(out$betadraw)
plot(b,type='l',col=4)
acf(b)
hist(b,breaks=100,col=3)
plot(density(b),type='l',col=4)





# sbeta = 0.01 works much better
R <- 10000
out=rhierBinLogit(Data=list(lgtdata=lgtdata),Mcmc=list(sbeta=0.005,R=R,keep=5))

b <- as.numeric(out$betadraw)
plot(b,type='l',col=4)
acf(b)
hist(b,breaks=100,col=3)
plot(density(b),type='l',col=4)









cat("Summary of Delta draws",fill=TRUE)
summary(out$Deltadraw,tvalues=as.vector(Delta))
cat("Summary of Vbeta draws",fill=TRUE)
summary(out$Vbetadraw,tvalues=as.vector(Vbeta[upper.tri(Vbeta,diag=TRUE)]))
