# sixsideddicerolls.R
# Roll m six-sided fair dice (start with m = 1 and increase m)

set.seed(123458)

par(mfrow=c(3,2))

# m = number of dice rolled
m <- 1
# n = number of times rolled
n <- 10
draws <- rep(0,n*m)
dim(draws) <- c(n,m)
for (i in 1:m) {
draws[,i] <- runif(n)
}

pr <- rep(0,6*m)

ps0 <- rep(0,m)
p1 <- (draws<=1/6)*1
p2 <- ((draws<=2/6) - p1)*2
p3 <- ((draws<=3/6) - (draws<=2/6))*3
p4 <- ((draws<=4/6) - (draws<=3/6))*4
p5 <- ((draws<=5/6) - (draws<=4/6))*5
p6 <- (draws>=5/6)*6
rolls <- p1 + p2 + p3 + p4 + p5 + p6
srolls <- apply(rolls,1,sum)
# b <- 10*m + 1
# hist(srolls,freq=F,breaks=b)

ff <- tabulate(srolls)/n   # also works without the 6 defining number of bins
plot(ff,type = "h", lwd=4, col=4,ylab="data frequencies",ylim=c(0,(max(ff)+0.02)))
abline(a=1/6,b=0)

