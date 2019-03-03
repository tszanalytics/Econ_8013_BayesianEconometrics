# tdraws.R

# random draws from the t-distribution
# Let t = (x - mx)/s, s = (SD(x)/sqrt(n)) 
# so x = s*t + mx

t1 <- rt(1000,1)
t3 <- rt(1000,3)
t5 <- rt(1000,5)
t10 <- rt(1000,10)

par(mfrow = c(2, 2), pty = "s")
hist(t1,breaks=100)
hist(t3,breaks=100)
hist(t5,breaks=100)
hist(t10,breaks=100)

# For a variable with mean = 5 and SD = 2
s <- 2
mx <- 5
x1 <- s*t1 + mx
x3 <- s*t3 + mx
x5 <- s*t5 + mx
x10 <- s*t10 + mx
par(mfrow = c(2, 2), pty = "s")
hist(x1,breaks=100,col=2)
hist(x3,breaks=100,col=3)
hist(x5,breaks=100,col=4)
hist(x10,breaks=100,col=5)


