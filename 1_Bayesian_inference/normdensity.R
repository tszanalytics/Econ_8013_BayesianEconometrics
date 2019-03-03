#
# Plot normal density with mean = m and variance = s2
# 

m = 0.5
s2 = 10
s22 = 30

R = 10000 # number of evaluations
i = seq(1,R)
x = m -4*sqrt(s2) + 8*sqrt(s2)*i/R
x2 = m -4*sqrt(s22) + 8*sqrt(s22)*i/R
px = dnorm(x,mean = m, sd = sqrt(s2))


plot(x,px,type = 'l',lwd=2,col=4)
legend("topleft",legend=c("Normal density"),col=4,lty=1,lwd=2,bty="n",cex=1.1)

# Now do it using the formula

normal1 = (1/sqrt(2*pi*s2))*exp(-((x-m)^2)/(2*s2))
plot(x,normal1,type = 'l',lwd=2,col=4)

pp = cbind(px,normal1)
matplot(x2,pp,,type = 'l',lwd=2,col=c(1,4))
legend("topleft",legend=c("Identical"),col=4,lty=1,lwd=2,bty="n",cex=1.1)

# Now change the variance of one of them

normal2 = (1/sqrt(2*pi*s22))*exp(-((x-m)^2)/(2*s22))

pp = cbind(px,normal2)
matplot(x,pp,,type = 'l',lwd=2,col=c(1,4),lty=1)
legend("topleft",legend=c("small var","large var"),col=c(1,4),lty=1,lwd=2,bty="n",cex=1.1)

# Using the larger variance as scale for horizontal axis and shift mean:
m2 = 5.0
px = dnorm(x2,mean = m2, sd = sqrt(s2))
px2 = dnorm(x2,mean = m, sd = sqrt(s22))
pp = cbind(px,px2)
matplot(x2,pp,,type = 'l',lwd=2,col=c(1,4),lty=1)
legend("topleft",legend=c("small var","large var"),col=c(1,4),lty=1,lwd=2,bty="n",cex=1.1)





