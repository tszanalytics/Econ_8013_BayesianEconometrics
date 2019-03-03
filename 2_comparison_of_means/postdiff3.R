# function to compare three distributions

postdiff3 = function(sm1,ssd1,sm2,ssd2,sm3,ssd3,n1,n2,n3,R) {

# 
#  purpose: posterior comparison of means analysis
#
#  sm1 = sample 1 mean, sm2 = sample 2 mean, sm3 = sample 3 mean
#  ssd1 = sample 1 SD, ssd2 = sample 2 SD, ssd3 = sample 3 SD
#  n1 = sample 1 size, n2 = sample 2 size, n3 = sample 3 size
#  R = number of MC draws

t7.5 = rt(R, df=(n1-1))   # sample 1
mc7.5 = sm1 + t7.5*ssd1/sqrt(n1)
t15 = rt(R, df=(n2-1))    # sample 2
mc15 = sm2 + t15*ssd2/sqrt(n2)
n = 14
t30 = rt(R, df=(n3-1))    # sample 3
mc30 = sm3 + t30*ssd3/sqrt(n3)

# posteriors for each dose level
#plot(density(mc7.5),col=4,xlim=c(2,7),ylim=c(0,2.2),main="Posterior half-life for child-adol.-adult",xlab="half-life for dose 30")
#legend("topright",legend=c("Sample 1","Sample 2","Sample 3"),col=c(4,3,2),lty=1,lwd=2,bty="n",cex=1.1)
#lines(density(mc15),col=3)
#lines(density(mc30),col=2)

#print("summary stats for samples 1, 2 and 3:")
#mean(mc7.5); sd(mc7.5)
#quantile(mc7.5,probs=c(0.005,0.025,0.5,0.975,0.995))
#mean(mc15); sd(mc15)
#quantile(mc15,probs=c(0.005,0.025,0.5,0.975,0.995))
#mean(mc30); sd(mc30)
#quantile(mc30,probs=c(0.005,0.025,0.5,0.975,0.995))

# differences
#d7515 = mc7.5 - mc15
#d1530 = mc15 - mc30
#d7530 = mc7.5 - mc30


#quantile(d7515,probs=c(0.005,0.025,0.5,0.975,0.995))
#quantile(d1530,probs=c(0.005,0.025,0.5,0.975,0.995))
#quantile(d7530,probs=c(0.005,0.025,0.5,0.975,0.995))

#plot(density(d7530),col=4,xlim=c(-3,3),ylim=c(0,1.3),main="Posterior difference in half-life for child-adol.-adult",xlab="Change in half-life for dose 30")
#legend("topright",legend=c("child-adult","adolescent-adult","child-adolescent"),col=c(4,3,2),lty=1,lwd=2,bty="n",cex=1.1)
#lines(density(d1530),col=3)
#lines(density(d7515),col=2)
#abline(v=0)

result = list(sdraws1 = mc7.5, sdraws2 = mc15, sdraws3 = mc30)
return(result)
}


