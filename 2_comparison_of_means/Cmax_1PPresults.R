# Cmax 1PP data

# READ in data and define variables: 

df = read.csv("Cmax_1PP.csv", header = TRUE) ### Steady state Cmax 1-PP data
str(df)
summary(df)


dose = df$Dose                 # dosages
cmcm = df$Cmax_CHILD.n.12.     # mean Cmax Child n = 12
cmcsd = df$Cmax_CHILD_SD       # SD Cmax Child n = 12
cmtm =  df$Cmax_ADOLES.n.12.   # mean Cmax Adolescent n = 12
cmtsd = df$Cmax_ADOLES_SD      # SD Cmax Adolescent n = 12
cmam =  df$Cmax_ADULTS.n.14.   # mean Cmax Adult n = 14
cmasd = df$Cmax_ADULTS_SD      # SD Cmax Adult n = 14


cbind(dose,cmcm,cmcsd,cmtm,cmtsd,cmam,cmasd)

# plots of posteriors for 3 sample means comparison
par(mfrow=c(1,2))
source(file="postdiff3.R")
# posteriors for each dose level

#### dose = 7.5
hflf = postdiff3(sm1=cmcm[1],ssd1=cmcsd[1],sm2=cmtm[1],ssd2=cmtsd[1],sm3=cmam[1],ssd3=cmasd[1],
              n1=12,n2=12,n3=14,R=100000)
mean(hflf$sdraws1); sd(hflf$sdraws1)
mean(hflf$sdraws2); sd(hflf$sdraws2)
mean(hflf$sdraws3); sd(hflf$sdraws3)

plot(density(hflf$sdraws1),col=4,xlim=c(1,12),ylim=c(0,0.9),main="Posterior Cmax for child-adol.-adult",cex.main=1.0,xlab="Cmax for dose 7.5")
legend("topright",legend=c("Child","Adolescent","Adult"),col=c(4,3,2),lty=1,lwd=2,bty="n",cex=0.8)
lines(density(hflf$sdraws2),col=3)
lines(density(hflf$sdraws3),col=2)


print("summary stats for Child, Adolescent and Adult:")
mean(hflf$sdraws1); sd(hflf$sdraws1)
quantile(hflf$sdraws1,probs=c(0.005,0.025,0.5,0.975,0.995))
mean(hflf$sdraws2); sd(hflf$sdraws2)
quantile(hflf$sdraws2,probs=c(0.005,0.025,0.5,0.975,0.995))
mean(hflf$sdraws3); sd(hflf$sdraws3)
quantile(hflf$sdraws3,probs=c(0.005,0.025,0.5,0.975,0.995))


# differences
ds1s2 = hflf$sdraws1 - hflf$sdraws2
ds2s3 = hflf$sdraws2 - hflf$sdraws3
ds1s3 = hflf$sdraws1 - hflf$sdraws3

cbind(mean(ds1s2), sd(ds1s2))
quantile(ds1s2,probs=c(0.005,0.025,0.5,0.975,0.995))
cbind(mean(ds2s3), sd(ds2s3))
quantile(ds2s3,probs=c(0.005,0.025,0.5,0.975,0.995))
cbind(mean(ds1s3), sd(ds1s3))
quantile(ds1s3,probs=c(0.005,0.025,0.5,0.975,0.995))


plot(density(ds1s2),col=4,xlim=c(-2.5,8),ylim=c(0,0.6),main="Posterior difference in 
             Cmax for child-adol.-adult",cex.main = 1.0,xlab="Change in Cmax for dose 7.5")
legend("topright",legend=c("child-adoles.","adoles.-adult","child-adult"),col=c(4,3,2),
                                                           lty=1,lwd=2,bty="n",cex=0.8)
lines(density(ds2s3),col=3)
lines(density(ds1s3),col=2)
abline(v=0)


# p-values and odds analytically - MUST ASSUME SAME SAMPLE SIZE FOR BOTH SAMPLES
# For differing sample sizes, do this numerically with MC draws
ts1s2 = mean(ds1s2)/sd(ds1s2)
ts2s3 = mean(ds2s3)/sd(ds2s3)
ts1s3 = mean(ds1s3)/sd(ds1s3)

rbind(cbind("       "," Child-Adoles."," Adoles.-Adult", "Child-Adult"),
cbind("p-value",round(1-pt(abs(ts1s2),df=11),3)*2, round(1-pt(abs(ts2s3),df=12),3)*2, round(1-pt(abs(ts1s3),df=12),3)*2),
cbind("odds vs. 0",round(dt(0,df=11)/dt(ts1s2,df=11),2),round(dt(0,df=12)/dt(ts2s3,df=12),2),
     round(dt(0,df=12)/dt(ts1s3,df=12),2)))

# 90% Credible intervals to help check p-values and odds
quantile(ds1s2,probs=c(0.05,0.95))
quantile(ds2s3,probs=c(0.05,0.95))
quantile(ds1s3,probs=c(0.05,0.95))



#### dose = 15
hflf = postdiff3(sm1=cmcm[2],ssd1=cmcsd[2],sm2=cmtm[2],ssd2=cmtsd[2],sm3=cmam[2],ssd3=cmasd[2],
              n1=12,n2=12,n3=14,R=100000)
mean(hflf$sdraws1); sd(hflf$sdraws1)
mean(hflf$sdraws2); sd(hflf$sdraws2)
mean(hflf$sdraws3); sd(hflf$sdraws3)

plot(density(hflf$sdraws1),col=4,xlim=c(4.5,20),ylim=c(0,0.5),main="Posterior Cmax: dose = 15", 
                                      xlab="Cmax for dose 15",cex.main = 1.0)
legend("topright",legend=c("Child","Adolescent","Adult"),col=c(4,3,2),lty=1,lwd=2,bty="n",cex=0.8)
lines(density(hflf$sdraws2),col=3)
lines(density(hflf$sdraws3),col=2)


print("summary stats for Child, Adolescent and Adult:")
mean(hflf$sdraws1); sd(hflf$sdraws1)
quantile(hflf$sdraws1,probs=c(0.005,0.025,0.5,0.975,0.995))
mean(hflf$sdraws2); sd(hflf$sdraws2)
quantile(hflf$sdraws2,probs=c(0.005,0.025,0.5,0.975,0.995))
mean(hflf$sdraws3); sd(hflf$sdraws3)
quantile(hflf$sdraws3,probs=c(0.005,0.025,0.5,0.975,0.995))


# differences
ds1s2 = hflf$sdraws1 - hflf$sdraws2
ds2s3 = hflf$sdraws2 - hflf$sdraws3
ds1s3 = hflf$sdraws1 - hflf$sdraws3


plot(density(ds1s2),col=4,xlim=c(-2.5,12),ylim=c(0,0.3),main="Posterior difference in 
             Cmax: dose = 15",xlab="Change in Cmax for dose 15",cex.main = 1.0)
legend("topright",legend=c("child-adoles.","adoles.-adult","child-adult"),col=c(4,3,2),
                                                           lty=1,lwd=2,bty="n",cex=0.8)
lines(density(ds2s3),col=3)
lines(density(ds1s3),col=2)
abline(v=0)


cbind(mean(ds1s2), sd(ds1s2))
quantile(ds1s2,probs=c(0.005,0.025,0.5,0.975,0.995))
cbind(mean(ds2s3), sd(ds2s3))
quantile(ds2s3,probs=c(0.005,0.025,0.5,0.975,0.995))
cbind(mean(ds1s3), sd(ds1s3))
quantile(ds1s3,probs=c(0.005,0.025,0.5,0.975,0.995))


# p-values and odds analytically - MUST ASSUME SAME SAMPLE SIZE FOR BOTH SAMPLES
# For differing sample sizes, do this numerically with MC draws
ts1s2 = mean(ds1s2)/sd(ds1s2)
ts2s3 = mean(ds2s3)/sd(ds2s3)
ts1s3 = mean(ds1s3)/sd(ds1s3)

rbind(cbind("       "," Child-Adoles."," Adoles.-Adult", "Child-Adult"),
cbind("p-value",round(1-pt(abs(ts1s2),df=11),3)*2, round(1-pt(abs(ts2s3),df=12),3)*2, round(1-pt(abs(ts1s3),df=12),3)*2),
cbind("odds vs. 0",round(dt(0,df=11)/dt(ts1s2,df=11),2),round(dt(0,df=12)/dt(ts2s3,df=12),2),
     round(dt(0,df=12)/dt(ts1s3,df=12),2)))

# 90% Credible intervals to help check p-values and odds
quantile(ds1s2,probs=c(0.05,0.95))
quantile(ds2s3,probs=c(0.05,0.95))
quantile(ds1s3,probs=c(0.05,0.95))


#### dose = 30
hflf = postdiff3(sm1=cmcm[3],ssd1=cmcsd[3],sm2=cmtm[3],ssd2=cmtsd[3],sm3=cmam[3],ssd3=cmasd[3],
              n1=12,n2=12,n3=14,R=200000)
mean(hflf$sdraws1); sd(hflf$sdraws1)
mean(hflf$sdraws2); sd(hflf$sdraws2)
mean(hflf$sdraws3); sd(hflf$sdraws3)

plot(density(hflf$sdraws1),col=4,xlim=c(9,30),ylim=c(0,0.26),main="Posterior Cmax: dose = 30", 
                                      xlab="Cmax for dose 30",cex.main = 1.0)
legend("topright",legend=c("Child","Adolescent","Adult"),col=c(4,3,2),lty=1,lwd=2,bty="n",cex=0.8)
lines(density(hflf$sdraws2),col=3)
lines(density(hflf$sdraws3),col=2)


print("summary stats for Child, Adolescent and Adult:")
mean(hflf$sdraws1); sd(hflf$sdraws1)
quantile(hflf$sdraws1,probs=c(0.005,0.025,0.5,0.975,0.995))
mean(hflf$sdraws2); sd(hflf$sdraws2)
quantile(hflf$sdraws2,probs=c(0.005,0.025,0.5,0.975,0.995))
mean(hflf$sdraws3); sd(hflf$sdraws3)
quantile(hflf$sdraws3,probs=c(0.005,0.025,0.5,0.975,0.995))


# differences
ds1s2 = hflf$sdraws1 - hflf$sdraws2
ds2s3 = hflf$sdraws2 - hflf$sdraws3
ds1s3 = hflf$sdraws1 - hflf$sdraws3


plot(density(ds1s2),col=4,xlim=c(-3.5,17),ylim=c(0,0.2),main="Posterior difference in 
             Cmax: dose = 30",xlab="Change in Cmax for dose 30",cex.main = 1.0)
legend("topright",legend=c("child-adoles.","adoles.-adult","child-adult"),col=c(4,3,2),
                                                           lty=1,lwd=2,bty="n",cex=0.8)
lines(density(ds2s3),col=3)
lines(density(ds1s3),col=2)
abline(v=0)


cbind(mean(ds1s2), sd(ds1s2))
quantile(ds1s2,probs=c(0.005,0.025,0.5,0.975,0.995))
cbind(mean(ds2s3), sd(ds2s3))
quantile(ds2s3,probs=c(0.005,0.025,0.5,0.975,0.995))
cbind(mean(ds1s3), sd(ds1s3))
quantile(ds1s3,probs=c(0.005,0.025,0.5,0.975,0.995))

# p-values and odds analytically - MUST ASSUME SAME SAMPLE SIZE FOR BOTH SAMPLES
# For differing sample sizes, do this numerically with MC draws
ts1s2 = mean(ds1s2)/sd(ds1s2)
ts2s3 = mean(ds2s3)/sd(ds2s3)
ts1s3 = mean(ds1s3)/sd(ds1s3)

rbind(cbind("       "," Child-Adoles."," Adoles.-Adult", "Child-Adult"),
cbind("p-value",round(1-pt(abs(ts1s2),df=11),3)*2, round(1-pt(abs(ts2s3),df=12),3)*2, round(1-pt(abs(ts1s3),df=12),3)*2),
cbind("odds vs. 0",round(dt(0,df=11)/dt(ts1s2,df=11),2),round(dt(0,df=12)/dt(ts2s3,df=12),2),
     round(dt(0,df=12)/dt(ts1s3,df=12),2)))

# 90% Credible intervals to help check p-values and odds
quantile(ds1s2,probs=c(0.05,0.95))
quantile(ds2s3,probs=c(0.05,0.95))
quantile(ds1s3,probs=c(0.05,0.95))





