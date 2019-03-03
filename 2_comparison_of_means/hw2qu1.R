# homework 2 qu.1

RDAdata = read.csv(file="Hahn_Table54p107.csv")

group1 = RDAdata$group1
group2 = RDAdata$group2

x = group1
y = group2

## function to compute posterior draws
source(file="posterior_NIGprior.R")
M = 10^6
mux_draws = posterior_NIGprior(x, M=M)
muy_draws = posterior_NIGprior(y, M=M)

plot(density(muy_draws),col=4,xlim=c(0.0,0.6))
lines(density(mux_draws),col=3)             #,ylim=c(0.0,20.0))
abline(v=0.0,lwd=2)
#abline(v=mean(x))
# posterior summary statistics
mean(mux_draws)
sd(mux_draws)
# 0.99, 095 HPD intervals and median:
quantile(mux_draws, probs=c(0.005,0.025,0.5,0.975,0.995))

diff = mux_draws - muy_draws
lines(density(diff),col=2)

# summary stats for difference
mean(diff)
sd(diff)
# 0.99, 095 HPD intervals and median:
quantile(diff, probs=c(0.005,0.025,0.5,0.975,0.995))

# posterior odds
source(file="postoddsmc.R")
postoddsmc(diff)

# Bayesian p-value = area in tail (pick the appropriate tail!)
d0 = ifelse(diff< 0, 1, 0)
pval = sum(d0)/length(diff)
pval



########### Examination of priors (used in class) below:

# same prior for each group
# mu ~ Normal(0, 10^6)
# tau = (1/s^2) ~ Gamma(0.001, 0.002)

priormu = rnorm(mean=0, sd = sqrt(10^6),100000)
plot(density(priormu))

par(mfrow=c(1,2))
priortau = rgamma(1000000, shape = 3, scale = 0.001)
plot(density(priortau),main="tau")

priorsig2 = 1/priortau
plot(density(priorsig2), main = "sigma sq.")

max(priorsig2)
summary(priortau)
min(priortau)

t1 = rgamma(10000,0.001, 0.001)
plot(density(t1))



