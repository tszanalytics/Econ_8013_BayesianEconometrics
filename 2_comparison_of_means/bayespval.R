# bayespval.R

bayespval = function(mcsample,nullh=0.0) {

# numerical Bayesian p-value computation from MC sample
# use: set bayespval(mcsample,testPoint=nullvalue)
#      mcsample = sample from the posterior
#      nullh = value under the null hypothesis (to be tested)

# purpose: evaluate the p-value from the posterior against a point null hypothesis
#          give a sample from the posterior density (from MC or MCMC)

# output:  Bayesian posterior p-value
#          from an MC sample from the posterior density

# Jeff Mills, 2017

DenB1 <- density(mcsample,n=1024)  #### posterior density
# find null point
index <- which(DenB1$x-nullh==min(abs(DenB1$x-nullh)))
if (length(index) == 0) {
  index <- which(-DenB1$x+nullh==min(abs(DenB1$x-nullh)))
}
mm = length(DenB1$x)
pval = 2*sum(DenB1$y[index:mm])/sum(DenB1$y)
if(pval >= 1) {
	pval = 2*(1-sum(DenB1$y[index:mm])/sum(DenB1$y))
}

return(list(pvalue = pval))
}

