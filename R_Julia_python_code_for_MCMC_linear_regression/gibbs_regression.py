# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 17:43:29 2016

@author: millsjf
"""

# Regression MCMC (Gibbs)
import numpy as np
import matplotlib.pyplot as plt

# dgp
n = 100
b = [1,1]
k = len(b)
mu = 0
sigma = 2
e = np.random.normal(mu, sigma, 100) 
x = np.random.normal(mu, sigma, 100) 
# y = [b[0] + b[1]*x[i] + e[i] for i in range(n)]

y = b[0] + b[1]*x + e

# plot the data
plt.plot(x, y, 'ro')
plt.show()

# OLS regression estimation
import statsmodels.api as sm
from patsy import dmatrices

# automatically includes an intercept this way:
y, X = dmatrices('y ~ x', return_type='dataframe')
reg = sm.OLS(y,X).fit()
print(reg.params)
print(reg.bic)
print(reg.summary())

####  TO GET ALL POSSIBLE EXTRACTIONS
# dir(reg)

# Gibbs0
# coeff. standard errors from regression
print(reg.bse)
s2 = reg.ssr/(n-len(b))  # estimated sigma^2

# prior
b0        = np.zeros(k)
b0.shape = (k,1)
V0        = 100*np.identity(k)
iV0       = np.linalg.inv(V0)
iV0b0     = iV0@b0
b0iV0b0   = b0.T@iV0b0
d0 = 0.01333
nu0 = 5


# starting values
s2draw = np.array([s2])
XpXi = np.linalg.inv(np.dot(X.T, X))
xpxdiag = np.sqrt(np.diag(XpXi))   # for uninformative - uniform prior
nu = nu0 + n

d0yy = d0 + np.dot(y.T,y) + b0iV0b0

iV1 = iV0 + np.dot(X.T,X)
V1 = np.linalg.inv(iV1)  # for 'informative' conjugate prior
v1diag = np.sqrt(np.diag(V1))
b1 = V1@(iV0b0 + np.dot(X.T,y))






iter = 5000
burnin = 500
M = iter + burnin
### which is faster .append, or fill in using index?
##  seems that concatenate may be fastest:
# s2draw = np.array([0])
# s2draw = np.concatenate((s2draw, s20))


# bdraw = np.array([])
## bdraw.append() = [1,2]
bdraw = np.zeros((M,k))
# bdraw[:,0]

# mcmc loop
d1 = d0yy - b1.T@iV1@b1
for i in range(1,M):
    # compute new variances for bdraws
    # use xpxdiag and reg.params for uninformative prior
    # use v1diag and b1 for conjugate prior (default is uninformative)
    seb = np.sqrt(s2draw[i-1])*v1diag
    # draw betas
    bdrawi = np.random.normal(b1.T, seb.T, k)
    bdrawi.shape = (k,1)
    
    # compute new shape and scale given bdraw
    # res = y - np.dot(X,bdrawi)
    # rrs = np.dot(res.T,res)  
    # draw sigma2
    d2 = d1
    d1 = d0yy - bdrawi.T@iV1@bdrawi
    if(d1<=0):      ### occasionally get d1 <0 ?!
        d1 = d2
    tau = np.random.gamma(shape=nu/2, scale=2/d1)
    s2drawi = 1/tau
    # R code:
#d1 = d0 + t(y)%*%y + b0iV0b0 - t(b1)%*%iV1%*%b1
#    tau = rgamma(1,shape = nu/2,rate = d1/2)
#    s2 = 1/tau
#    bs[iter,]   = t(b)
#    taus[iter]   = tau
    s2drawi = np.array([s2drawi])
# s2draw = np.concatenate((s2draw, s20))
    s2draw = np.concatenate((s2draw, s2drawi))
    # bdraw = np.concatenate((bdraw, bdrawi))
    bdraw[i,:] = bdrawi.T

print("MCMC estimate for b0 = ", np.mean(bdraw[burnin:M,0]))
print("MCMC estimate for b1 = ", np.mean(bdraw[burnin:M,1]))
print("MCMC estimate for s2 = ",np.mean(s2draw[burnin:M]))
print()
print("OLS b estimates = ",reg.params)
print("OLS s2 estimate = ",s2)

# plot MCMC chain
import matplotlib.pyplot as plt
# import scipy.special as sps

# plot of chain
plt.plot(bdraw[burnin:M,0])

# histogram of density
plt.hist(bdraw[burnin:M,1], 50, normed=True)

mu = np.mean(bdraw[burnin:M,1])
sigma = np.std(bdraw[burnin:M,1])
count, bins, ignored = plt.hist(bdraw[burnin:M,1], 50, normed=True)
plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu)**2 / (2 * sigma**2) ), linewidth=2, color='r')
plt.show()




# must define scale and shape parameters for the posterior of sigma
count, bins, ignored = plt.hist(s2draw[burnin:M], 50, normed=True)
# y = bins**(shape-1)*(np.exp(-bins/scale) / (sps.gamma(shape)*scale**shape))
# plt.plot(bins, y, linewidth=2, color='r')
# plt.show()
