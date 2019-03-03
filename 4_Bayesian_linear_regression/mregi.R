mregi = function(y,x,b0=0,Vb0=10^100,lambda=0.1333,nu=3){
#
# purpose: compute standard linear regression using 
# relatively informative normal-gamma prior
#
# usage: out = mreg(y,x)
# This creates data list out with coeff estimates, etc., 
# type str(out) to see its structure
#
# arguments:    	y = vector of dep var
# 			x = array of indep vars

# output:
# list containing marginal posterior means, s.e's, 95% CIs 
# summary stats, AIC, BIC, odds test of zero coeffs
# plots of marginal posteriors


X = cbind(1,x)
k = ncol(X)
n = nrow(X)
In = diag(1,n)

# Prior hyperparameters
if(length(b0)==1) b0 = rep(b0,k)
if(length(Vb0)==1) Vb0 = diag(Vb0,k)
# diffuse prior:
# b0 = rep(0,k)
# Vb0   = diag(10^100,k)
# lambda = 0.1333
# nu     = 3

# Marginal posterior distributions of sigma2 and beta
d0   = (nu+n)/2

# posterior means and SDs
XXI    = solve(t(X)%*%X+solve(Vb0))
b1   = XXI%*%(t(X)%*%y+solve(Vb0)%*%b0)
d1  = (nu*lambda+t(y)%*%y+t(b0)%*%solve(Vb0)%*%b0-t(b1)%*%solve(XXI)%*%b1)/2
s21   = d1/(d0+1)
sds2  = s21/sqrt(d0+2)
nu1lambda1=nu*lambda+t(y-X%*%b0)%*%(In-X%*%solve(solve(Vb0)+t(X)%*%X)%*%t(X))%*%(y-X%*%b0)
sdb=sqrt((as.numeric(nu1lambda1)/(nu+n))*diag(XXI))

ym = mean(y)
tss = t(y-ym)%*%(y-ym)
yp = X%*%b1
res = y - yp
rss = t(res)%*%res
rsq = 1-rss/tss
aic = log(rss/n) + 2*(k+1)/n
bic = log(rss/n) + log(n)*(k+1)/n
df = n-k

matout = round(rbind(cbind(b1,sdb),c(s21,sds2)),5)
cat(" posterior Mean (col 1) & Std Dev (col 2) for beta (rows 1 to k) & sigma (row k+1)",fill=TRUE)
print(matout)

cat(" R-squared, AIC, BIC")
print(cbind(rsq,aic,bic))

list(resid=res,yhat=yp,rss=rss,rsquare=rsq,aic=aic,bic=bic,df=df,
			b=b1,bse=sdb,s2=s21,ses2=sds2)

}
