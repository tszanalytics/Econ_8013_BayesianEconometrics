acregalt = function(y,x,R,rmax=F){
#
#
# purpose: compute AR(1) errors linear regression using default
# relatively uninformative normal-gamma prior
#
# usage: out = acreg(y,x,R)
# This creates data list out with coeff estimates, etc., 
# type str(out) to see its structure
#
# arguments:    	y = vector of dep var
# 			x = array of indep vars
#                 R = number of MCMC iterations >= 2000
#              rmax = T for |rho| <1, rho unrestricted otherwise (default)   

# output:
# MCMC sample for regression coeffs, equation variance and the AR(1) coefficient, rho

x = cbind(1,x)
k = ncol(x)
n = nrow(x)
In = diag(1,n)


# prior
b0        = rep(0,k)
V0        = diag(100,k)
nu0       = 0.005
s02       = 1.0
rho0      = 0.5
Vrho      = 1000
iV0       = solve(V0)
iV0b0     = iV0%*%b0
nu0n      = nu0+n
nu0s02    = nu0*s02
iVrho     = 1/Vrho
iVrhorho0 = iVrho*rho0


# Initial values
# --------------
ols   = lm(y~x-1)
b     = ols$coef
s2    = mean(ols$res^2)
ols   = lm(ols$res[2:n]~ols$res[1:(n-1)]-1)
r     = as.numeric(ols$coef)
s2    = as.numeric(mean(ols$res^2))
pars = c(b,s2,r)
O      = matrix(0,n,n)
for (iter in 1:5){
  for (i in 1:n)
    O[i,] = r^abs(((1-i):(n-i)))
  O     = O/(1-r^2)
  iO    = solve(O)
  xiOx  = t(x)%*%iO%*%x
  xiOy  = t(x)%*%iO%*%y
  b     = solve(xiOx)%*%xiOy
  e     = y - x%*%b
  s2    = mean(e^2)*(1-r^2)
  r     = lm(e[2:n]~e[1:(n-1)]-1)$coef
  pars  = rbind(pars,c(b,s2,r))
}
b = pars[5,1:k]
s2 = pars[5,3]
r  = pars[5,4]

# MCMC scheme
set.seed(1268)
burn  = 1000
M     = R
LAG   = 1
niter = burn+LAG*M
s2s   = rep(0,niter)
rs    = rep(0,niter)
bs    = matrix(0,niter,k)
for (iter in 1:niter){
#  print(iter)
  yr <- y[2:n] - r*y[1:(n-1)]
  xr <- x[2:n,] - r*x[1:(n-1),]
  er <- yr-xr%*%b
  s2 = 1/rgamma(1,(nu0n)/2,(nu0s02+t(er)%*%er)/2)  
  V1 = solve(iV0+t(xr)%*%xr/s2)
  b1 = V1%*%(iV0b0+t(xr)%*%yr/s2)
  b  = b1+t(chol(V1))%*%rnorm(k)
  e  = y-x%*%b
  e1 = e[2:n]
  e2 = e[1:(n-1)]
  Vr = 1/(iVrho+sum(e2^2)/s2)
  r1 = Vr*(iVrhorho0+sum(e1*e2)/s2)


if(rmax == T) {
# restricting |rho|<1
  A  = pnorm(-1,r1,sqrt(Vr))
  B  = pnorm(1,r1,sqrt(Vr))
  u  = runif(1)
  r  = r1+sqrt(Vr)*qnorm(u*B+(1-u)*A)
} else {
r  = r1+sqrt(Vr)*rnorm(1)   # without restricting rho
}
  
  s2s[iter]   = s2
  bs[iter,]   = t(b)
  rs[iter]    = r
}
s2s    = s2s[seq(burn+1,niter,by=LAG)]
bs     = bs[seq(burn+1,niter,by=LAG),1:k]
rs     = rs[seq(burn+1,niter,by=LAG)]


list(bs=bs,s2s=s2s,rs=rs)

}




