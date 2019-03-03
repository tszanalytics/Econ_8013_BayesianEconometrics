mcprobit <- function(X,y,b0,nit){
  N <- length(y)
  k <- length(b0)
  B <- matrix(0,nit,k)
  llh <- rep(0,N)
  if(any((y!=0)&(y!=1))){error("dependent var must be 0's and 1's")}
  for( id in 1:nit){
    z <- runif(N)
    yhat <- X %*% b0
    thresh <- pnorm(yhat)
    z <- ((y==1)*thresh*z+(y==0)*(thresh+(1-thresh)*z))
    z <- -qnorm(z)+yhat
    est <- lm.fit(X,z)
    b0 <- est$coefficients
    R <- qr.R(est$qr)
    b0 <- b0+solve(R,rnorm(k))
    B[id,] <- b0
    llh[id] <- sum(y*log(pnorm(yhat))+(1-y)*log(1-pnorm(yhat)))
  }
  return(list(B=B,llh=llh))
}
