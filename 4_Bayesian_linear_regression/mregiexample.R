# mregiexample.R

# example of use of mregi 

# some fake data:
x <- rnorm(100)
y <- 1.0 + 2.0*x + 2.0*rnorm(100)
summary(lm(y~x))

source(file="mregi.R")

# default uninformative prior
outu <- mregi(y,x)

# informative prior
b0 <- c(4.0,5.0)
Vb0 <- diag(0.01,2)
out <- mregi(y,x,b0=b0,Vb0=Vb0)

