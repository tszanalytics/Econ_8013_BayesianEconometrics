# to pause between plots:
par(ask=T)

# MCMC sample from other code (such as gibssegch3.R or Bayesregexample.R):
m <- 20   # size of small subsample to plot

b <- out$betadraw

b1 <- b[200:(201+m),1] 
b2 <- b[200:(201+m),2] 

plot(b1,type='l')
plot(b2,type='l')


##### start of gibbsgraph function
# create function for interactive - run from command prompt
gibbsgraph = function(b1,b2) {
library(MASS)
bb1 <- embed(b1,2)
bb2 <- embed(b2,2)

xlim <- c(min(b1),max(b1))
ylim <- c(min(b2),max(b2))

### plot(mean(b1),mean(b2),xlim=xlim,ylim=ylim,col=2,lwd=2)  

### try to get contour line
# estimate non-parameteric density surface via kernel smoothing
x <- b1; y <- b2
den<-kde2d(x, y, n=n)
# store z values of density estimate in an extra variable
den.z <-den$z

# this is the critical block, which I still do not comprehend in detail
z <- array()
for (i in 1:n){
        z.x <- max(which(den$x < x[i]))
        z.y <- max(which(den$y < y[i]))
        z[i] <- den$z[z.x, z.y]
}

# store class/level borders of confidence interval in variables
confidence.border <- quantile(z, probs=c(0.1,0.5,0.8), na.rm = TRUE) # 

plot(x,y,col="grey")
par(new=TRUE)
contour(den, levels=confidence.border, col = c(2,3,5), labcex=0.8, lwd=2, add = TRUE,labels=c("0.1","0.5","0.8"),drawlabels=T)

# contour(den, levels=confidence.border, col = 2, add = TRUE)


for (i in 1:20) {

x <- c(bb1[i,2],bb1[i,1])
y <- c(bb2[i,1],bb2[(i+1),2])
lines(x,y,type='l',col="blue",lwd=2)
line <- readline()

x <- c(bb1[i,1],bb1[(i+1),2])
y <- c(bb2[(i+1),2],bb2[(i+1),1])
lines(x,y,type='l',col="blue",lwd=2)
line <- readline()
}

}
########### end of function

#*** now call the function in R console: gibbsgraph(b1,b2) ***




###############################################
### How to do contour plots example below
set.seed(1388)
library(MASS)
n <- 2000
x <- rnorm(n); y <- 0.5*x + 0.5*rnorm(n)
plot(x,y,col="grey")

den<-kde2d(x, y, n=n)
# store z values of density estimate in an extra variable
den.z <-den$z

# this is the critical block, which I still do not comprehend in detail
z <- array()
for (i in 1:n){
        z.x <- max(which(den$x < x[i]))
        z.y <- max(which(den$y < y[i]))
        z[i] <- den$z[z.x, z.y]
}

# store class/level borders of confidence interval in variables
confidence.border <- quantile(z, probs=c(0.1,0.5,0.8), na.rm = TRUE) 

par(new=TRUE)
contour(den, levels=confidence.border, col = c(2,3,4), labcex=0.8, lwd=2, add = TRUE,labels=c("0.1","0.5","0.8"),drawlabels=T)






### try to get contour line
# estimate non-parameteric density surface via kernel smoothing
x <- b1; y <- b2
den<-kde2d(x, y, n=n)
# store z values of density estimate in an extra variable
den.z <-den$z

# this is the critical block, which I still do not comprehend in detail
z <- array()
for (i in 1:n){
        z.x <- max(which(den$x < x[i]))
        z.y <- max(which(den$y < y[i]))
        z[i] <- den$z[z.x, z.y]
}

# store class/level borders of confidence interval in variables
confidence.border <- quantile(z, probs=1-0.9, na.rm = TRUE) # +-1sd

plot(x,y)
par(new=TRUE)
contour(den, levels=confidence.border, col = "red", add = TRUE)



# try following
library(MASS)

# parameters:
n<-100

# generate samples:
set.seed(138813)
#seed <- .Random.seed
x<-rnorm(n); y<-rnorm(n)

# estimate non-parameteric density surface via kernel smoothing
den<-kde2d(x, y, n=n)
# store z values of density estimate in an extra variable
den.z <-den$z

# this is the critical block, which I still do not comprehend in detail
z <- array()
for (i in 1:n){
        z.x <- max(which(den$x < x[i]))
        z.y <- max(which(den$y < y[i]))
        z[i] <- den$z[z.x, z.y]
}

# store class/level borders of confidence interval in variables
confidence.border <- quantile(z, probs=1-0.6827, na.rm = TRUE) # +-1sd

plot(x,y)
par(new=TRUE)
contour(den, levels=confidence.border, col = "red", add = TRUE)



# library(plot3D)
xy <- cbind(b[,1],b[,2])

    nbins <- 100
    x.bin <- seq(floor(min(xy[,1])), ceiling(max(xy[,1])), length=nbins)
    y.bin <- seq(floor(min(xy[,2])), ceiling(max(xy[,2])), length=nbins)

    freq <-  as.data.frame(table(findInterval(xy[,1], x.bin),findInterval(xy[,2], y.bin)))
    freq[,1] <- as.numeric(freq[,1])
    freq[,2] <- as.numeric(freq[,2])

    freq2D <- diag(nbins)*0
    freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]

    # par(mfrow=c(1,2))
    image(x.bin, y.bin, freq2D) #, col=topo.colors(max(freq2D)))
 

plot(0.5,0.0,xlim=c(0,1),ylim=c(-0.1,0.1),col=2,lwd=2)
   contour(x.bin, y.bin, freq2D, add=TRUE, col=3)


require(grDevices) # for colours
x <- -6:16
op <- par(mfrow = c(2, 2))
contour(outer(x, x), method = "edge", vfont = c("sans serif", "plain"))
z <- outer(x, sqrt(abs(x)), FUN = "/")
image(x, x, z)
contour(x, x, z, col = "pink", add = TRUE, method = "edge",
        vfont = c("sans serif", "plain"))
contour(x, x, z, ylim = c(1, 6), method = "simple", labcex = 1)
contour(x, x, z, ylim = c(-6, 6), nlev = 20, lty = 2, method = "simple")
par(op)



