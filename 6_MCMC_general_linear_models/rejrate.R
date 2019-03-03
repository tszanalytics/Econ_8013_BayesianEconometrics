rejrate = function(x) {
# rejrate.R

# x = MCMC sample
# calculate MCMC chain rejection rate
# output = rejection rate (# of sample values repeated in the chain)

m <- length(x)
rej <- 0
for (i in 2:m) {
d <- ifelse(x[i] == x[i-1], 1, 0)
rej <- rej + d
}
rrate <- rej/m
rrate
}