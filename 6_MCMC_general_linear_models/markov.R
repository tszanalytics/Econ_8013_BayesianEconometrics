# markov.R
# Markov chain example
P  <- matrix(c(0.75,0.25,0.125,0.875), nrow = 2, ncol=2, byrow=TRUE)
n = 20 # calculate P^n
Z = P
for (i in 1:(n-1)) {
Z = Z%*%P
cat("iteration",i,Z,"\n")
}
Z