#hahnp37.R

# Hahn (2014) p.37-38 difference in proportions example

n.iter <- 10000

pi.w <- rbeta(n.iter,15,24)
pi.m <- rbeta(n.iter,4,20)
diff <- pi.w - pi.m

hist(diff,50)


