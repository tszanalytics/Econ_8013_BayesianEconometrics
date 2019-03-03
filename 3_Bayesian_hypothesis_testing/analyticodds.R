# analyticodds.R

# evaluate odds for Normal comparing 0 and max
max <- 1.0  # put value of max = mean here
null <- 0.0
numN <- dnorm(max,mean=max,sd=1)  # set mean and sd
denN <- dnorm(null, mean=max,sd=1)

poddsN <- numN/den

# Same for a t density with 8 df
numt <- dt(max,df=8)
dent <- dt(null,df=8)

poddst <- numt/dent
