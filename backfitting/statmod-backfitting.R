# StatMod2 - Backfitting HW - Exercices 4

library(stats)

setwd("~/Google Drive/2. SPRING 2015/STAT MOD 2 - Prof Scott/backfitting")
data <- read.csv("air.csv")
data <- data[order(data[, 1]), ]

y <- as.matrix(data$Ozone)
x <- data[,c("Solar.R","Wind","Temp")]
n <- length(data[,1])
p <- dim(x)[2]
times <- 20

# Initialize alpha, f, and mses.
a <- mean(y)
f <- array(0, dim=c(p, n))
mses <- NULL

PredictOutput <- function(pairs, x.star, bandwidth) {
  output <- 0
  weight.sum <- 0
  len <- length(pairs[,2])
  for (i in 1:len) {
    x.i <- pairs[i,1]
    y.i <- pairs[i,2]
    weight.i <- dnorm((x.i-x.star)/bandwidth)/bandwidth
    weight.sum <- weight.sum + weight.i
    output.i <- weight.i*y.i
    output <- output + output.i
  }
  # Normalize weights by dividing by weight sum.
  output <- output/weight.sum
  return (output)
}

Backfit <- function(q) {
  mses <- c(1000)
  delta <- 100
  t <- 1
  
  # Do t iterations, until convergence criterion.
  while (delta > 0.05) {
  
    # Do for each partial residual (kth predictor).
    #for (k in c(1, 2, 3)) {
    #for (k in c(3, 2, 1)) {
    for (k in c(2, 1, 3)) {
      
      # Remove one column, to make p-1 by 1 matrix.
      f.minus.k <- as.matrix(f[-k,])
      
      # Prepare value to smooth; the Yi - a - sum(...).
      to.smooth <- rep(0, n)
      for (i in 1:n) {
        to.smooth[i] <- y[i] - a - sum(f.minus.k[,i])
      }
      
      # Fit kth residual against kth column of x.
      partial.data <- cbind(y, x[k], to.smooth)
      partial.data <- partial.data[order(partial.data[,2]),]
      predicted <- rep(0, n)
      for (i in 1:n) {
        x.star <- partial.data[i,2]
        # Do "2%" smoothing, by making bandwidth 1/50 of x range.
        b <- (max(x[k])-min(x[k]))/50
        predicted[i] <- PredictOutput(partial.data[,c(2,3)], x.star, b)
      }
      #predicted <- lowess(partial.data[,c(2,3)], f=.1)$y
      
      # Plot stuff.
      if (k==q) {
        if (t==1) {
          plot(t(x[k]), to.smooth, main=bquote("Fitting Residuals on X"[.(k)]),
               xlab=bquote("X"[.(k)]~" Value"), ylab="Residual Value")
        }
        # Plot lines, ignoring first 10 (burn-in).
        if (t>0) {
          lines(partial.data[,2], predicted, col=t)
        }
      }
      
      # Attach predicted values to a resorted (by y) matrix.
      to.sort <- cbind(partial.data, predicted)
      resorted <- to.sort[order(to.sort[,1]),]
      # Record K-th predictor for t-th iteration.
      f[k,] <- resorted$predicted - mean(resorted$predicted)
  
    }

    # Record MSE for this run, and compute delta to test convergence.
    mses[t+1] <- mean((y - a - colSums(f))^2)
    delta <- abs(mses[t+1]-mses[t])
    print(bquote("MSE delta=abs("~.(mses[t+1])~"-"~.(mses[t])~")="~.(delta)))
    t <- t+1
  }
}

par(mfrow=c(3,1))
for (q in 1:p) {
  mses <- Backfit(q)
}

print (length(mses))
print (summary(mses))


