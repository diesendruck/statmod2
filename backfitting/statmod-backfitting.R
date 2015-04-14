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
  # Predicts probability density for given x.star, based on pairs of data and
  # a weight function. Specifically, doing sum(weight.i*y.i), where each weight
  # comes from WeightFunction.
  #
  # Args:
  #   pairs: The pairs of X's and corresponding Y's, i.e. the data.
  #   x.star: This is the point of interest, whose probability under the
  #     density you want to find.
  #   bandwidth: A neighborhood around the requested x.star, used to average 
  #     the area under the smoother function.
  #
  # Returns:
  #   output: The predicted value for the function of interest, at x.star.
  output <- 0
  weight.sum <- 0
  len <- length(pairs[,2])
  for (i in 1:len) {
    x.i <- pairs[i,1]
    y.i <- pairs[i,2]
    weight.i <- WeightFunction(x.i, x.star, bandwidth)
    weight.sum <- weight.sum + weight.i
    output.i <- weight.i*y.i
    output <- output + output.i
  }
  # Normalize weights by dividing by weight sum.
  output <- output/weight.sum
  return (output)
}

WeightFunction <- function(x.i, x.star, bandwidth) {
  # Gives weight of the Y associated with x.i, according to a Gaussian kernel.
  #
  # Args:
  #   x.i: Nearby X value to x.star, used to determine a weight for its Y.
  #   x.star: The center of the smoother function.
  #   bandwidth: A neighborhood around the requested x.star, used to average 
  #     the area under the smoother function.
  #
  # Returns:
  #   weight: The weight of the y that corresponds to x.i.
  density <- ((2*pi)^(-1/2))*exp(-1/2*((x.i-x.star)/bandwidth)^2)
  weight <- density/bandwidth
  return (weight)
}

# Do t iterations.
for (t in 1:times) {
  # Rather than convergence, I record all MSEs for review, below.
  
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
    partial.data <- cbind(y, x[k],to.smooth)
    partial.data <- partial.data[order(partial.data[,2]),]
    predicted <- rep(0, n)
#     for (i in 1:n) {
#       x.star <- partial.data[i,1]
#       predicted[i] <- PredictOutput(partial.data, x.star, 10)
#     }
    predicted <- lowess(partial.data[,c(2,3)], f=.1)$y
   
    # Plot stuff.
    if (k==2) {
      if (t==1) {
        plot(t(x[k]), to.smooth, main=k, ylim=c(-50,100))
      }
      lines(partial.data[,2], predicted, col=t)
    }
    
    # Attach predicted values to a resorted (by y) matrix.
    to.sort <- cbind(partial.data, predicted)
    resorted <- to.sort[order(to.sort[,1]),]
    # Record K-th predictor for t-th iteration.
    f[k,] <- resorted$predicted - mean(resorted$predicted)
#     if (k==3) {
#       print(f[k,1:5])
#     }
  }

  mses[t] <- (1/n)*sum((y - a - colSums(f))^2)
  print (mses[t])
}

print (summary(mses))


