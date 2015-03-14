# Stat Mod 2 - Exercises 2 - Smoothing (p.3)
# 
# Performs curve fitting and calculates MSE of estimated density functions.
# Original data function and the kernel smoother are hard-coded in functions
# "GeneratePairs", "WeightFunction", and "ComputeMSE".
#
# Author: Maurice Diesendruck
# Last updated: 2015-02-17

Main <- function() {
  # Runs all operations.
  #
  # Args:
  #   NA
  #
  # Returns:
  #   NA
  fitted.curves <- FitCurves(1000, c(0.1, 1, 4))
  mses <- ComputeMSE(fitted.curves)
  mses
}

FitCurves <- function(sample.size, bandwidths) {
  # Wraps curve fitting steps into one, and allows choice of three bandwidths.
  #
  # Args:
  #   sample.size: Scalar number that determines size of sample.
  #   bandwidths: Vector of three bandwidths.\
  #
  # Returns:
  #   fitted.curves: List of X-values, three sets of predicted Y's based on
  #     different bandwidths, and means for original Xs and Ys.
  # Generate new pairs of original data.
  data <- GeneratePairs(sample.size)
  pairs <- data[[1]]
  plot(pairs[,2], pairs[,1], col="gray", xlab="Sample X-Values",
       ylab="Sample Y-Values", main="Curve Fitting with Normal Kernel
       Smoothing and Varying Bandwidths")
  
  # Predict output and plot.
  x <- seq(min(pairs[,2]), max(pairs[,2]), .1)
  par(bty="l",las=1)
  output1 <- PredictOutput(pairs, x, bandwidths[1])
  output2 <- PredictOutput(pairs, x, bandwidths[2])
  output3 <- PredictOutput(pairs, x, bandwidths[3])
  lines(x, output1, col="steelblue1", lwd="2")
  lines(x, output2, col="dodgerblue3", lwd="2")
  lines(x, output3, col="grey20", lwd="2")
  legend("bottomright", title="Bandwidth", legend=bandwidths, bty="n",
         lty=c(1,1),lwd=c(3,3), col=c("steelblue1","dodgerblue3","grey20"))
  
  fitted.curves <- list(x, output1, output2, output3, data[[2]])
  return (fitted.curves)
}

GeneratePairs <- function(n) {
  # Takes set of non-linear X's and applies the y = x^2 + e transformation,
  # where error values "e" are Normal(0,1). Then subracts means from each.
  #
  # Args:
  #   n: Number of pairs to be generated.
  #
  # Returns:
  #   data: List, with mean-subtracted pairs of inputs X and outputs Y, and
  #     with vector of means.
  x.input <- rnorm(n, 5, 2)
  errors <- rnorm(n, 0, 1)
  y.output <- x.input^2 + errors
  
  # Subract out means.
  ms.input <- x.input - mean(x.input)
  ms.output <- y.output - mean(y.output)
  
  ms.pairs <- cbind(ms.output, ms.input)
  pairs <- cbind(y.output, x.input)
  data <- list(ms.pairs, c(mean(x.input), mean(y.output)), pairs)
  return (data)
}

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
    x.i <- pairs[i,2]
    y.i <- pairs[i,1]
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

ComputeMSE <- function(fitted.curves) {
  # Compares the predictions (based on the training set), to the true function
  # values (from the relationship y=x^2). Predictions are mean-subtracted, so
  # true function values (both X and Y) must be adjusted before comparison.
  # 
  # Args:
  #   fitted.curves: List of X-values, and three sets of predicted Y's based on
  #     different bandwidths.
  #
  # Returns:
  #   mses: Vector of mean squared error values for each of three functions.
  x.values <- fitted.curves[[1]]
  y.predictions <- list(fitted.curves[[2]],
                        fitted.curves[[3]],
                        fitted.curves[[4]])
  x.mean <- fitted.curves[[5]][1]
  y.mean <- fitted.curves[[5]][2]
  num.xvalues <- length(x.values)
  num.ypredictions <- length(y.predictions)
  mses <- c()
  
  # Compute sum of squared errors for all y prediction functions.
  for (i in 1:num.ypredictions) {
    sum.sq.err.i <- 0
    pred.func.i <- y.predictions[[i]]
    
    # For given function, compute squared errors, and add to sum.
    for (j in 1:num.xvalues) {
      # Mean of Y is added back to predicted, mean-subtracted Y outputs.
      predicted.value <- pred.func.i[j]+y.mean
      # Mean of X is added back to original mean-subtracted X inputs.
      expected.value <- (x.values[j]+x.mean)^2
      sq.err <- (predicted.value - expected.value)^2
      sum.sq.err.i <- sum.sq.err.i + sq.err
    }
    
    mse.element <- (1/num.xvalues)*sum.sq.err.i
    mses <- c(mses, mse.element)
  }
  
  return (mses)
}

Main()



