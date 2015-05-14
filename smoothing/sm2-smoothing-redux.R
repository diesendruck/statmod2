# Stat Mod 2 - Exercises 2 - Smoothing (p.4)
# 
# Performs curve fitting and calculates MSE of estimated density functions.
# Original data function and the kernel smoother are hard-coded in functions
# "GeneratePairs", "WeightFunction", and "ComputeMSE".
#
# Author: Maurice Diesendruck
# Last updated: 2015-02-17
library(ggplot2)
library(reshape)

Main <- function() {
  # Performs curve fitting and smoothing for four different training data sets,
  # each with three different bandwidths; and plots their mean squared errors.
  #
  # Args:
  #   NA
  #
  # Returns:
  #   NA
  fitted.curves1 <- FitCurves(1000, c(0.01, 0.1, 2), GenerateSmoothClean)
  mses1 <- ComputeMSE(fitted.curves1, flag="Smooth")
  
  fitted.curves2 <- FitCurves(1000, c(0.01, 0.1, 2), GenerateSmoothNoisy)
  mses2 <- ComputeMSE(fitted.curves2, flag="Smooth")
  
  fitted.curves3 <- FitCurves(1000, c(0.01, 0.1, 2), GenerateWigglyClean)
  mses3 <- ComputeMSE(fitted.curves3, flag="Wiggly")
  
  fitted.curves4 <- FitCurves(1000, c(0.01, 0.1, 2), GenerateWigglyNoisy)
  mses4 <- ComputeMSE(fitted.curves4, flag="Wiggly")
  
  # Plot results.
  d1 <- cbind(rep("Smooth Clean",3), c("0.01","0.1","2"), mses1)
  d2 <- cbind(rep("Smooth Noisy",3), c("0.01","0.1","2"), mses2)
  d3 <- cbind(rep("Wiggly Clean",3), c("0.01","0.1","2"), mses3)
  d4 <- cbind(rep("Wiggly Noisy",3), c("0.01","0.1","2"), mses4)
  df1 <- data.frame(rbind(d1, d2, d3, d4), stringsAsFactors=FALSE)
  names(df1) <- c("Case","Bandwidth","Value")
  df1[,3] <- sapply(df1[,3], as.numeric)
  ggplot(data=df1, aes(x=Case, y=Value, fill=Bandwidth, order=Bandwidth)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle("MSE Values for Smooth, Wiggly, Clean, and Noisy Data") + 
    theme(plot.title = element_text(lineheight=.8, face="bold"))
}

GenerateSmoothClean <- function(n) {
  # Takes set of non-linear X's, applies a function, and adds noise.
  #
  # Args:
  #   n: Number of pairs to be generated.
  #
  # Returns:
  #   data: List, with mean-subtracted pairs of inputs X and outputs Y, and
  #     with vector of means.
  x.input <- seq(0, 1, length=100)
  errors <- rnorm(n, 0, 0.01)
  y.output <- x.input^2 + errors
  
  # Subract out means.
  ms.input <- x.input - mean(x.input)
  ms.output <- y.output - mean(y.output)
  
  ms.pairs <- cbind(ms.output, ms.input)
  pairs <- cbind(y.output, x.input)
  data <- list(ms.pairs, c(mean(x.input), mean(y.output)), pairs)
  return (data)
}
GenerateSmoothNoisy <- function(n) {
  # Takes set of non-linear X's, applies a function, and adds noise.
  #
  # Args:
  #   n: Number of pairs to be generated.
  #
  # Returns:
  #   data: List, with mean-subtracted pairs of inputs X and outputs Y, and
  #     with vector of means.
  x.input <- seq(0, 1, length=100)
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
GenerateWigglyClean <- function(n) {
  # Takes set of non-linear X's, applies a function, and adds noise.
  #
  # Args:
  #   n: Number of pairs to be generated.
  #
  # Returns:
  #   data: List, with mean-subtracted pairs of inputs X and outputs Y, and
  #     with vector of means.
  x.input <- seq(0, 1, length=100)
  errors <- rnorm(n, 0, 0.01)
  y.output <- sin(2*pi*x.input) + errors
  
  # Subract out means.
  ms.input <- x.input - mean(x.input)
  ms.output <- y.output - mean(y.output)
  
  ms.pairs <- cbind(ms.output, ms.input)
  pairs <- cbind(y.output, x.input)
  data <- list(ms.pairs, c(mean(x.input), mean(y.output)), pairs)
  return (data)
}
GenerateWigglyNoisy <- function(n) {
  # Takes set of non-linear X's, applies a function, and adds noise.
  #
  # Args:
  #   n: Number of pairs to be generated.
  #
  # Returns:
  #   data: List, with mean-subtracted pairs of inputs X and outputs Y, and
  #     with vector of means.
  x.input <- seq(0, 1, length=100)
  errors <- rnorm(n, 0, 1)
  y.output <- sin(2*pi*x.input) + errors
  
  # Subract out means.
  ms.input <- x.input - mean(x.input)
  ms.output <- y.output - mean(y.output)
  
  ms.pairs <- cbind(ms.output, ms.input)
  pairs <- cbind(y.output, x.input)
  data <- list(ms.pairs, c(mean(x.input), mean(y.output)), pairs)
  return (data)
}

FitCurves <- function(sample.size, bandwidths, GenerateX) {
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
  data <- GenerateX(sample.size)
  pairs <- data[[1]]
  par(bty="l",las=1)
  plot(pairs[,2], pairs[,1], col="gray", xlab="Sample X-Values",
       ylab="Sample Y-Values", main="Curve Fitting with Normal Kernel
       Smoothing and Varying Bandwidths")
  
  # Predict output on test set "x" and plot. (Here, "x" is a uniformly 
  # distributed set of inputs on the range of the training set, designed to
  # show the full density of the estimation.)
  x <- seq(min(pairs[,2]), max(pairs[,2]), .1)
  output1 <- PredictOutput(pairs, x, bandwidths[1])
  output2 <- PredictOutput(pairs, x, bandwidths[2])
  output3 <- PredictOutput(pairs, x, bandwidths[3])
  lines(x, output1, col="steelblue1", lwd="2")
  lines(x, output2, col="dodgerblue3", lwd="2")
  lines(x, output3, col="grey20", lwd="2")
  legend("topright", title="Bandwidth", legend=bandwidths, bty="n",
         lty=c(1,1),lwd=c(3,3), col=c("steelblue1","dodgerblue3","grey20"))
  
  fitted.curves <- list(x, output1, output2, output3, data[[2]])
  return (fitted.curves)
}

PredictOutput <- function(pairs, x.star, bandwidth) {
  # Predicts output y.star for given x.star, based on pairs of data and
  # a weight function. Specifically, doing sum(weight.i*y.i), where each weight
  # comes from a weight function, which here is dnorm(mean=0,sd=1).
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
    weight.i <- dnorm((x.i-x.star)/bandwidth)
    weight.sum <- weight.sum + weight.i
    output.i <- weight.i*y.i
    output <- output + output.i
  }
  # Normalize weights by dividing by weight sum.
  output <- output/weight.sum
  return (output)
}


ComputeMSE <- function(fitted.curves, flag) {
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
      if (flag == "Smooth") {
        expected.value <- (x.values[j]+x.mean)^2
        sq.err <- (predicted.value - expected.value)^2
      } else if (flag == "Wiggly") {
        expected.value <- sin(2*pi*(x.values[j]+x.mean))
        sq.err <- (predicted.value - expected.value)^2
      } else {
        print("Didn't find correct flag in ComputeMSE.")
        quit(save="default")
      }
      sum.sq.err.i <- sum.sq.err.i + sq.err
    }
    
    mse.element <- (1/num.xvalues)*sum.sq.err.i
    mses <- c(mses, mse.element)
  }
  
  return (mses)
}

Main()




