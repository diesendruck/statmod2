# Stat Mod 2 - Exercises 2 - Smoothing (p.3)


GeneratePairs <- function(n) {
  # Takes set of non-linear X's and applies the y = x^2 + e transformation,
  # where error values "e" are Normal(0,1). Then subracts means from each.
  #
  # Args:
  #   n: Number of pairs to be generated.
  #
  # Returns:
  #   pairs: Pairs of inputs X and outputs Y, where Y is a non-linear function
  #     of X, with added error.
  x.input <- rnorm(n, 5, 2)
  errors <- rnorm(n, 0, 1)
  y.output <- x.input^2 + errors
  
  # Subract out means.
  ms.input <- x.input - mean(x.input)
  ms.output <- y.output - mean(y.output)
  
  ms.pairs <- cbind(ms.output, ms.input)
  pairs <- cbind(y.output, x.input)
  return (ms.pairs)
}

PredictOutput <- function(x.star, pairs, bandwidth) {
  # Predicts probability density for given x.star, based on pairs of data and
  # a weight function. Specifically, doing sum(weight.i*y.i), where each weight
  # comes from WeightFunction.
  #
  # Args:
  #   x.star: This is the point of interest, whose probability under the
  #     density you want to find.
  #   pairs: The pairs of X's and corresponding Y's, i.e. the data.
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

# Generate new pairs of original data.
par(mfrow=c(2,1))
pairs <- GeneratePairs(10000)
plot(pairs[,2], pairs[,1])
#hist(pairs[,2])
#hist(pairs[,1])

# Predict output.
x <- seq(-5, 5, .1)
output <- PredictOutput(x, pairs, 0.1)
plot(x, output)
