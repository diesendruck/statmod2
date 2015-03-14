library(mosaic)
library(mlbench)
library(mvnmle)

# PART A: Compute sigma.hat for ozone and compare it to parametric estimate.

# Load the data.
ozone <- data(Ozone, package='mlbench')

# Extract the relevant columns.
ozone <- na.omit(Ozone)[,4:13]
y <- ozone[,1]
x <- as.matrix(ozone[,2:10])
x <- cbind(1, x)

# Compute the estimator.
betahat <- solve(t(x) %*% x) %*% t(x) %*% y

# Compute beta covariance matrix by normal theory.
pred <- x %*% betahat
sum.squared.errors <- sum((y-pred)^2)
n <- length(y)
p <- length(betahat) - 1
beta.variance <- sum.squared.errors / (n-(p+1))
xtx <- t(x) %*% x
xtx.inv <- solve(xtx)
betacov <- beta.variance * xtx.inv

# Create a placeholder array for var-cov values for each bootstrapped sample.
boot1 <- matrix(0, nrow=1000, ncol=10)

# Go sample by sample in an explicit loop. Save each set of coefficients in a
# single layer of boot1.
for(i in 1:1000)
{
  # Resample, and set up temporary y and x variables.
  myindices <- sample(1:nrow(ozone), nrow(ozone), replace=TRUE)
  resample <- ozone[myindices,]
  
  # Run linear model on resampled set.
  lmtemp <- lm(V4~V5+V6+V7+V8+V9+V10+V11+V12+V13, data=resample)
  #lmtemp <- lm(y~x-1, data<-resample)
  
  # Compute beta covariance matrix using boostrapped betahat.
  boot1[i,] <- coef(lmtemp)
}

# Compute covariance of bootstrapped betas.
boot.cov <- cov(boot1)
# Compare to parametric estimate based on normal theory.
betacov

# PART B1: For specified mean vector mu and covariance matrix sigma, simulate
# multivariate normal random variables.

SampleMvn <- function(vector.size, dimension) {
  # Simulates multivariate normal random variable by construction of standard
  # normal variable "z", times "A" from covariance matrix sigma=AA', plus mean
  # vector "u".
  #
  # Args:
  #   vector.size: Length of multivariate normal vector.
  #   dimension: Dimension of square covariance matrix.
  #
  # Returns:
  #   x: Multivariate normal random variable
  sigma <- matrix(c(2, 1, 1, 1), nrow=dimension, ncol=dimension)
  A <- chol(sigma)
  z <- matrix(rnorm(vector.size), nrow=dimension, ncol=1)
  mu <- matrix(1, nrow=dimension, ncol=1)
  x <- A %*% z + mu
  return (x)
}

# Sample 300 multivariate random normals, as defined in SampleMvn.
mvn <- do(100) * SampleMvn(vector.size=2, dimension=2)

# PART B2: For sample of multivariate normals, estimate mean vector and
# covariance matrix by maximum likelihood.

SolveMle <- function(mvn.sample) {
  # Estimates optimal mean vector and covariance matrix, given multivariate
  # normal sample.
  #
  # Args:
  #   mvn.sample: Sample of multivariate normal variables.
  #
  # Returns:
  #   theta: Vector of estimates for mean and covariance.
  
  # Define log likelihood function to optimize.
  LogLik <- function(pars, data=mvn.sample) {
    
    # Pass in mean vector and covariance matrix element-wise.
    mu.o <- matrix(c(pars[1], pars[2]), nrow=2, ncol=1)
    sigma.o <- matrix(c(pars[3], pars[4], pars[5], pars[6]), nrow=2, ncol=2)
    
    # Compute helper variables.
    sigma.o.inv <- solve(sigma.o)
    det.sigma.o <- det(sigma.o)
    det.sigma.o.inv <- det(sigma.o.inv)
    n <- dim(data)[1]
    k <- dim(data)[2]
    
    # Compute sum for trace component.
    q <- c()
    for (i in 1:n) {
      q <- c(q, as.matrix(data[i,]-mu) %*% t(as.matrix(data[i,]-mu)))
    }
    
    # Compute value of log likelihood (modulo some constant).
    val.loglik <- (-n/2)*log(det.sigma.o) - (1/2)*sum(diag(sigma.o.inv * sum(q)))
    
    return (val.loglik)
  }
  
  # Use "optim" function to optimize over mean vector and covariance matrix.
  start.mu <- matrix(0.3, nrow=2, ncol=1)
  start.sigma <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow=2, ncol=2)
  optimization <- optim(c(start.mu, start.sigma), LogLik)
  
  return (optimization)
}

# Run "optim" optimization.
o <- SolveMle(mvn)
o

# Cheaply use pre-packaged function to get MLE.
# mlest(mvn)



