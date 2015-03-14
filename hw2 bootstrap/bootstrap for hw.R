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
  sigma <- matrix(c(2,1,1,1), nrow=dimension, ncol=dimension)
  A <- chol(sigma)
  z <- matrix(rnorm(vector.size), nrow=dimension, ncol=1)
  mu <- matrix(c(1,1), nrow=dimension, ncol=1)
  x <- A %*% z + mu
  return (x)
}

# Sample 300 multivariate random normals, as defined in SampleMvn.
mvn <- do(1000) * SampleMvn(vector.size=2, dimension=2)
plot(mvn, pch=16, col=4, cex=0.5, xlab="x", ylab="y")

# PART B2: For sample of multivariate normals, estimate mean vector and
# covariance matrix by maximum likelihood.

SolveMle <- function(mvn) {
  # Estimates optimal mean vector and covariance matrix, given multivariate
  # normal sample.
  #
  # Args:
  #   mvn: Sample of multivariate normal variables.
  #
  # Returns:
  #   val.loglik: Log likelihood value, used to find optimal mean and covariance.

  # Define log likelihood function to optimize.
  LogLik <- function(pars, mvn) {
    
    # Pass in mean vector and covariance matrix element-wise.
    mu.o <- matrix(c(pars[1], pars[2]), nrow=2, ncol=1)
    sigma.o <- matrix(c(pars[3], pars[4], pars[5], pars[6]), nrow=2, ncol=2)
    
    # Compute helper variables.
    sigma.o.inv <- solve(sigma.o)
    det.sigma.o <- det(sigma.o)
    det.sigma.o.inv <- det(sigma.o.inv)
    n <- dim(mvn)[1]
    k <- dim(mvn)[2]

    # USING FULL QUADRATIC    
    # Compute exponential sum.
    quad <- 0
    for (i in n) {
      quad <- quad + as.matrix(mvn[i,]-mu.o) %*% 
                     sigma.o.inv %*% 
                     as.matrix(t(mvn[i,]-mu.o))
    }
    
    # Compute value of log likelihood (modulo some constant).
    # https://www.cs.cmu.edu/~epxing/Class/10701-08s/recitation/gaussian.pdf
    # First line in "trace trick" section.
    val.loglik <- -1*(n/2)*log(det.sigma.o) - (1/2)*(quad)
    
    return (-val.loglik)

#     ### USE TRACE TRICK ###
#     new.loglik <- -1*(n/2)*log(det.sigma.o) - 
#                   (1/2)*sum(diag(t(as.matrix(mvn[i,]-mu.o)) %*% 
#                                    as.matrix(mvn[i,]-mu.o)) %*% 
#                                    sigma.o.inv)
#     return (new.loglik)
    
#     ### USE TRACE TRICK ###
#     #https://statistics.stanford.edu/sites/default/files/2009-10.pdf
#     # Compute sum for trace component.
#     q <- matrix(0, nrow=2, ncol=2)
#     for (i in 1:n) {
#       q <- q + t(as.matrix(mvn[i,]-mu.o)) %*% as.matrix(mvn[i,]-mu.o)
#     }
#     # Compute value of log likelihood (modulo some constant).
#     val.loglik <- (n/2)*log(det.sigma.o.inv) - (1/2)*sum(diag(sigma.o.inv %*% q/n))

  }
  
  # Use "optim" function to optimize over mean vector and covariance matrix.
  start.mu <- matrix(c(5,5), nrow=2, ncol=1)
  start.sigma <- matrix(c(3,1,1,2), nrow=2, ncol=2)
  optimization <- optim(c(start.mu, start.sigma), LogLik, mvn=mvn)#, 
                        #control=list(fnscale=-1000, parscale=c(1,1,.25,.25,.25,.25)))
  
  return (optimization)
}

# Run "optim" optimization.
o <- SolveMle(mvn)
o

# Cheaply use pre-packaged function to get MLE.
# mlest(mvn)

# PART B3: Bootstrap a given sample to estimate the sampling distribution of the MLE.
# First create master sample.
n <- 500
mvn <- do(n) * SampleMvn(vector.size=2, dimension=2)

# Set up matrix to store values of muhat and sigmahat for each resampling.
booted.mvn <- matrix(-99, nrow=n, ncol=6)

# Do each resample, and store values.
for(i in 1:n) {
  # Resample original mvn dataset.
  myindices <- sample(1:nrow(mvn), nrow(mvn), replace=TRUE)
  resample <- mvn[myindices,]
  
  # Estimate MLE on resampled set.
  mle <- mlest(resample)
  booted.mvn[i,] <- c(mle$muhat, mle$sigmahat)
}

# Average muhat and sigmahat values across resamplings.
bootstrapped.estimators <- colMeans(booted.mvn)
booted.muhat <- matrix(bootstrapped.estimators[1:2], 2, 1)
booted.sigmahat <- matrix(bootstrapped.estimators[3:6], 2, 2)

booted.muhat
booted.sigmahat

colMeans(mvn)
cov(mvn)
