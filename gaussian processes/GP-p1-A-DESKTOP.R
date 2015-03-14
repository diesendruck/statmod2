# StatMod 2 - Gaussian Processes - A - (p.1)

library(MASS)
library(graphics)
library(mvtnorm)
library(psych)
library(Matrix)
require(grDevices)

Main <- function() {
  # Runs main functions.
  #
  # Args:
  #   NA: None.
  #
  # Returns:
  #   NA: Plots results.
  
  # ------ Regular Gaussian Processes ------ #
  
#   for (i in seq(0.1, 1, 0.3)) {
#     for (j in seq(0.1, 1, 0.1)) {
#       par(mfrow=c(2,2))
#       inputs <- seq(0, 1, 0.01)
#       hyperparams <- c(i, j, 0.0000001)
#       z1 <- GaussianProcess(inputs, hyperparams, SqExp, Y=NULL)[[1]]
#       z2 <- GaussianProcess(inputs, hyperparams, Matern52, Y=NULL)[[1]]
#       #z3 <- GaussianProcess(inputs, hyperparams, LinearPlane, Y=NULL)[[1]]
#     }
#   }
  
  # ------ Gaussian Processes (posterior) Part C - (p.3) ------ #
  
  #par(mfrow=c(2,1))
  d <- GetData()
  hyperparams <- c(120, 68, 0) # Paramaters solved by optim.
  #hyperparams <- c(1, 0.1, 0) # Parameters chosen by Mo.
  inputs <- d$x
  Y <- d$y
  n <- length(Y)
  results <- GaussianProcess(inputs, hyperparams, Matern52, Y)
  posterior.mean.and.cov <- results[[1]]
  post.c <- results[[2]]
  
  # Do local linear polynomial smoothing on posterior data.
  post.means <- posterior.mean.and.cov[,1]
  post.covs <- posterior.mean.and.cov[,2]
  post.data <- data.frame(cbind(inputs, post.means))
  names(post.data) <- c("x", "y")
  hbs <- PerformLOOCV(post.data, bandwidth=NULL)
  H <- hbs[[1]]
  bw <- hbs[[2]]
  yhat <- H %*% d$y
  r <- CalcResid(post.data, yhat, H, bw, opt=post.covs)
  
  # ------ Gaussian Processes (marginal likelihood) Part E - (p.3) ------ #
  
  d <- GetData()
  inputs <- d$x
  Y <- d$y
  n <- length(Y)
  
  # Find ideal tuning parameters under marginal likelihood. "OPTIM" approach.
  LML <- function(params) {
    hyperparams <- c(params[1], params[2], 0)
    c <- Cov(inputs, inputs, hyperparams, Matern52)
    hbs <- PerformLOOCV(d, bandwidth=NULL)
    H <- hbs[[1]]
    sig.sq.hat <- as.matrix(Diagonal(n, x=(t(Y)%*%(diag(n)-H)%*%Y)/
                                          (n-2*tr(H)+tr(t(H)%*%H))))
    v <- dmvnorm(Y, sigma=c+sig.sq.hat, log=T)
    return (-v)
  }
  optim(c(5, 1), LML)

  # ------ Gaussian Processes (posterior contours) Part F - (p.3) ------ #
  
  pressure.hyperparams <- c(.5,5,150000,0)
  temperature.hyperparams <- c(4,.5,20,0)
  pr.means <- WeatherContours("Pressure", pressure.hyperparams)[[1]]
  te.means <- WeatherContours("Temperature", temperature.hyperparams)[[1]]
  
}

WeatherContours <- function(name, hyperparams) {
  # Makes contours for weather data.
  #
  # Args:
  #   name: The name of the variable.
  #   hyperparams: A vector of hyperparameters: b1, b2, t1.sq, t2.sq.
  #
  # Returns:
  #   z: List of (1) matrix of posterior means by lon/lat, (2) matrix of
  #     posterior covariances by lon/lat, and (3) overall covariance matrix.
  
  #setwd("~/Google Drive/2. SPRING 2015/STAT MOD 2 - Prof Scott/gaussian processes/")
  setwd("~/Documents")
  d <- read.csv("weather.csv", header=T, sep=",")
  d <- d[order(d[, 1]), ]
  b1 <- hyperparams[1]
  b2 <- hyperparams[2]
  tau1sq <- hyperparams[3]
  tau2sq <- hyperparams[4]
  lon <- d$lon
  lat <- d$lat
  
  if (name=="Pressure") {
    names(d) <- c("y", "temperature", "lon", "lat")
  } else if (name=="Temperature") {
    names(d) <- c("pressure", "y", "lon", "lat")
  }
  
  y <- d$y
  
#   SqExpNonIso <- function(b1,b2,tau1sq,tau2sq){
#     function(lon, lat){
#       nonIsoDist <- as.matrix(dist(lon, diag=T, upper=T)/b1 + dist(lat, diag=T, upper=T)/b2)
#       kron.delta <- diag(nrow=length(lon))
#       c <-tau1sq*exp(-.5*nonIsoDist^2) + tau2sq*kron.delta
#       return (c)
#     }
#   }
  SqExpNonIso <- function(b1,b2,tau1sq,tau2sq){
    function(lon, lat){
      nonIsoDist.sq <- as.matrix((dist(lon, diag=T, upper=T))^2/b1^2 + 
                                 (dist(lat, diag=T, upper=T))^2/b2^2)
      kron.delta <- diag(nrow=length(lon))
      c <-tau1sq*exp(-.5*nonIsoDist.sq) + tau2sq*kron.delta
      return (c)
    }
  }
  
  # Here, sigma^2 = var(y), rather than the (error^2)/(n-2tr(H)+tr(H'H)); ???
  error.variance <- as.matrix(Diagonal(length(y), x=var(y)))
  CovFun <- SqExpNonIso(b1, b2, tau1sq, tau2sq)
  c <- CovFun(lon, lat)
  post.c <- c + error.variance
  
  # For each point, compute posterior mean and covariance.
  grid.lon <- seq(-132, -113, 1)
  grid.lat <- seq(39, 53, 1)
  n <- length(y)
  n1 <- length(grid.lon)
  n2 <- length(grid.lat)
  post.mean <- matrix(NA, nrow=n1, ncol=n2)
  post.cov <- matrix(NA, nrow=n1, ncol=n2)
  
  for (i in 1:n1) {
    for (j in 1:n2) {
      # Compute posterior mean = u_x* + K(X*,X) %*% (K(X,X)+sigsq*I)^(-1) %*% (y-u_y)
      predictive.c <- as.matrix(CovFun(c(grid.lon[i], lon), c(grid.lat[j], lat)))
      k.xs.x <- as.matrix(predictive.c[1,2:(n+1)])
      k.xs.xs <- as.matrix(CovFun(c(grid.lon[i], grid.lon[i]),
                                  c(grid.lat[j], grid.lat[j])))[1,1]
      post.mean[i, j] <- t(k.xs.x)%*%ginv(post.c)%*%(y)
      post.cov[i, j] <- k.xs.xs - t(k.xs.x)%*%ginv(post.c)%*%k.xs.x
    }
  }
  z <- list(post.mean, post.cov, post.c)
  
  filled.contour(grid.lon, grid.lat, post.mean, color=terrain.colors,
                 plot.title=title(main=paste(name, ", Params: ",
                                             hyperparams[1],
                                             ", ", hyperparams[2],
                                             ", ", hyperparams[3],
                                             ", ", hyperparams[4], sep=""),
                                  xlab="Lon", ylab="Lat"),
                 key.title=title(main="Dif"))
  
  return (z)  
}

GaussianProcess <- function(inputs, hyperparams, kernel, Y) {
  # Runs main operations to perform Gaussian Process (GP), and plots results.
  #
  # Args:
  #   inputs: A vector of scalar points.
  #   hyperparams: A vector of three hyperparameters, (b, tau1^2, tau2^2)
  #   Y: A flag with Y values, used to add diagonal error variance to 
  #     covariance matrix, and to compute posterior mean.
  #
  # Returns:
  #   z: List of the values produced by the Gaussian Process, and the
  #     covariance matrix "c". If Y is provided for posterior mean and
  #     covariance, then z is a matrix with a column for posterior mean, and
  #     another for posterior covariance.
  
  # Compute covariance matrix and use it to get GP values/output.
  c <- Cov(inputs, inputs, hyperparams, kernel)
  
  # ------ Posterior Calculations (Y!=NULL) ------ #
  
  if (!is.null(Y)) {
    # Compute K(X,X)+sigsq*I.
    d <- as.data.frame(cbind(inputs, Y))
    names(d) <- c("x","y")
    hbs <- PerformLOOCV(d, bandwidth=NULL)
    H <- hbs[[1]]
    sig.sq.hat <- as.matrix(Diagonal(n, x=(t(Y)%*%(diag(n)-H)%*%Y)/
                                          (n-2*tr(H)+tr(t(H)%*%H))))
    post.c <- c + sig.sq.hat
  
    # Make storage for posterior means and covariances.
    posterior.mean.and.cov <- matrix(0, nrow=length(inputs), ncol=2)
    
    # For each point, compute posterior mean and covariance.
    n <- length(inputs)
    for (i in 1:n) {
      # Compute posterior mean = u_x* + K(X*,X) %*% (K(X,X)+sigsq*I)^(-1) %*% (y-u_y)
      k.xs.x <- as.matrix(Cov(inputs[i], inputs, hyperparams, kernel))
      post.mean.i <- Y[i] + k.xs.x%*%ginv(post.c)%*%(Y-mean(Y)) # With mu's added.
      #post.mean.i <- k.xs.x%*%ginv(post.c)%*%(Y) # Experimenting with not adding mu's.
      post.cov.i <- Cov(inputs[i], inputs[i], hyperparams, kernel) - 
        k.xs.x%*%ginv(post.c)%*%t(k.xs.x)
      posterior.mean.and.cov[i,] <- c(post.mean.i, post.cov.i)
    }
    z <- list(posterior.mean.and.cov, post.c)
    return (z)
  }
  
  # ------ Regular Calculations (Y=NULL) ------ #
  
  if (is.null(Y)) {
    n <- length(inputs)
    u <- rnorm(n)
    decomp <- svd(c)
    g <- decomp$u %*% sqrt(diag(decomp$d)) %*% u
    
    # Plot the output.
    plot(inputs, g, main=paste("b=", hyperparams[1],
                               " t1.sq=", hyperparams[2],
                               " t2.sq=", hyperparams[3]))
    lines(inputs, g)
    hist(g, 40)
    z <- list(g, c)
    return (z)
  }

}

LogMarLik <- function(y, c) {
  # Computes log marginal likelihood.
  #
  # Args:
  #   y: The vector of Y values.
  #   c: The covariance matrix, calculated as part of the Gaussian Process.
  #
  # Returns:
  #   lml.value: The value of the log marginal likelihood.
  lml.value <- (-1/2)*t(y)%*%c%*%y - (1/2)*log(det(c)) - length(y)/2*log(2*pi)
  return (lml.value)
}

Cov <- function(ins, inputs, hyperparams, kernel) {
  # Computes the covariance of a finite set of points.
  #
  # Args:
  #   ins: The vector of points compared to inputs for covariance.
  #   inputs: The vector of scalar points of the data.
  #   hyperparams: A vector of three hyperparameters, (b, tau1^2, tau2^2)
  #
  # Returns:
  #   cov.matrix: A symmetrix, positive semi-definite covariance matrix. The
  #     dimensions of this matrix are: length(ins) x length(inputs)
  b <- hyperparams[1]
  t1.sq <- hyperparams[2]
  t2.sq <- hyperparams[3]
  
  # Compute covariance with specific covariance function.
  cov.matrix <- matrix(NA, nrow=length(ins), ncol=length(inputs))
  n1 <- length(ins)
  n2 <- length(inputs)
  for (i in 1:n1) {
    for (j in 1:n2) {
      cov.matrix[i, j] <- kernel(ins[i], inputs[j], b, t1.sq, t2.sq)
    }
  }
  return (cov.matrix)
  
  # Svelte, but hard to understand. "Tell me about it..."
  #   cov <- as.matrix(t(sapply(
  #     ins, function(k) kernel(k, inputs, b, t1.sq, t2.sq))))
  #   return (cov)
}

GetData <- function() {
  # Prepares data from file.
  #
  # Args:
  #   NA: None.
  #
  # Returns:
  #   data: A matrix of data.
  #setwd("~/Google Drive/2. SPRING 2015/STAT MOD 2 - Prof Scott/local polynomial regression")
  setwd("~/Documents")
  data <- read.csv("utilities.csv", header=T, sep=",")
  
  # Add column for average daily bill.
  data <- within(data, {
    y <- gasbill/billingdays
    x <- temp
  })
  
  data <- data[order(data[, 1]), ]
  attach(data)
  plot(x, y)
  
  return (data)
}

PerformLOOCV <- function(data, bandwidth) {
  # Computes hat matrix and LOOCV score for a given bandwidth; or, if no
  # bandwidth is chosen, finds optimal bandwidth (one with lowest LOOCV score)
  # and then computes hat matrix.
  # 
  # Args:
  #   data: The full data set.
  #   bandwidth: The bandwidth for the model. Find optimal by leaving NULL.
  #
  # Returns:
  #   hat.bw.score: A list of hat matrix, bandwidth, and loocv score.
  x <- data[, c("x")]
  y <- data[, c("y")]
  
  hat.matrix <- function(bw) t(sapply(x, function(z) Weight(z, bw)))
  loocv.score <- function(bw) {
    H <- hat.matrix(bw)
    yhat <- H%*%y
    score <- sum(((y-yhat)/(1-diag(H)))^2)
    return (score)
  }
  
  if (!is.null(bandwidth)) {
    bw <- bandwidth
  } else {
    bw <- optimize(loocv.score, interval=c(0, max(x)-min(x)))$minimum
  }
  
  H <- hat.matrix(bw)
  score <- loocv.score(bw)
  hat.bw.score <- list(H, bw, score)
  return (hat.bw.score)
}

Weight <- function(z, bandwidth, kernel=dnorm) {
  # For estimation at point z, gives the Gaussian-smoothed weights that will be
  # applied to every x (including itself). This is equivalent to a row of the
  # Hat Matrix.
  #
  # Args:
  #   z: Vector of X values to receive weights.
  #   bandwidth: A neighborhood around each requested x.star, used to average 
  #     the area under the smoother function.
  #   kernel: The function used to define the type of smoothing around z.
  #
  # Returns:
  #   weight: The weight of the y that corresponds to x.i.
  dif <- x-z
  weight <- kernel(dif/bandwidth)/bandwidth
  s1 <- sum(weight*dif)
  s2 <- sum(weight*dif^2)
  weight <- weight*(s2-dif*s1)
  weight <- weight/sum(weight)
  return (weight)
}

CalcResid <- function(data, yhat, H, bw, opt) {
  # Calculates and plots residuals for yhat.
  # 
  # Args:
  #   data: The full data set.
  #   yhat: The predictions generated by multiplying the hat matrix with y.
  #   H: The hat matrix.
  #   bw: The bandwidth.
  #   opt: An optional parameter, used here to pass bootstrapped variances.
  #
  # Returns:
  #   r: Vector of residuals.
  x <- data[, c("x")]
  y <- data[, c("y")]
  n <- length(y)
  
  # Compute residuals for y prediction.
  r <- y-yhat
  
  # Compute variance using formula (probably p.6 of Ex3).
  if (!is.null(opt)) {
    variance <- opt
  } else {
    variance <- sum(r*r)/(n-2*sum(diag(H))-sum(diag(t(H)%*%H)))
  }
  se <- sqrt(variance*rowSums(H^2))
  
  # Plot residuals.
  plot(x, r, xlab="Temp", ylab="Res. on Avg. Daily Bill",
       main=paste("Residuals for Bandwidth =", round(bw, 2)))
  
  # Plot prediction with 95% confidence interval lines.
  plot(x, y, col="gray", xlab="Temp", ylab="Y",
       main=paste("Local Polynomial Regression with 95% CI, Bandwidth =",
                  round(bw, 2)))
  t <- qt(c(0.025, 0.975), 116)
  lines(x, yhat, col="grey20", lwd="2")
  lines(x, yhat+t[1]*se, col="dodgerblue3", lwd="2", lty=3)
  lines(x, yhat+t[2]*se, col="steelblue1", lwd="2", lty=3)
  
  return (r)
}

# Define available covariance functions, and a helper Delta function.
SqExp <- function(i1, i2, b, t1.sq, t2.sq) {
  d <- i2 - i1
  val <- t1.sq * exp((-1/2)*(d/b)^2) + (t2.sq)*Delta(i2, i1)
  return (val)
}
Matern52 <- function(i1, i2, b, t1.sq, t2.sq) {
  d <- sqrt(sum((i2 - i1)^2))
  val <- t1.sq*(1 + sqrt(5)*d/b + 5/3*(d^2)/(b^2)) * exp(-sqrt(5)*d/b) + (t2.sq)*Delta(i2, i1)
  return (val)
}
LinearPlane <- function(i1, i2, b, t1.sq, t2.sq) {
  val <- i1*i2
  return (val)
}
Delta <- function(x1, x2) {
  ifelse(x1==x2, 1, 0)
} 

Main()

