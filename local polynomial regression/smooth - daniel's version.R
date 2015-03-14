# Daniel Lewis Mitchell
# Feb 25, 2015

# read data
setwd("~/Google Drive/2. SPRING 2015/STAT MOD 2 - Prof Scott/local polynomial regression")
utils <- read.csv("utilities.csv")

# create response and explanatory variables
utils <- within(utils, {
  y <- gasbill / billingdays
  x <- temp
})

# sort data (for lines)
utils <- utils[order(utils$x),]

# local linear estimation, if bandwidth is not given it chosen by LOOCV
smooth <- function(y, x, kernel=dnorm, bandwidth=NULL){

  # compute normalized weights at z (this is one row of H)
  weight <- function(z,h){
    u <- x-z
    w <- kernel(u/h)/h
    s1 <- sum(w*u)
    s2 <- sum(w*u*u)
    w <- w*(s2-u*s1)
    w/sum(w)
  }

  # compute weights for all x (this is H)
  hat.matrix <- function(h) t(sapply(x, function(z) weight(z,h)))

  # loocv objective function
  loocv.error <- function(h){
    H <- hat.matrix(h)
    yhat <- H %*% y
    sum(((y - yhat)/(1-diag(H)))^2)
  }

  # bandwidth and predictions
  if (!is.null(bandwidth)){
    h <- bandwidth
  } else {
    h <- optimize(loocv.error, c(0,max(x)-min(x)))$minimum
  }
  H <- hat.matrix(h)
  yhat <- H%*%y

  # confidence bands
  r <- y - yhat
  sigma <- sum(r*r)/(length(y)-2*sum(diag(H))-sum(diag(t(H)%*%H)))
  sigma <- sqrt(sigma*rowSums(H*H))

  # return a smooth object
  structure(list(x=x,y=y,yhat=yhat,h=h,mse=mean((y-yhat)^2), sd=sigma), class="smooth")
}

# plot a smooth object with confidence bands and residuals
plot.smooth <- function(A){
  layout(matrix(1:2,ncol=1))
  with(A, {
    # scatter plot
    plot(x, y, pch=16, cex=0.5, xlab="", ylab="")
    title(sprintf("h=%.3g, mse=%.3g", h, mse))
    # fitted values and confidence bands
    lines(x, yhat, lwd=2, col=4)
    lines(x, yhat+2*sd, lty=2, col=2)
    lines(x, yhat-2*sd, lty=2, col=2)
    # residual plot
    plot(yhat, y-yhat, pch=16, cex=0.5, xlab="", ylab="", main="residuals vs fitted values")
  })
  layout(matrix(1))
}

with(utils, plot(smooth(y,x)))
