# StatMod2 - Gibbs - Exercises 4

library(MCMCpack)
library(BayesBridge)
library(MASS)

# Get data and initialize constants.
data(diabetes, package="BayesBridge")
y <- diabetes$y
y <- y - mean(y);
X <- diabetes$x
XtX <- t(X)%*%X
Xty <- t(X)%*%y
n <- dim(X)[1]
p <- dim(X)[2]
I.p <- diag(p)
a=1; b=1; c=1; d=1;

# Compare to basic OLS.
fit <- lm(y~ X-1)
summary(fit)

init <- function() { # initialize parameters
  sigmasq <- 1/(rgamma(1, a/2, b/2))
  tausq <- 1/(rgamma(1, c/2, d/2))
  return (theta=list(tausq=tausq, sigmasq=sigmasq))
}

gibbs <- function(n.iter, verbose) {
  THETA <- NULL
  theta <- init()  
  tausq <- theta$tausq
  sigmasq <- theta$sigmasq
  a <- theta$a
  b <- theta$b
  c <- theta$c
  d <- theta$d
  
  for (i in 1:n.iter) { # loop over iterations
    beta <- sample.beta(sigmasq, tausq)
    sigmasq <- sample.sigmasq(beta)
    tausq <- sample.tausq(beta)
    if (verbose > 0){
      if (i %% 10 == 0) {
        # print short summary, for every 10th iteration
        print(paste("Iteration: ", i))
        print(round(z, 2))
        print(round(q, 2))
        print(round(r, 2))
        print(round(pistar, 2))
      }  
    }
    ## Save iteration. First p columns are Beta; followed by sigmasq, tausq.
    THETA <- rbind(THETA, c(beta, sigmasq, tausq))
  }
  return(list(THETA=THETA))
}

sample.beta <- function(sigmasq, tausq) {
  cov <- ginv(XtX/sigmasq+I.p/tausq)
  mean <- cov%*%(Xty/sigmasq)
  beta <- mvrnorm(1, mu=mean, Sigma=cov)
  return (beta)
}
sample.sigmasq <- function(beta) {
  shape <- n/2 + a
  rate <- t(y-X%*%beta)%*%(y-X%*%beta)/2 + b
  sigmasq <- rgamma(1, shape=shape, rate=rate)
  return (1/sigmasq)
}
sample.tausq <- function(beta) {
  shape <- p/2 + c
  rate <- t(beta)%*%(beta)/2 + d
  tausq <- rgamma(1, shape=shape, rate=rate)
  return (1/tausq)
}

ex <- function() { 
  ## RUN these lines to get the plots
  n.iter <- 1000; verbose=0
  gbs <- gibbs(n.iter, verbose)
  THETA <- gbs$THETA
  its <- 1:n.iter
  
  par(mfrow=c(5,2))
  for (i in 1:p) {
    ## trajectory plot of Beta_1 (out of Beta_1,...,Beta_p).
    plot(its, THETA[,p],xlab="ITER",ylab="Beta",bty="l",type="l",
         main=paste("Trajectories of Beta", i))
  }
    
  ## trajectory plot of SigmaSq
  plot(its, THETA[,p+1],xlab="ITER",ylab="SigmaSq",bty="l",type="l",
       main="Trajectories of SigmaSq")
  ## trajectory plot of TauSq
  plot(its, THETA[,p+2],xlab="ITER",ylab="TauSq",bty="l",type="l",
       main="Trajectories of TauSq")
  
}

ex()





