# StatMod2 - Multinomial-Dirichlet

library(MCMCpack)

#####
# QUESTION 1
# Simulate school test scores: The school-wide results appear unimodal and
# slightly skewed left.
simulate.scores <- function() {
  y.rem <- rnorm(100, 55, 15)
  y.avg <- rnorm(400, 70, 10)
  y.hon <- rnorm(150, 85, 5)
  school.wide <- c(y.rem, y.avg, y.hon)
  hist(school.wide)
}
simulate.scores()

#####
# QUESTION 2
# Write functions to draw from full conditional posterior distributions of
# gamma.i and w.
sample.gamma <- function(y.i, w) {
  if (w==w1) {
    p <- w*dnorm(y.i, 55, 15)
  } else if (w==w2) {
    p <- w*dnorm(y.i, 70, 10)
  } else if (w==w3) {
    p <- w*dnorm(y.i, 85, 5)
  }
  return (p)
}
sample.w <- function(group.counts=c(n1, n2, n3), dir.params=c(a1, a2, a3)) {
  rdirichlet(1, c(group.counts[1]+dir.params[1],
                  group.counts[2]+dir.params[2],
                  group.counts[3]+dir.params[3]))
}
w <- sample.w(group.counts, dir.params)
w1 <- w[1]
w2 <- w[2]
w3 <- w[3]
g <- sample.gamma(90, w2)

#####
# QUESTION 3 - on paper

#####
# GIBBS SAMPLER

require("gtools") ## found this package on-line
## for generating from a Dirichlet distribution
require("coda") ## for convergence diagonstics
## read in data 
x <- scan("mid1-2.dta") 
n <- length(x) 
J <- 5 
xgrid <- seq(from=0,to=1,length=100) 
## hyperparameters 
tau2 <- 1 
a <- 5 # 1/sig2 ~ Ga(5, 0.05) s.t. E(1/sig2) = 100 
b <- 0.05 
alpha <- 1

sample.s <- function(w,mu,sig2) { 
  ## sample s[i] from p(s[i] | ...) 
  n <- length(x) 
  sd <- sqrt(sig2) 
  s <- rep(0,n) # initialize 
  for(i in 1:n){
    pr <- w*dnorm(x[i],m=mu,sd=sd) 
    s[i] <- sample(1:J,1,replace=T,prob=pr) 
  } 
  return(s) 
}

sample.mu <- function(s,w,sig2) { 
  ## sample mu[j] mu <- rep(0,J) # initialize 
  for (j in 1:J) {
    Aj <- which(s==j) 
    nj <- length(Aj) 
    if (nj==0) { 
      m <- 0; V <- tau2
    } else { 
      xbar <- mean(x[Aj]) 
      V <- 1/(1/tau2 + nj/sig2) 
      m <- V*xbar*nj/sig2 
    } 
    mu[j] <- rnorm(1,m=m,sd=sqrt(V))
  } 
  return(mu) 
}

sample.w <- function(s,sig2) {
  ## sample w 
  a <- rep(0,J) # initialize 
  for(j in 1:J) { 
    Aj <- which(s==j) 
    nj <- length(Aj) 
    a[j] <- alpha+nj 
  } 
  w <- rdirichlet(1, a) 
  w <- c(w) 
  return(w)
}

sample.sig <- function(s,w,mu) {
  ## sample sig2 
  S2 <- sum((x-mu[s])^2) 
  a1 <- a+n/2 
  b1 <- b+S2/2 
  sig2 <- 1.0/rgamma(1,shape=a1,rate=b1) 
  return(sig2)
}

f <- function(xi, w,mu,sig2) { 
  y <- 0 
  for(j in 1:J){
    y <- y+w[j]*dnorm(xi,m=mu[j],sd=sqrt(sig2)) 
  } 
  return(y)
}

init <- function() { 
  ## use some exploratory data analyssi ## for initial values of the parameters
  hc <- hclust(dist(x), "ave") # initialize s with a 
  s <- cutree(hc,k=J) # hierarchical tree, cut for J=5 clusters 
  ## this is not important -- you could instead initialize arbitrarily.. 
  ## google "hierarchical clustering" to learn more. 
  mu <- lapply(split(x,s), mean) # iniatialize mu as mean
                                 ## for each cluster 
  mu <- unlist(mu) # to make it a vector 
  sig2 <- var(x-mu[s]) 
  w <- table(s)/n 
  w <- c(w) 
  return(th=list(s=s,mu=mu,sig2=sig2,w=w)) 
}
    
gibbs <- function(n.iter=1000) { 
  ## initialize the parameters 
  th <- init() 
  s <- th$s; mu <- th$mu; sig2 <- th$sig2; w <- th$w
  ## set up lists to save simulations 
  slist <- NULL 
  flist <- NULL 
  siglist <- NULL 
  wlist <- NULL 
  lplist <- NULL
  
  for(iter in 1:n.iter){
    s <- sample.s(w,mu,sig2)
    mu <- sample.mu(s,w,sig2) 
    w <- sample.w(s,sig2) 
    sig2 <- sample.sig(s,w,mu) 
    ## save summaries of current simulation 
    lp <- pth(w,mu,sig2) 
    lplist <- c(lplist,lp) 
    slist <- rbind(slist,s) 
    flist <- rbind(flist, f(xgrid,w,mu,sig2)) 
    siglist <- c(siglist,sig2) 
    wlist <- rbind(wlist,w) 
  } 
  return(list(s=slist,f=flist,sig=siglist,w=wlist,lp=lplist)) 
}




