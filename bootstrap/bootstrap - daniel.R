
#---- Bootstrap ---------------------------------------------------------------#

sigma.bootstrap <- function(y, X, R=1000){
  sigma <- function(){
    k <- sample(length(y), replace=TRUE)
    vcov(lm(y ~ X - 1, subset=k))
  }
  apply(replicate(R, sigma()), 1:2, mean)
}

sigma.boot <- sigma.bootstrap(y, X)
signif(cbind(boot=sqrt(diag(sigma.boot)), lm=sqrt(diag(betacovlm))), 3)

#---- Multivariate normal -----------------------------------------------------#

rmvnorm <- function(n, mu=rep(0,nrow(sigma)), sigma=diag(length(mu))){
  L <- chol(sigma)
  Z <- matrix(rnorm(n * length(mu)), nrow=length(mu))
  t(L %*% Z + mu)
}

mu <- 1:2
sigma <- matrix(c(1,0.5,0.5,1),2,2)
plot(rmvnorm(1000, mu, sigma), pch=16, col=4, cex=0.5, xlab="x", ylab="y")
abline(v=mu[1], h=mu[2])

mle <- function(Y, mu=rep(0,ncol(Y)), sigma=diag(ncol(Y))){
  p <- ncol(Y)
  
  disassemble <- function(mu, sigma){
    sd <- sqrt(diag(sigma))
    cor <- cov2cor(sigma)
    c(mu, sd, cor[upper.tri(cor)])
  }
  
  assemble <- function(x){
    mu <- x[1:p]
    sd <- x[1:p+p]
    rho <- x[-(1:(p+p))]
    cor <- diag(p)
    cor[upper.tri(cor)] <- rho
    cor[lower.tri(cor)] <- t(cor)[lower.tri(cor)]
    list(mu=mu, sigma=cor*outer(sd,sd))
  }
  
  objective <- function(mu, sigma){
    omega <- solve(sigma)
    d <- log(det(sigma))
    sum(apply(Y, 1, function(y){ d + t(y-mu) %*% omega %*% (y-mu) }))
  }
  
  assemble(optim(disassemble(mu,sigma),
                 function(x){ do.call(objective, assemble(x)) })$par)
}


mle.out <- lapply(c(100,1000,10000), function(n){ mle(rmvnorm(n, mu, sigma)) })

mle.bootstrap <- function(Y, R=100){
  ml <- function(){
    k <- sample(nrow(Y), replace=TRUE)
    mle(Y[k,])
  }
  replicate(R, ml())
}

mle.bootstrap.out <- mle.bootstrap(rmvnorm(100, mu, sigma))
mus <- do.call(rbind, mle.bootstrap.out["mu",])
plot(mus, pch=16, col=4, cex=0.5, xlab="x", ylab="y")
abline(v=mu[1], h=mu[2])

sigmas <- simplify2array(mle.bootstrap.out["sigma",])
par(mfcol=c(2,1))
hist(sqrt(sigmas[1,1,]), xlab="Sigma11", main="")
hist(sqrt(sigmas[2,2,]), xlab="Sigma22", main="")
par(mfcol=c(1,1))
