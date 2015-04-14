# StatMod2 - Cheese - JAGS Version

library(rjags)
library(boot)
library(VGAM)
library(MASS)

# Fetch and prepare data.
path <- paste("~/Google Drive/2. SPRING 2015/STAT MOD 2 - Prof Scott/",
              "hierarch-shrinkage/cheese", sep="")
setwd(path)
data <- read.csv("cheese.csv")
data$log.Q <- log(data$vol)
data$log.P <- log(data$price)
attach(data)
n <- dim(data)[1]
store.names <- unique(data$store)

# Execute Gibbs sampler for all stores.
Execute <- function() {
  num.stores <- length(store.names)
  # Make container for store name, and means of 4 vars: log.a, b, g, d.
  PARAMS.BY.STORE <- matrix(NA, nrow=num.stores, ncol=5)
  for (i in 1:num.stores) {
    name <- toString(store.names[i])
    store.data <- data[which(data$store==name),]
    PARAMS.BY.STORE[i,] <- c(name, Gibbs(store.data, name))
  }
  return (PARAMS.BY.STORE)
}

Gibbs <- function(store.data, name) {
  # Set up JAGS model over 3 variables.
  jagsmodel <- jags.model(file="cheese-jags-model.txt",
                           data=list(LOG.Q=store.data$log.Q,
                                     log.P=store.data$log.P,
                                     disp=store.data$disp,
                                     cross=store.data$log.P*store.data$disp,
                                     N=length(store.data$log.Q)),
                           n.chains=2)
  # Fit jagsmodel.
  jagsfit <- jags.samples(jagsmodel,
                          variable.names=c("beta", "sigsq", "tau1", "tau2",
                                           "tau3", "tau4"),
                          n.iter=100000, thin=1000)
  
  # Create traceplots for beta chains.
  par(mfrow=c(4, 2))
  ref <- c("log.a", "b", "g", "d")
  for (i in 1:4) {
    # jagsfit$var takes three arguments: parameter, row, chain.
    plot(jagsfit$beta[i,,1], type="l", col="red", ylab="",
         main=bquote(.(name) ~ .(ref[i]) == .(i)))
    lines(jagsfit$beta[i,,2], type="l", col="blue")
    plot(density(jagsfit$beta[i,,]), type="l",
         main=bquote(.(name) ~ .(ref[i]) == .(i)))
  }
  
#   # Traceplots for tau's.
#   # jagsfit$var takes three arguments: parameter, row, chain.
#   par(mfrow=c(2,1))
#   plot(jagsfit$tau1[1,,1], type="l", col="red", ylab="",
#        main=bquote(.(name) ~ .(ref[i]) == .(i)))
#   lines(jagsfit$tau1[1,,2], type="l", col="blue")
#   plot(density(jagsfit$tau1[1,,]), type="l",
#        main="tau1")
#   
#   plot(jagsfit$tau2[1,,1], type="l", col="red", ylab="",
#        main=bquote(.(name) ~ .(ref[i]) == .(i)))
#   lines(jagsfit$tau2[1,,2], type="l", col="blue")
#   plot(density(jagsfit$tau2[1,,]), type="l",
#        main="tau2")
#   
#   plot(jagsfit$tau3[1,,1], type="l", col="red", ylab="",
#        main=bquote(.(name) ~ .(ref[i]) == .(i)))
#   lines(jagsfit$tau3[1,,2], type="l", col="blue")
#   plot(density(jagsfit$tau3[1,,]), type="l",
#        main="tau3")
#   
#   plot(jagsfit$tau4[1,,1], type="l", col="red", ylab="",
#        main=bquote(.(name) ~ .(ref[i]) == .(i)))
#   lines(jagsfit$tau4[1,,2], type="l", col="blue")
#   plot(density(jagsfit$tau4[1,,]), type="l",
#        main="tau4")
  
  mean.params <- rowMeans(jagsfit$beta[,,1])
  return (mean.params)
}

PARAMS.BY.STORE <- Execute()

# Inspect results.
ref <- c("name", "log.a", "b", "g", "d")
par(mfrow=c(4, 1))
for (i in 2:5){
  hist(as.numeric(PARAMS.BY.STORE[,i]), breaks=30,
       main=bquote("Histogram of " ~ .(ref[i])))
}

