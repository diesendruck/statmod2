# StatMod2 - Polls

library(MASS)

#####
# SUMMARIZE.
# Albert and Chib (1993) present exact Bayesian methods for posterior sampling
# of binary response data. Where a probit regression would normally be used,
# they augment the model by introducing a left- or right-truncated normal
# variable Z, that corresponds to 1- or 0-valued y response, respectively. 
# In doing so, the model enables standard linear model results for normal-normal
# posteriors, which can be used to find full conditionals for Gibbs sampling.

#####
# LOAD DATA.
# Load and investigate data.
data <- read.csv("polls.csv")
attach(data)
#aggregate(data, by=list(state=state), FUN=mean)
#subset(data, is.na(bush) & state=="CA", select=c(bush, state))

#####
# CLEAN DATA.
# Remove NA's.
d <- subset(data, !is.na(bush))
d <- subset(d, select=c(bush, state, edu, age, female, black, weight))

# Add indicators for factor levels of edu and age.
unique.edus <- unique(d$edu)
unique.ages <- unique(d$age)
unique.states <- sort(unique(d$state))
for (i in 1:length(unique.edus)) {
  d[,c(paste(unique.edus[i]))] <- ifelse(d$edu==paste(unique.edus[i]), 1, 0)
}
for (i in 1:length(unique.ages)) {
  d[,c(paste(unique.ages[i]))] <- ifelse(d$age==paste(unique.ages[i]), 1, 0)
}
for (i in 1:length(unique.states)) {
  d[,c(paste(unique.states[i]))] <- ifelse(d$state==paste(unique.states[i]), 1, 0)
}

attach(d)
# Keep only numerical values for design matrix, X.
d <- d[,c("bush", "female", "black", "weight", "NoHS", "HS", "SomeColl", "Bacc",
          "18to29", "30to44", "45to64", "65plus", as.character(unique.states))]
d$edu.checksum <- ifelse(rowSums(d[,c(5:8)])>1, "1", "0")
d$age.checksum <- ifelse(rowSums(d[,c(9:12)])>1, "1", "0")
d$state.checksum <- ifelse(rowSums(d[,c(13:61)])>1, "1", "0")
if (sum(as.numeric(d$edu.checksum))>0) print("ERROR in Edu Indicators!")
if (sum(as.numeric(d$age.checksum))>0) print("ERROR in Age Indicators!")
if (sum(as.numeric(d$state.checksum))>0) print("ERROR in State Indicators!")

#####
# SET UP X AND Y MATRICES.
Y <- as.matrix(d$bush)
X <- as.matrix(d[,c(2:61)])

# OLS
f1 <- lm(Y~X)
f2 <- lm(data$bush ~ data$state+data$edu+data$age+data$female+data$black+data$weight)

#####
# DO GIBBS SAMPLER.

hist(rgamma(10000, n.i/(1+2*n.i), 1/2))

Gibbs <- function(data) {
  
  n.iter <- 2000
  
  B0 <- ginv(t(X)%*%(X))%*%t(X)%*%Y
  
  # Prepare per-state data.
  state.names <- names(d)[13:61]
  num.states <- length(state.names)
  
  # Create containers for chains and variables.
  PARAMS.BY.STORE <- array(rep(NA, num.stores*4*n.iter),
                           c(num.stores, 4, n.iter))
  MU <- matrix(NA, nrow=n.iter+1, ncol=4)
  MU[1,] <- c(mu0, mu1, mu2, mu3)
  TAUSQ <- matrix(NA, nrow=n.iter+1, ncol=4)
  TAUSQ[1,] <- c(t0sq, t1sq, t2sq, t3sq)
  
  # Do Gibbs many times.
  
  for (iter in 1:n.iter) {
    
    # Sample Beta's (be0, be1, be2, be3) for each store.
    # Use Normal-Normal full conditional posterior.
    for (i in 1:num.stores) {
      name <- toString(store.names[i])
      store.data <- data[which(data$store==name),]
      ni <- dim(store.data)[1]
      y <- store.data$log.Q
      X <- as.matrix(cbind(rep(1, ni),
                           store.data[,c("log.P", "disp", "cross")]))
      XtX <- t(X)%*%X
      Xty <- t(X)%*%y
      latest.mu <- MU[iter,]
      diag.tausqs <- diag(c(t0sq, t1sq, t2sq, t3sq))
      b <- sample.beta(sigma2, diag.tausqs, XtX, Xty, latest.mu)
      PARAMS.BY.STORE[i,,iter] <- b
    }
    
    # Use Normal-Normal full conditional posterior to update Mu's, given Beta's.
    beta.means <- colMeans(PARAMS.BY.STORE[,,iter])
    mu0 <- sample.mu(m0, v0, num.stores, beta.means[1], t0sq)
    mu1 <- sample.mu(m1, v1, num.stores, beta.means[2], t1sq)
    mu2 <- sample.mu(m2, v2, num.stores, beta.means[3], t2sq)
    mu3 <- sample.mu(m3, v3, num.stores, beta.means[4], t3sq)
    MU[iter+1,] <- c(mu0, mu1, mu2, mu3)
    
    # Use Normal-InvGamma full conditional posterior to update Tau's, given Beta's.
    t0sq <- sample.tausq(beta.means[1], MU[iter,][1], a0, b0)
    t1sq <- sample.tausq(beta.means[2], MU[iter,][2], a1, b1)
    t2sq <- sample.tausq(beta.means[3], MU[iter,][3], a2, b2)
    t3sq <- sample.tausq(beta.means[4], MU[iter,][4], a3, b3)
    TAUSQ[iter+1,] <- c(t0sq, t1sq, t2sq, t3sq)
  }  
  
  return (list(PARAMS.BY.STORE=PARAMS.BY.STORE, MU=MU, TAUSQ=TAUSQ))
}

sample.beta <- function(sigma2, diag.tausqs, XtX, Xty, latest.mu) {
  cov <- ginv(XtX/sigma2 + solve(diag.tausqs))
  mean <- cov%*%(Xty/sigma2 + solve(diag.tausqs)%*%latest.mu)
  beta <- mvrnorm(1, mu=mean, Sigma=cov)
  return (beta)
}
sample.mu <- function(mu.pr.mean, mu.pr.var, num.stores, beta.mean, beta.var) {
  var <- ginv(1/mu.pr.var + num.stores/beta.var)
  mean <- var*(mu.pr.mean/mu.pr.var + num.stores*beta.mean/beta.var)
  mu0 <- rnorm(1, mean=mean, sd=sqrt(var))
  return (mu0)
}
sample.tausq <- function(be0, mu0, pr.a, pr.b) {
  shape <- pr.a + 1/2
  rate <- 1/2*(be0-mu0)^2 + pr.b
  tausq <- rgamma(1, shape=shape, rate=rate)
  return (1/tausq)
}

results <- Gibbs(data)
P <- results$PARAMS.BY.STORE
MU <- results$MU
TAUSQ <- results$TAUSQ

# Show traceplots for Mu's and Tau's.
par(mfrow=c(2, 2))
its <- 1:dim(MU)[1]
its <- its[seq(1001, dim(MU)[1], 10)]
count <- dim(MU)[2]
for (p in 1:count) {
  plot(its, MU[,p][seq(1001, dim(MU)[1], 10)], xlab="Iteration",
       ylab=bquote("Value of Mu"~.(p)), type="l")
}
for (p in 1:count) {
  plot(its, TAUSQ[,p][seq(1001, dim(MU)[1], 10)], xlab="Iteration",
       ylab=bquote("Value of Tausq"~.(p)), type="l",
       ylim=c(0, 30))
}

# Show boxplots for betas of all stores together.
par(mfrow=c(2,2))
for (p in 1:count) {
  data.to.plot <- t(P[,p,seq(1001, dim(MU)[1]-1, 10)])
  m <- melt(data.to.plot)[c(2, 3)]
  names(m) <- c("store", "value")
  bymedian <- with(m, reorder(m$store, m$value, median))
  boxplot(m$value ~ bymedian, data=m, xlab="Store",
          ylab=bquote("Value of Beta"~.(p)),
          main=bquote("Beta"~.(p)~"Values by Store"))
}



