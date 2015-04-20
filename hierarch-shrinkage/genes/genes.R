# StatMod2 - Hierarchical Models and Shrinkage - Genes

# Get data.
library(lattice)
path <- paste("~/Google Drive/2. SPRING 2015/STAT MOD 2 - Prof Scott/",
              "hierarch-shrinkage/genes", sep="")
setwd(path)
data <- read.csv("droslong.csv")
data <- data[order(data[,1], data[,2]),]
head(data, 100)

# Quickly plot.
xyplot(log2exp~time | gene, data=data)
xyplot(log2exp~time | group, data=data)

# Average technical replicants, and run frequentist linear model fit for each,
# to get estimates for prior on group mean and group variance.
DATA.AVG.TR <- NULL
lm.coefs <- NULL
gene.names <- unique(data$gene)
for (name in gene.names) {
  # Subset of data for each gene.
  d <- data[which(data$gene==name),]
  gene.time.avgs <- aggregate(d$log2exp~d$time, d, mean)
  d <- cbind(d$gene[1], d$group[1], gene.time.avgs)
  names(d) <- c("gene", "group", "time", "log2exp")
  DATA.AVG.TR <- rbind(DATA.AVG.TR, d)
  
  # Also do a Frequentist exploration and store lm coefs.
  fit.coefs <- lm(d$log2exp~d$time)$coefficients
  lm.coefs <- rbind(lm.coefs, c(name, fit.coefs))
}
# Got linear model coefficients for each gene. Now average them for the prior
# on group mean and group variance.
c <- lm.coefs[,c(2,3)]
class(c)  <- "numeric"
lm.output <- colMeans(c)

D <- DATA.AVG.TR
xyplot(log2exp~time | gene, data=D)

sigma2 <- 1/rgamma(1, 0.5, 0.5)
a0=0.5; b0=0.5; # Hyper params for Tau's, which are ~ InvGamma.
a1=0.5; b1=0.5;
m0=lm.output[1]; v0=1000; # Hyper params for Mu's, which are ~ Normal.
m1=lm.output[2]; v1=1000;

Gibbs <- function(D, n.iter=1000) {
  t0sq <- 1/rgamma(1, a0, b0) # Initialize Tau's.
  t1sq <- 1/rgamma(1, a1, b1)
  mu0 <- rnorm(1, m0, v0) # Initialize Mu's.
  mu1 <- rnorm(1, m1, v1)
  
  # Prepare per-store data.
  gene.names <- unique(D$gene)
  num.genes <- length(gene.names)
  
  # Create containers for chains and variables.
  PARAMS.BY.GENE <- array(rep(NA, num.genes*2*n.iter),
                           c(num.genes, 2, n.iter))
  MU <- matrix(NA, nrow=n.iter+1, ncol=2)
  MU[1,] <- c(mu0, mu1)
  TAUSQ <- matrix(NA, nrow=n.iter+1, ncol=2)
  TAUSQ[1,] <- c(t0sq, t1sq)
  
  # Do Gibbs many times.
  
  for (iter in 1:n.iter) {
    
    # Sample Beta's (be0, be1, be2, be3) for each store.
    # Use Normal-Normal full conditional posterior.
    for (i in 1:num.genes) {
      name <- toString(gene.names[i])
      gene.data <- data[which(D$gene==name),]
      ni <- dim(gene.data)[1]
      y <- gene.data$log2exp
      X <- as.matrix(cbind(rep(1, ni),
                           gene.data[,c("time")]))
      XtX <- t(X)%*%X
      Xty <- t(X)%*%y
      latest.mu <- MU[iter,]
      diag.tausqs <- diag(c(t0sq, t1sq))
      b <- sample.beta(sigma2, diag.tausqs, XtX, Xty, latest.mu)
      PARAMS.BY.GENE[i,,iter] <- b
    }
    
    # Use Normal-Normal full conditional posterior to update Mu's, given Beta's.
    beta.means <- colMeans(PARAMS.BY.GENE[,,iter])
    mu0 <- sample.mu(m0, v0, num.genes, beta.means[1], t0sq)
    mu1 <- sample.mu(m1, v1, num.genes, beta.means[2], t1sq)
    MU[iter+1,] <- c(mu0, mu1)
    
    # Use Normal-InvGamma full conditional posterior to update Tau's, given Beta's.
    t0sq <- sample.tausq(beta.means[1], MU[iter,][1], a0, b0)
    t1sq <- sample.tausq(beta.means[2], MU[iter,][2], a1, b1)
    TAUSQ[iter+1,] <- c(t0sq, t1sq)
  }  
  
  return (list(PARAMS.BY.GENE=PARAMS.BY.GENE, MU=MU, TAUSQ=TAUSQ))
}

sample.beta <- function(sigma2, diag.tausqs, XtX, Xty, latest.mu) {
  cov <- ginv(XtX/sigma2 + solve(diag.tausqs))
  mean <- cov%*%(Xty/sigma2 + solve(diag.tausqs)%*%latest.mu)
  beta <- mvrnorm(1, mu=mean, Sigma=cov)
  return (beta)
}
sample.mu <- function(mu.pr.mean, mu.pr.var, num.genes, beta.mean, beta.var) {
  var <- ginv(1/mu.pr.var + num.genes/beta.var)
  mean <- var*(mu.pr.mean/mu.pr.var + num.genes*beta.mean/beta.var)
  mu0 <- rnorm(1, mean=mean, sd=sqrt(var))
  return (mu0)
}
sample.tausq <- function(be0, mu0, pr.a, pr.b) {
  shape <- pr.a + 1/2
  rate <- 1/2*(be0-mu0)^2 + pr.b
  tausq <- rgamma(1, shape=shape, rate=rate)
  return (1/tausq)
}


#############################################

results <- Gibbs(data)
P <- results$PARAMS.BY.GENE
MU <- results$MU
TAUSQ <- results$TAUSQ
gene.names <- unique(D$gene)

# Show traceplots for Mu's and Tau's.
par(mfrow=c(2, 1))
its <- 1:dim(MU)[1]
# Burn first 10th and thin by 5.
it <- seq(floor(0.1*dim(MU)[1]), dim(MU)[1], 5)
its <- its[it]
count <- dim(MU)[2]
for (p in 1:count) {
  plot(its, MU[,p][it], xlab="Iteration",
       ylab=bquote("Value of Mu"~.(p)), type="l")
}
for (p in 1:count) {
  plot(its, TAUSQ[,p][it], xlab="Iteration",
       ylab=bquote("Value of Tausq"~.(p)), type="l")
}

# Show boxplots for betas of all stores together.
par(mfrow=c(2,1))
for (p in 1:count) {
  data.to.plot <- t(P[,p,it])
  m <- melt(data.to.plot)[c(2, 3)]
  names(m) <- c("gene", "value")
  bymedian <- with(m, reorder(m$gene, m$value, median))
  boxplot(m$value ~ bymedian, data=m, xlab="Gene",
          ylab=bquote("Value of Beta"~.(p)),
          main=bquote("Beta"~.(p-1)~"Values by Gene"),
          yaxt="n")
  if (p==1) {
    axis(2, at=seq(6, 16, by=2), las=2)
  } else if (p==2) {
    axis(2, at=seq(-3, 3, by=0.2), las=2)
  }
}






