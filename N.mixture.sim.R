rm(list=ls())

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-p))

###
### Simulate data according to model N.mixture.MCMC
###

m <- 2 # Number of sites
J <- rpois(m,10) # Number of replicate observations per site; irregular sampling effort
lambda <- rpois(m,100) # Intensity of Poisson distribution
N <- rpois(m,lambda) # True population size 

# Create design matrix; pad with NAs to account for irregular sampling effort
W <- array(NA,dim=c(max(J),2,m))
for(i in 1:m){
  W[1:J[i],,i] <- cbind(1,rgamma(J[i],12,2)) # Detection covariate, e.g., flight duration
}
W <- apply(W,2,I) # Flatten array; i.e., convert 3D array to 2D matrix
hist(W[,2])

# Detection probability
alpha <- matrix(c(-5,0.75),2,1) # Detection parameters
curve(expit(cbind(1,x)%*%alpha),xlim=c(0,max(W[,2],na.rm=TRUE)),add=FALSE)
p <- matrix(expit(W%*%alpha),m,max(J),byrow=TRUE) # Detection probability

# Observed counts; pad with NAs to account for irregular sampling effort
Y <- matrix(NA,m,max(J))
for(i in 1:m){
  Y[i,1:J[i]] <- rbinom(J[i],N[i],p[i,])  
}

# Compare observed versus true detection probability
p.hat <- Y/N # Observed detection probability
matplot(p.hat,pch=19,cex=1,col=rgb(0,0,0,0.25),xlab="Site")
matplot(p,pch=19,cex=0.5,col=rgb(1,0,0,0.25),add=TRUE) # True detection probability

# Plot true abundance and observed counts
plot(N,pch=19,ylim=c(min(Y,na.rm=TRUE),max(N))) # True abundance
points(rep(1:m,max(J)),Y,cex=(p+1)^2,pch=19,col=rgb(1,0,0,0.25)) # Observed counts


###
### Fit model using MCMC algorithm
###

source("/Users/brost/Documents/git/N-mixture/N.mixture.MCMC.R")
# hist(rgamma(1000,5,0.01))
priors <- list(r=5,q=0.001,tau=2)
tune <- list(N=10,alpha=0.01)
start <- list(N=N,alpha=alpha,lambda=lambda)
out1 <- N.mixture.MCMC(Y,W,priors,tune,start,n.mcmc=3000)

matplot(out1$N,type="l",col=1:m,lty=1);abline(h=N,col=1:m,lty=2)
apply(out1$N,2,quantile,c(0.025,0.975))
N
matplot(out1$alpha,type="l",col=1:2);abline(h=alpha,col=1:2,lty=2)
apply(out1$alpha,2,quantile,c(0.025,0.975))
alpha
matplot(out1$lambda,type="l",col=1:m,lty=1);abline(h=lambda,col=1:m,lty=2)
apply(out1$lambda,2,quantile,c(0.025,0.975))
lambda
