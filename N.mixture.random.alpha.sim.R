### Note that the corresponding model, with alpha (i.e., coefficients that quantify the affect of covariates on
### detection variable, alpha[i]~N(mu.alpha,tau^2*I)) that varies by sites, is a work progress. 
### It currently does not estimate mu.alpha properly.



rm(list=ls())

library(mvtnorm)
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-p))

###
### Simulate data according to model N.mixture.re.MCMC
###

m <- 5 # Number of sites
J <- rpois(m,30) # Number of replicate observations per site; irregular sampling effort
lambda <- rpois(m,300) # Intensity of Poisson distribution
N <- rpois(m,lambda) # True population size 

# Create design matrix; pad with NAs to account for irregular sampling effort
W <- array(NA,dim=c(max(J),2,m))
for(i in 1:m){
  W[1:J[i],,i] <- cbind(1,rgamma(J[i],12,2)) # Detection covariate, e.g., flight duration
}
hist(W[,2,])


# Detection probability
mu.alpha <- matrix(c(-5,0.75),2,1) # Detection parameters
tau <- 0.25
alpha <- rmvnorm(m,mu.alpha,tau^2*diag(length(mu.alpha)))
curve(expit(cbind(1,x)%*%mu.alpha),xlim=c(0,max(W[,2,],na.rm=TRUE)),add=FALSE)
for(i in 1:m) curve(expit(cbind(1,x)%*%alpha[i,]),add=TRUE,lty=2,col=i)

p <- matrix(NA,m,max(J),byrow=TRUE) # Detection probability
for(i in 1:m){
  p[i,] <- expit(W[,,i]%*%alpha[i,])
}

# Observed counts; pad with NAs to account for irregular sampling effort
Y <- matrix(NA,m,max(J))
for(i in 1:m){
  Y[i,1:J[i]] <- rbinom(J[i],N[i],p[i,])  
}

W <- apply(W,2,I) # Flatten array; i.e., convert 3D array to 2D matrix

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

source("/Users/bmb/Documents/git/Nmixture/N.mixture.random.alpha.MCMC.R")
# hist(rgamma(1000,5,0.01))
priors <- list(r=5,q=0.001,zeta=2,a=0,b=5,mu.0=rep(0,2))
tune <- list(N=15,alpha=0.01,tau=0.25)
start <- list(N=N,alpha=alpha,lambda=lambda,mu.alpha=mu.alpha,tau=tau)
out1 <- N.mixture.random.alpha.MCMC(Y,W,priors,tune,start,n.mcmc=10000)

matplot(out1$N,type="l",col=1:m,lty=1);abline(h=N,col=1:m,lty=2)
apply(out1$N,2,quantile,c(0.025,0.975))
N
matplot(out1$mu.alpha,type="l",col=1:2);abline(h=mu.alpha,col=1:2,lty=2)
idx <- 1
matplot(out1$alpha[,,idx],type="l",col=1:2);abline(h=alpha[idx,],col=1:2,lty=2)
plot(out1$tau,type="l");abline(h=tau,col=2,lty=2)
matplot(out1$lambda,type="l",col=1:m,lty=1);abline(h=lambda,col=1:m,lty=2)
apply(out1$lambda,2,quantile,c(0.025,0.975))
lambda
