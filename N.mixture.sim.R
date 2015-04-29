rm(list=ls())

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-p))

###
### Simulate data according to model N.mixture.MCMC
###

m <- 3 # Number of sites
J <- 20 # Number of replicate observations
lambda <- rep(150,m) # Intensity of Poisson distribution
lambda <- rpois(m,100)
N <- rpois(m,lambda) # True population size 

W <- array(1,dim=c(J,2,m)) # Design matrix for detection probability
# hist(rgamma(1000,12,2))
W[,2,] <- rgamma(m*J,12,2) # Detection covariate, e.g., flight duration
hist(W[,2,])

alpha <- matrix(c(-5,0.75),2,1) # Detection parameters
curve(expit(cbind(1,x)%*%alpha),xlim=c(0,max(W[,2,])),add=FALSE)
p <- t(apply(W,3,function(x) expit(x%*%alpha))) # Detection probability
Y <- t(sapply(1:m,function(x) rbinom(J,N[x],p[x,]))) # Observed counts

# Compare observed versus true detection probability
p.hat <- t(sapply(1:m,function(x) Y[x,]/N[x])) # Observed detection probability
matplot(p.hat,pch=19,cex=1,col=rgb(0,0,0,0.25),xlab="Site")
matplot(p,pch=19,cex=0.5,col=rgb(1,0,0,0.25),add=TRUE) # True detection probability

# Plot true true abundance and observed counts
plot(N,pch=19,ylim=c(min(Y),max(N))) # True abundance
matplot(Y,pch=19,cex=0.5,col=rgb(1,0,0,0.25),add=TRUE) # Observed counts


###
### Fit model using MCMC algorithm
###

W.mat <- apply(W,2,I) # Convert W from 3-D array to 2-D matrix

source("N.mixture.MCMC.R")
# hist(rgamma(1000,5,0.01))
priors <- list(r=5,q=0.001,tau=2)
tune <- list(N=10,alpha=0.01)
start <- list(N=N,alpha=alpha,lambda=lambda)
out1 <- N.mixture.MCMC(Y,W.mat,priors,tune,start,n.mcmc=3000)

matplot(out1$N,type="l",col=1:m,lty=1);abline(h=N,col=1:m,lty=2)
apply(out1$N,2,mean);apply(out1$N,2,quantile,c(0.025,0.975))
N                           
matplot(out1$alpha,type="l",col=1:2);abline(h=alpha,col=1:2,lty=2)
apply(out1$alpha,2,mean)
matplot(out1$lambda,type="l",col=1:m,lty=1);abline(h=lambda,col=1:m,lty=2)
apply(out1$lambda,2,mean);apply(out1$lambda,2,quantile,c(0.025,0.975))
lambda


###
### Fit model using R package 'unmarked'
###

# library(unmarked)
# visitMat <- matrix(as.character(1:J), m, J, byrow=TRUE)
# umf <- unmarkedFramePCount(y=Y, obsCovs=list(visit=visitMat))
# summary(umf)
# 
# fm1 <- pcount(~1 ~1, umf, K=200)
# fm1
# plogis(coef(fm1, type="det")) # Should be close to p
