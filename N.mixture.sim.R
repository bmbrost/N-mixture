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
p <- apply(W,3,function(x) expit(x%*%alpha)) # Detection probability
Y <- t(sapply(1:m,function(x) rbinom(J,N[x],p[,x]))) # Observed counts

###
### Fit model using MCMC algorithm
###

source("N.mixture.MCMC.R")
# hist(rgamma(1000,5,0.01))
out1 <- N.mixture.MCMC(Y,W,priors=list(r=5,q=0.001,tau=2),tune=list(N=5,alpha=0.01),n.mcmc=10000)

matplot(t(out1$N),type="l");abline(h=N,col=2,lty=2)
matplot(out1$alpha,type="l");abline(h=alpha,col=2,lty=2)
matplot(t(out1$lambda),type="l");abline(h=lambda,col=2,lty=2)

mean(out1$N[-(1:100)])
mean(out1$lambda[-(1:100)])
mean(out1$p[-(1:100)])


###
### Fit model using R package 'unmarked'
###

library(unmarked)
visitMat <- matrix(as.character(1:J), m, J, byrow=TRUE)
umf <- unmarkedFramePCount(y=Y, obsCovs=list(visit=visitMat))
summary(umf)

fm1 <- pcount(~1 ~1, umf, K=200)
fm1
plogis(coef(fm1, type="det")) # Should be close to p
