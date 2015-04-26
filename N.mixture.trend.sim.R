rm(list=ls())

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-p))

###
### Simulate data according to model N.mixture.trend.MCMC
###

T <- 8 # Number of years of data
m <- 10 # Number of replicate observations
lambda <- 1200 # Intensity of Poisson distribution for t=1
theta <- 0.987 # Population growth rate
# z <- log(lambda)
  
# Simulate true abundance with trend
N <- numeric(T)
N[1] <- rpois(1,lambda) # Initial population size at t=1
for(i in 2:T) N[i] <- rpois(1,theta*N[i-1]) # Intensity with linear trend
plot(N,type="l")

W <- array(1,dim=c(m,2,T)) # Design matrix for detection probability
# hist(rgamma(1000,12,2))
W[,2,] <- rgamma(m*T,12,2) # Detection covariate, e.g., flight duration
hist(W[,2,])

alpha <- matrix(c(-5,0.75),2,1) # Detection parameters
curve(expit(cbind(1,x)%*%alpha),xlim=c(0,max(W[,2,])),add=FALSE) # Response curve for p
p <- t(apply(W,3,function(x) expit(x%*%alpha))) # Detection probability
Y <- t(sapply(1:T,function(x) rbinom(m,N[x],p[x,]))) # Observed counts, T*m matrix

# Compare observed versus true detection probability
p.hat <- t(sapply(1:T,function(x) Y[x,]/N[x])) # Observed detection probability
matplot(p.hat,pch=19,cex=1,col=rgb(0,0,0,0.25),xlab="Time")
matplot(p,pch=19,cex=0.5,col=rgb(1,0,0,0.25),add=TRUE) # True detection probability

# Plot true true abundance and observed counts
plot(N,type="l",ylim=c(min(Y),max(N))) # True abundance
matplot(Y,pch=19,cex=0.5,col=rgb(1,0,0,0.25),add=TRUE) # Observed counts

###
### Fit model using MCMC algorithm
###

source("N.mixture.trend.MCMC.R")
# hist(rgamma(1000,2,0.001),breaks=100)

priors <- list(r=2,q=0.001,tau=2,sigma=2)
tune <- list(N=30,alpha=0.003,theta=0.05)
start <- list(N=N,alpha=alpha,lambda=lambda,theta=theta)
out1 <- N.mixture.trend.MCMC(Y,W,priors,tune,start,n.mcmc=3000)

matplot(out1$alpha[,],type="l");abline(h=alpha,col=2,lty=2)
apply(out1$alpha,2,mean);apply(out1$alpha,2,quantile,c(0.025,0.975))
matplot(out1$N,type="l",col=1:T,lty=1);abline(h=N,col=1:T,lty=2)
apply(out1$N,2,mean)
N
plot(out1$lambda,type="l");abline(h=lambda,col=2,lty=2)
plot(out1$theta,type="l");abline(h=theta,col=2,lty=2)
