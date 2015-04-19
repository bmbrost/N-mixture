rm(list=ls())

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-p))

###
### Simulate data according to model N.mixture.MCMC
###

m <- 1 # Number of sites
T <- 20 # Number of years of data
J <- 20 # Number of replicate observations
lambda <- numeric(T)
lambda[1] <- 1200 # Intensity of Poisson distribution for t=1
theta <- 0.987 # Population growth rate
z <- log(lambda)

for(i in 2:T) lambda[i] <- theta*lambda[i-1] # Intensity with linear trend
plot(lambda,type="l")
for(i in 2:T) z[i] <- log(theta)+z[i-1] # Log-transformed intensity with linear trend
plot(z,type="l")

N <- rpois(T,exp(z)) # True population size
plot(N,type="l")

W <- array(1,dim=c(J,2,T)) # Design matrix for detection probability
# hist(rgamma(1000,12,2))
W[,2,] <- rgamma(J*T,12,2) # Detection covariate, e.g., flight duration
hist(W[,2,])

alpha <- matrix(c(-5,0.75),2,1) # Detection parameters
curve(expit(cbind(1,x)%*%alpha),xlim=c(0,max(W[,2,])),add=FALSE)
p <- apply(W,3,function(x) expit(x%*%alpha)) # Detection probability
Y <- t(sapply(1:T,function(x) rbinom(J,N[x],p[,x]))) # Observed counts, T*J matrix


###
### Fit model using MCMC algorithm
###

source("N.mixture.trend.MCMC.R")
# hist(rgamma(1000,5,0.01))

priors <- list(r=5,q=0.001,tau=2,sigma=2)
tune <- list(N=5,alpha=0.01,theta=0.01)
out1 <- N.mixture.trend.MCMC(Y,W,priors=priors,tune=tune,n.mcmc=10000)

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
