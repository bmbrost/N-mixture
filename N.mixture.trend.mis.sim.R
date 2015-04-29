rm(list=ls())

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-p))

###
### Simulate data according to model N.mixture.trend.MCMC
###

T <- 10 # Number of years of data
t <- 1:T
m <- 4 # Number of replicate observations
lambda <- 366 # Intensity of Poisson distribution for t=1
theta <- 1-0.4/100 # Population growth rate of -0.4%/year; can use for identity link for Poisson rate

# Simulate true abundance with trend
N <- numeric(T)
N[1] <- rpois(1,lambda) # Initial population size at t=1
# for(i in 2:T) N[i] <- rpois(1,theta*N[i-1]) # Intensity with linear trend; indentity link
for(i in 2:T) N[i] <- rpois(1,exp(log(theta)+log(N[i-1]))) # Intensity with linear trend; log link
plot(N,type="l")

W <- array(1,dim=c(m,2,T)) 
# hist(rgamma(1000,18,3))
W[,2,] <- rgamma(m*T,18,3) # Detection covariate, e.g., flight duration
hist(W[,2,])

alpha <- matrix(c(-5,0.7),2,1) # Detection parameters
# alpha <- matrix(c(-4,0.1),2,1) # Detection parameters
curve(expit(cbind(1,x)%*%alpha),xlim=c(0,max(W[,2,])),add=FALSE) # Response curve for p
p <- t(apply(W,3,function(x) expit(x%*%alpha))) # Detection probability
Y <- t(sapply(1:T,function(x) rbinom(m,N[x],p[x,]))) # Observed counts, T*m matrix

# Compare observed versus true detection probability
p.hat <- t(sapply(1:T,function(x) Y[x,]/N[x])) # Observed detection probability
matplot(p.hat,pch=19,cex=1,col=rgb(0,0,0,0.25),xlab="Time")
matplot(p,pch=19,cex=0.5,col=rgb(1,0,0,0.25),add=TRUE) # True detection probability

# Plot true true abundance and observed counts
plot(N,type="b",ylim=c(min(Y),max(N)),pch=19) # True abundance
points(rep(1:T,m),c(Y),pch=19,cex=(c(p)+(1-min(p)))^2,col=rgb(1,0,0,0.25)) # Observed counts

# Remove observations to test posterior predictive distribution for missing data
mis.idx <- 9 # Row of observations to remove
t <- t[-mis.idx]
W <- W[,,-mis.idx]
Y <- Y[-mis.idx,]

plot(N,type="b",ylim=c(min(Y),max(N)),pch=19) # True abundance
points(rep(t,m),c(Y),pch=19,cex=(c(p)+(1-min(p)))^2,col=rgb(1,0,0,0.25)) # Observed counts


###
### Fit model using MCMC algorithm
###

source("/Users/brost/Documents/git/Nmixture/N.mixture.trend.mis.MCMC.R")
# hist(rgamma(1000,9,0.025),breaks=100)

W.mat <- apply(W,2,I) # Convert W from 3-D array to 2-D matrix

priors <- list(r=9,q=0.025,tau=2,sigma=0.1)
tune <- list(N=40,alpha=0.003,theta=0.05)
# start <- list(N=N,alpha=alpha,lambda=lambda,theta=theta)
start <- list(N=round(apply(Y,1,max)*1.1),alpha=c(-5,1),lambda=350,theta=1)
out1 <- N.mixture.trend.mis.MCMC(t,Y,W.mat,priors,tune,start,n.mcmc=5000)

matplot(out1$alpha[,],type="l");abline(h=alpha,col=2,lty=2)
apply(out1$alpha,2,mean);apply(out1$alpha,2,quantile,c(0.025,0.975))
matplot(out1$N[-c(1:1000),],type="l",col=1:T,lty=1);abline(h=N,col=1:T,lty=2)
apply(out1$N,2,mean);apply(out1$N,2,quantile,c(0.025,0.975))
N
plot(out1$lambda,type="l");abline(h=lambda,col=2,lty=2)
plot(out1$theta,type="l");abline(h=theta,col=2,lty=2)
quantile(out1$theta,c(0.025,0.975))

idx <- 5
hist(out1$N[-(1:1000),idx],col="gray",breaks=100)
abline(v=N[idx],lty=2,col=2)

