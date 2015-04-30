rm(list=ls())

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-p))

###
### Simulate data according to model N.mixture.trend.MCMC
###

T <- 10 # Number of time periods
m <- rpois(T,10) # Number of replicate observations per time period; irregular sampling effort
lambda <- 366 # Intensity of Poisson distribution for t=1
theta <- 1-0.4/100 # Population growth rate of -0.4%/year; can use for identity link for Poisson rate

# Simulate true abundance with trend
N <- numeric(T)
N[1] <- rpois(1,lambda) # Initial population size at t=1
# for(i in 2:T) N[i] <- rpois(1,theta*N[i-1]) # Intensity with linear trend; indentity link
for(i in 2:T) N[i] <- rpois(1,exp(log(theta)+log(N[i-1]))) # Intensity with linear trend; log link
plot(N,type="l")

W <- array(NA,dim=c(max(m),2,T)) 
# hist(rgamma(1000,18,3))
for(i in 1:T){
  W[1:m[i],,i] <- cbind(1,rgamma(m[i],18,3)) # Detection covariate, e.g., flight duration
}
W <- apply(W,2,I) # Flatten array; i.e., convert 3D array to 2D matrix
hist(W[,2])

# Detection probability
alpha <- matrix(c(-5,0.70),2,1) # Detection parameters
curve(expit(cbind(1,x)%*%alpha),xlim=c(0,max(W[,2],na.rm=TRUE)),add=FALSE)
p <- matrix(expit(W%*%alpha),T,max(m),byrow=TRUE) # Detection probability

# Observed counts; pad with NAs to account for irregular sampling effort
Y <- matrix(NA,T,max(m))
for(i in 1:T){
  Y[i,1:m[i]] <- rbinom(m[i],N[i],p[i,])  
}

# Compare observed versus true detection probability
p.hat <- Y/N # Observed detection probability
matplot(p.hat,pch=19,cex=1,col=rgb(0,0,0,0.25),xlab="Site")
matplot(p,pch=19,cex=0.5,col=rgb(1,0,0,0.25),add=TRUE) # True detection probability

# Plot true abundance and observed counts
plot(N,pch=19,ylim=c(min(Y,na.rm=TRUE),max(N))) # True abundance
points(rep(1:T,max(m)),Y,cex=(p+1)^2,pch=19,col=rgb(1,0,0,0.25)) # Observed counts

###
### Simulate missing data
###

idx <- 9 # Row of Y to remove
Y[idx,] <- NA
W[((idx-1)*max(m)+1):((idx)*max(m)),] <- NA # Remove corresponding W


###
### Fit model using MCMC algorithm
###

source("/Users/brost/Documents/git/Nmixture/N.mixture.trend.MCMC.R")
# hist(rgamma(1000,9,0.025),breaks=100)
priors <- list(r=9,q=0.025,tau=2,sigma=0.1)
tune <- list(N=40,alpha=0.003,theta=0.05)
# start <- list(N=N,alpha=alpha,lambda=lambda,theta=theta)
N.start <- round(sapply(1:T,function(x) # Moving average for missing time periods
  ifelse(sum(Y[x,],na.rm=TRUE)>0,max(Y[x,],na.rm=TRUE), mean(c(Y[x-1,],Y[x+1,]),na.rm=TRUE)))*1.1)
start <- list(N=N.start,alpha=c(-5,1),lambda=350,theta=1)
out1 <- N.mixture.trend.MCMC(Y,W,priors,tune,start,n.mcmc=5000)

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

