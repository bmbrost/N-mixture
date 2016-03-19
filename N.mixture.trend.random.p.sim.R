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
hist(W[,2,])

# Detection probability
zeta <- 1
eps <- matrix(NA,T,max(m)) # Random effect on p
for(i in 1:T) eps[i,1:m[i]] <- rnorm(m[i],0,zeta)
alpha <- matrix(c(-5,0.75),2,1) # Detection parameters
curve(expit(cbind(1,x)%*%alpha),xlim=c(0,max(W[,2,],na.rm=TRUE)),add=FALSE,col=2,lwd=2)
for(i in 1:(max(m)*T)) curve(expit(cbind(1,x)%*%alpha+c(eps)[i]),add=TRUE,lty=1,col=rgb(0,0,0,0.2))

p <- matrix(NA,T,max(m),byrow=TRUE) # Detection probability
for(i in 1:T){
  p[i,] <- expit(W[,,i]%*%alpha+eps[i,])
}

# Observed counts; pad with NAs to account for irregular sampling effort
Y <- matrix(NA,T,max(m))
for(i in 1:T){
  Y[i,1:m[i]] <- rbinom(m[i],N[i],p[i,])  
}

W <- apply(W,2,I) # Flatten array; i.e., convert 3D array to 2D matrix
eps <- c(t(eps)) # Vectorize eps

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

source("/Users/brost/Documents/git/N-mixture/N.mixture.trend.random.p.MCMC.R")
# hist(rgamma(1000,9,0.025),breaks=100)
priors <- list(r=9,q=0.025,tau=5,sigma=0.1,a=0,b=3)
tune <- list(N=40,alpha=0.003,theta=0.05,zeta=0.25,eps=0.5)
# start <- list(N=N,alpha=alpha,lambda=lambda,theta=theta)
N.start <- round(sapply(1:T,function(x) # Moving average for missing time periods
  ifelse(sum(Y[x,],na.rm=TRUE)>0,max(Y[x,],na.rm=TRUE), mean(c(Y[x-1,],Y[x+1,]),na.rm=TRUE)))*1.1)
start <- list(N=N.start,alpha=c(-5,1),lambda=350,theta=1,eps=eps,zeta=zeta)
out1 <- N.mixture.trend.random.p.MCMC(Y,W,priors,tune,start,n.mcmc=50000)

idx <- 40000:50000
matplot(out1$alpha[idx,],type="l");abline(h=alpha,col=2,lty=2)
apply(out1$alpha[idx,],2,mean);apply(out1$alpha[idx,],2,quantile,c(0.025,0.975))
matplot(out1$N[idx,],type="l",col=1:T,lty=1);abline(h=N,col=1:T,lty=2)
apply(out1$N[idx,],2,mean);apply(out1$N[idx,],2,quantile,c(0.025,0.975))
N
plot(out1$lambda[idx],type="l");abline(h=lambda,col=2,lty=2)
plot(out1$theta[idx],type="l");abline(h=theta,col=2,lty=2)
quantile(out1$theta[idx],c(0.025,0.975))
plot(out1$zeta[idx],type="l");abline(h=zeta,col=2,lty=2)

idx <- 1
plot(out1$eps[,idx],type="l");abline(h=eps[idx],col=2,lty=2)

idx <- 5
hist(out1$N[-(1:1000),idx],col="gray",breaks=100)
abline(v=N[idx],lty=2,col=2)

beanplot(c(out1$N[idx,])~rep(1:T,each=length(idx)),ll=0,beanlinewd=0,bw=10,log="",
   col=rgb(0,0,1,0.25),las=1,add=FALSE,ylim=c(0,max(out1$N[idx,])*1.1))
points(N,pch=19)
points(rep(1:T,max(m)),Y,cex=(p+1)^2,pch=19,col=rgb(1,0,0,0.25)) # Observed counts
