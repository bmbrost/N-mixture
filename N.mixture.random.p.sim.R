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
zeta <- 1
eps <- matrix(NA,m,max(J)) # Random effect on p
for(i in 1:m) eps[i,1:J[i]] <- rnorm(J[i],0,zeta)
alpha <- matrix(c(-5,0.75),2,1) # Detection parameters
curve(expit(cbind(1,x)%*%alpha),xlim=c(0,max(W[,2,],na.rm=TRUE)),add=FALSE,col=2,lwd=2)
for(i in 1:(max(J)*m)) curve(expit(cbind(1,x)%*%alpha+c(eps)[i]),add=TRUE,lty=1,col=rgb(0,0,0,0.2))

p <- matrix(NA,m,max(J),byrow=TRUE) # Detection probability
for(i in 1:m){
  p[i,] <- expit(W[,,i]%*%alpha+eps[i,])
}

# Observed counts; pad with NAs to account for irregular sampling effort
Y <- matrix(NA,m,max(J))
for(i in 1:m){
  Y[i,1:J[i]] <- rbinom(J[i],N[i],p[i,])  
}

W <- apply(W,2,I) # Flatten array; i.e., convert 3D array to 2D matrix
eps <- c(t(eps)) # Vectorize eps

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

source("/Users/brost/Documents/git/N-mixture/N.mixture.random.p.MCMC.R")
# hist(rgamma(1000,5,0.01))
priors <- list(r=5,q=0.001,tau=2,a=0,b=3)
tune <- list(N=5,alpha=0.005,eps=0.75,zeta=0.05)
start <- list(N=N,alpha=alpha,lambda=lambda,zeta=zeta,eps=eps)
out1 <- N.mixture.random.p.MCMC(Y,W,priors,tune,start,n.mcmc=10000)

matplot(out1$N,type="l",col=1:m,lty=1);abline(h=N,col=1:m,lty=2)
apply(out1$N,2,quantile,c(0.025,0.975))
N
matplot(out1$alpha,type="l",col=1:2);abline(h=alpha,col=1:2,lty=2)
matplot(out1$lambda,type="l",col=1:m,lty=1);abline(h=lambda,col=1:m,lty=2)
apply(out1$lambda,2,quantile,c(0.025,0.975))
plot(out1$zeta,type="l");abline(h=zeta,col=2,lty=2)
idx <- 50
plot(out1$eps[,idx],type="l");abline(h=eps[idx],col=2,lty=2)
