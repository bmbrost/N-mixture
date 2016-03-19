rm(list=ls())

###
### Simulate data according to model base.N.mixture.MCMC
###

m <- 1 # Number of sites
J <- rep(100,m) # Number of replicate observations
lambda <- rep(150,m) # Intensity of Poisson distribution
p <- rep(0.75,m) # Detection probability
N <- rpois(m,lambda) # True population size 
Y <- t(sapply(1:m,function(x) rbinom(J[x],N[x],p[x]))) # Observed counts


###
### Fit model using MCMC algorithm
###

source("/Users/brost/Documents/git/N-mixture/base.N.mixture.MCMC.R")
# hist(rgamma(1000,2,0.01))
priors <- list(a=1,b=1,r=2,q=0.01)
tune <- list(N=2)
start <- list(N=N,lambda=lambda,p=p)
out1 <- base.N.mixture.MCMC(Y,priors,tune,start,n.mcmc=50000)

matplot(out1$N,type="l");abline(h=N,col=2,lty=2)
matplot(out1$p,type="l");abline(h=p,col=2,lty=2)
matplot(out1$lambda,type="l");abline(h=lambda,col=2,lty=2)