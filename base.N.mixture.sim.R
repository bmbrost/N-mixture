###
### Simulate data according to model base.N.mixture.MCMC
###

rm(list=ls())

m <- 1 #Number of sites
J <- 100 #Number of replicate observations
lambda <- 150
p <- 0.75
N <- rpois(m,lambda)
Y <- matrix(0,m,J)
Y[] <- c(sapply(1:m,function(x) rbinom(J,N[x],p[x])))

library(unmarked)



source("base.N.mixture.MCMC.R")
hist(rgamma(1000,15,0.1))
out1 <- base.N.mixture.MCMC(Y,priors=list(p=c(1,1),lambda=c(15,0.1)),n.mcmc=1000)

plot(1:1000,out1$N,type="l");abline(h=N)
mean(out1$N[-(1:100)])
plot(1:1000,out1$lambda,type="l");abline(h=lambda)
mean(out1$lambda[-(1:100)])
plot(1:1000,out1$p,type="p");abline(h=p)
mean(out1$p[-(1:100)])













#####
##### Source MCMC function for N-mixture
#####

source("binom.beta.pois.mcmc.R")

#####
##### Fit model with R&D data 
#####

mcmc.out=binom.beta.pois.mcmc(c(15,11,12,5,12),50,5000)

#####
##### View Trace Plots 
#####

layout(matrix(1:6,3,2))
plot(mcmc.out$N.save[1,],type="l",main="",ylab=bquote(N[1]))
plot(mcmc.out$N.save[2,],type="l",main="",ylab=bquote(N[2]))
plot(mcmc.out$N.save[3,],type="l",main="",ylab=bquote(N[3]))
plot(mcmc.out$N.save[4,],type="l",main="",ylab=bquote(N[4]))
plot(mcmc.out$N.save[5,],type="l",main="",ylab=bquote(N[5]))
plot(mcmc.out$phi.save,type="l",main="",ylab=bquote(phi))

#####
##### View Posterior Results  
#####

mean(mcmc.out$phi.save[-(1:1000)])
quantile(mcmc.out$phi.save[-(1:1000)],c(0.025,0.975))

apply(mcmc.out$N.save[,-(1:1000)],1,mean)
apply(mcmc.out$N.save[,-(1:1000)],1,quantile,c(0.025,0.975))

layout(matrix(1:2,1,2))
hist(mcmc.out$phi.save[-(1:1000)],breaks=30,col=8,xlim=c(0,1),prob=TRUE,xlab=bquote(phi),ylab=substitute(paste("[ ",phi," | y ]")),main="")
hist(mcmc.out$N.total.save[-(1:1000)],breaks=1000,col=8,prob=TRUE,xlab="N",ylab="[ N | y ]",main="")