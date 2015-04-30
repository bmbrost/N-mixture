### Note that this model, with alpha (i.e., coefficients that quantify the affect of covariates on
### detection variable, alpha[i]~N(mu.alpha,tau^2*I)) that varies by sites, is a work progress. 
### It currently does not estimate mu.alpha properly.






N.mixture.random.alpha.MCMC <- function(Y,W,priors,tune,start,n.mcmc=1000){

  ###
  ### Brian M. Brost (25APR2015)
  ###
  ### N-mixture model for multiple sites with detection modeled by covariates and random effect on alpha
  ###
  ### Model statement: (i indexes site, j indexes observation)
  ### Y[i,j]~Binom(N[i],p[i,j])
  ### N[i]~Pois(lambda[i])
  ### lambda[i]~Gamma(r,q)
  ### logit(p[i,j])~W[j,,i]%*%alpha
  ### alpha[t] ~ N(mu_alpha,tau^2*I)
  ### mu_alpha~N(0,zeta^2*I)
  ###
  ### Function arguments: 
  ### Y=m*J matrix, where m is the number of sites and J is the maximum number of
  ###   observations across all sites
  ### W=J*length(alpha)*m array; each 'slice' of W is a design matrix for detection model at a site
  ### priors=parameters of prior distributions for lambda and alpha
  ### tune=tuning parameters for N and alpha
  ### start=starting values for N,lambda, and alpha
  ###
    
  ###
  ###  Libraries and subroutines
  ###
  
  expit <- function(x) 1/(1+exp(-x))
  logit <- function(x) log(x/(1-p))
  
  
  ###
  ###  Setup Variables 
  ###

# browser()
  m <- nrow(Y) # Number of sites
  qY <- ncol(Y)
  J <- apply(Y,1,function(x) sum(!is.na(x))) # Number of observations per site
  y.max <- apply(Y,1,max,na.rm=TRUE)
  qW <- ifelse(is.null(ncol(W)),1,ncol(W))
  keep <- list(N=0,alpha=0,tau=0)

  N.save <- matrix(0,n.mcmc,m)
  lambda.save <- matrix(0,n.mcmc,m)
  mu.alpha.save <- matrix(0,n.mcmc,qW)
  alpha.save <- array(0,dim=c(n.mcmc,qW,m))  
  tau.save <- numeric(n.mcmc)  

  ###
  ###  Tuning and starting values 
  ###
  
  N.tune <- seq(-1*tune$N,tune$N,1)

  alpha <- start$alpha
  N <- start$N
  lambda <- start$lambda
  mu.alpha <- start$mu.alpha  
  tau <- start$tau

  alpha.tmp <- sapply(1:qW,function(x) rep(alpha[,x],each=qY))
  p <- matrix(expit(apply(W*alpha.tmp,1,sum)),m,qY,byrow=TRUE)
  

  ###
  ###  Begin MCMC loop
  ###
  
  for(k in 1:n.mcmc){
    if(k%%1000==0) cat(k,"");flush.console()
        
    ###
    ###  Sample lambda 
    ###
    
    lambda <- rgamma(m,shape=N+priors$r,rate=1+priors$q)
    
    
    ###
    ###  Sample alpha[t] 
    ###
    
    alpha.star <- matrix(rnorm(qW*m,alpha,tune$alpha),m,qW)
    alpha.star.tmp <- sapply(1:qW,function(x) rep(alpha.star[,x],each=qY))
    p.star <- matrix(expit(apply(W*alpha.star.tmp,1,sum)),m,qY,byrow=TRUE)
    for(i in 1:m){
      mh.star.alpha <- sum(dbinom(Y[i,],N[i],p.star[i,],log=TRUE),na.rm=TRUE)+
        sum(dnorm(alpha.star[i,],mu.alpha,tau,log=TRUE))
      mh.0.alpha <- sum(dbinom(Y[i,],N[i],p[i,],log=TRUE),na.rm=TRUE)+
        sum(dnorm(alpha[i,],mu.alpha,tau,log=TRUE))
      if(exp(mh.star.alpha-mh.0.alpha)>runif(1)){
        alpha[i,] <- alpha.star[i,]
        keep$alpha <- keep$alpha+1
      }
    }
    alpha.tmp <- sapply(1:qW,function(x) rep(alpha[,x],each=qY))
    p <- matrix(expit(apply(W*alpha.tmp,1,sum)),m,qY,byrow=TRUE)

    ###
    ###  Sample tau
    ###

    tau.star <- rnorm(1,tau,tune$tau)
    if(tau.star>priors$a&tau.star<priors$b){
      mh.star.tau <- sum(dnorm(t(alpha),mu.alpha,tau.star,log=TRUE),na.rm=TRUE)
      mh.0.tau <- sum(dnorm(t(alpha),mu.alpha,tau,log=TRUE),na.rm=TRUE)
      if(exp(mh.star.tau-mh.0.tau)>runif(1)){
        tau <- tau.star
        keep$tau <- keep$tau+1
      } 
    }
    
    
    ###
    ###  Sample mu.alpha 
    ###
#     browser()
    
    A <- m*solve(tau^2*diag(qW))+solve(zeta^2*diag(qW))
    b <- rowSums(apply(alpha,1,function(x) x%*%solve(tau^2*diag(qW))))+solve(zeta^2*diag(qW))
#       +priors$mu.0%*%solve(zeta^2*diag(qW))
    mu.alpha <- rnorm(2,solve(A)%*%t(b),diag(solve(A)))
#     mu.alpha <- start$mu.alpha
    
#     alpha[1,]%*%solve(tau^2*diag(qW))  
#   
#     exp(0.5*(2*3))*exp(0.5*(2*4))
#     exp(0.5*(2*(3+4)))
#   m 
#   m*t(mu.alpha)%*%solve(tau^2*diag(2))%*%mu.alpha
#   t(mu.alpha)%*%(m*solve(tau^2*diag(2)))%*%mu.alpha
#     
    ###
    ###  Sample N 
    ###
    
    N.star <- N + sample(N.tune,m,replace=TRUE)    
    idx <- which(N.star>y.max)
      for(i in idx){
        mh.star.N <- sum(dbinom(Y[i,],N.star[i],p[i,],log=TRUE),na.rm=TRUE)+dpois(N.star[i],lambda[i],log=TRUE)
        mh.0.N <- sum(dbinom(Y[i,],N[i],p[i,],log=TRUE),na.rm=TRUE)+dpois(N[i],lambda[i],log=TRUE)  
        if(exp(mh.star.N-mh.0.N)>runif(1)){
          N[i] <- N.star[i]
          keep$N <- keep$N+1
        }
      }  
    
    
    ###
    ###  Save Samples 
    ###
    
    alpha.save[k,,] <- t(alpha)
    mu.alpha.save[k,] <- mu.alpha
    N.save[k,] <- N
    lambda.save[k,] <- lambda
    tau.save[k] <- tau
  }
  
  ###
  ###  Write Output 
  ###
  
  keep$tau <- keep$tau/n.mcmc
  keep$N <- keep$N/(n.mcmc*m)
  keep$alpha <- keep$alpha/(n.mcmc*m)
  cat(paste("\nN acceptance rate:",round(keep$N,2)))  
  cat(paste("\nalpha acceptance rate:",round(keep$alpha,2)))  
  cat(paste("\ntau acceptance rate:",round(keep$tau,2)))  
  list(mu.alpha=mu.alpha.save,alpha=alpha.save,N=N.save,lambda=lambda.save,tau=tau.save,keep=keep,n.mcmc=n.mcmc)
}
