N.mixture.MCMC <- function(Y,W,priors,tune,start,n.mcmc=1000){

  ###
  ### Brian M. Brost (25APR2015)
  ###
  ### N-mixture model for multiple sites with detection modeled by covariates
  ###
  ### Model statement: (i indexes site, j indexes observation)
  ### Y[i,j]~Binom(N[i],p[i,j])
  ### N[i]~Pois(lambda[i])
  ### lambda[i]~Gamma(r,q)
  ### logit(p[i,j])~W[j,,i]%*%alpha
  ### alpha ~ N(mu,tau^2*I)
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
  keep <- list(N=0,alpha=0)

  N.save <- matrix(0,n.mcmc,m)
  lambda.save <- matrix(0,n.mcmc,m)
  alpha.save <- matrix(0,n.mcmc,qW)
  
  ###
  ###  Tunings and starting values 
  ###
  
  N.tune <- seq(-1*tune$N,tune$N,1)

  alpha <- start$alpha
  N <- start$N
  lambda <- start$lambda
  p <- matrix(expit(W%*%alpha),m,qY,byrow=TRUE)
  
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
    ###  Sample alpha 
    ###

    alpha.star <- rnorm(qW,alpha,tune$alpha*I(qW))
    p.star <- matrix(expit(W%*%alpha.star),m,qY,byrow=TRUE)
    mh.star.alpha <- sum(dbinom(Y,N,p.star,log=TRUE),na.rm=TRUE)+sum(dnorm(alpha.star,0,priors$tau,log=TRUE))
    mh.0.alpha <- sum(dbinom(Y,N,p,log=TRUE),na.rm=TRUE)+sum(dnorm(alpha,0,priors$tau,log=TRUE))
    if(exp(mh.star.alpha-mh.0.alpha)>runif(1)){
      alpha <- alpha.star
      p <- matrix(expit(W%*%alpha),m,qY,byrow=TRUE)
      keep$alpha <- keep$alpha+1
    }

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
    
    alpha.save[k,] <- alpha
    N.save[k,] <- N
    lambda.save[k,] <- lambda
  }
  
  ###
  ###  Write Output 
  ###
  
  keep$N <- keep$N/(n.mcmc*m)
  keep$alpha <- keep$alpha/n.mcmc
  cat(paste("\nN acceptance rate:",round(keep$N,2)))  
  cat(paste("\nalpha acceptance rate:",round(keep$alpha,2)))  
  list(alpha=alpha.save,N=N.save,lambda=lambda.save,keep=keep,n.mcmc=n.mcmc)
}
