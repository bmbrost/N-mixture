N.mixture.trend.MCMC <- function(Y,W,priors=list(r=15,q=0.1,tau=2,sigma=2),tune=list(N=5,alpha=0.01,theta=0.01),n.mcmc=1000){

  ###
  ### Brian M. Brost (17APR2015)
  ###
  ### N-mixture model for multiple sites with detection covariates (Royle 2004)
  ###
  ### Model statement:
  ### Y[i,t]~Binom(N[t],p[i,t])
  ### N[1]~Pois(lambda)
  ### N[t]=theta*N[t-1]
  ### lambda~Gamma(r,q)
  ### theta~N(mu,sigma^2)
  ### logit(p[i,j])~W[j,,i]%*%alpha
  ###
  ### Function arguments: 
  ### Y=m*J matrix, where m is the number of sites and J is the maximum number of
  ###   observations across all sites
  ### W=J*length(alpha)*m array; each 'slice' of W is a design matrix for model on detection
  ### priors=parameters of prior distributions for lambda and alpha
  ###
    
  ###
  ###  Setup Variables 
  ###

  #   browser()
  T <- nrow(Y) # Number of sites
  J <- ncol(Y) # Number of replicate counts
  y <- rowSums(Y) # Total of observed counts by site
  y.min <- apply(Y,1,min)
  keep <- list(N=0,alpha=0)
  qW <- ncol(W)
  W.mat <- apply(W,2,I) # Convert W from 3-D array to 2-D matrix
  
  N.save <- matrix(0,m,n.mcmc)
  lambda.save <- matrix(0,m,n.mcmc)
  alpha.save <- matrix(0,n.mcmc,qW)
  theta.save <- numeric(n.mcmc)
  
  ###
  ###  Priors and starting values 
  ###
  
  N.tune <- seq(-1*tune$N,tune$N,1)

  alpha <- rnorm(qW,0,priors$tau)
  alpha <- c(-5,0.75)
  N <- apply(Y,1,max)+1
  lambda <- N
  
  ###
  ###  Begin MCMC loop
  ###
  
  for(k in 1:n.mcmc){
    if(k%%1000==0) cat(k,"");flush.console()
        
    ###
    ###  Sample alpha 
    ###
# browser()
  
    alpha.star <- rnorm(qW,alpha,tune$alpha*I(qW))
    p.star <- expit(W.mat%*%alpha.star)
    N.tmp <- rep(N,each=J)  
    mh.star.alpha <- sum(dbinom(c(t(Y)),N.tmp,p.star,log=TRUE))+sum(dnorm(alpha.star,0,priors$tau,log=TRUE))
    mh.0.alpha <- sum(dbinom(c(t(Y)),N.tmp,p,log=TRUE))+sum(dnorm(alpha,0,priors$tau,log=TRUE))
    if(exp(mh.star.alpha-mh.0.alpha)>runif(1)){
      alpha <- alpha.star
      p <- expit(W.mat%*%alpha)
      keep$alpha <- keep$alpha+1
    }

    ###
    ###  Sample N for t=1 
    ###    
    
    N.star <- N[1] + sample(N.tune,1,replace=TRUE)
    if(N.star>=0 & N.star>y.min[1]){
      mh.star.N <- sum(dbinom(Y[1,],N.star,p,log=TRUE))+dpois(N.star[1],lambda[1],log=TRUE)
      mh.0.N <- sum(dbinom(Y[1,],N,p,log=TRUE))+dpois(N[1],lambda[1],log=TRUE)  
      if(exp(mh.star.N-mh.0.N)>runif(1){
        N[1] <- N.star
        keep$N <- keep$N+1
      }    
    }      

    ###
    ###  Sample N for t=2:T 
    ###    

    N.star <- N + sample(N.tune,m,replace=TRUE)

    idx <- which(N.star>=0 & N.star>y.min)
    N.star.tmp <- rep(N.star[idx],each=J)
    N.tmp <- rep(N[idx],each=J)
    n.tmp <- length(idx)
    #     browser()
    if(n.tmp>0){
      mh.star.N <- sum(dbinom(c(t(Y[idx,])),N.star.tmp,p,log=TRUE))+dpois(N.star[idx],lambda[idx],log=TRUE)
      mh.0.N <- sum(dbinom(c(t(Y[idx,])),N.tmp,p,log=TRUE))+dpois(N[idx],lambda[idx],log=TRUE)  
      idx <- idx[exp(mh.star.N-mh.0.N)>runif(n.tmp)]
      N[idx] <- N.star[idx]
      keep$N <- keep$N+length(idx)
    }      


    ####
    ####  Sample lambda 
    ####
    
    lambda <- rgamma(m,shape=N+priors$r,rate=1+priors$q)
    
    ####
    ####  Save Samples 
    ####
    
    alpha.save[k,] <- alpha
    N.save[,k] <- N
    lambda.save[,k] <- lambda
  }
  
  ####
  ####  Write Output 
  ####
  
  keep$N <- keep$N/(n.mcmc*m)
  keep$alpha <- keep$alpha/n.mcmc
  cat(paste("\nN acceptance rate:",keep$N))  
  cat(paste("\nalpha acceptance rate:",keep$alpha))  
  list(alpha=alpha.save,N=N.save,lambda=lambda.save,keep=keep,n.mcmc=n.mcmc)
}
