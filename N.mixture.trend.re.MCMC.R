N.mixture.trend.MCMC <- function(Y,W,priors,tune,start,n.mcmc=1000){  
  
  ###
  ### Brian M. Brost (25APR2015)
  ###
  ### N-mixture model with temporal trend
  ###
  ### Model statement: (i indexes observation and t indexes time)
  ### Y[i,t]~Binom(N[t],p[i,t])
  ### N[1]~Pois(lambda)
  ### N[t]~Pois(exp(log(theta)+log(N[t-1])) #log link
  ### lambda~Gamma(r,q)
  ### theta~N(mu_theta,sigma^2)
  ### logit(p[i,t])~W[i,,t]%*%alpha
  ### alpha~N(mu_alpha,tau^2*I)
  ###
  ### Note: an identity link for N[t] could also be used. Code for this alternative, i.e., 
  ###   N[t]~Pois(theta*N[t-1]), is 'commented' out below in updates for N[t], N[T], and theta 
  ###
  ### Note: this algorithm allows estimation of the posterior predictive distribution (i.e., smoother
  ###   distribution) for time periods in which observations are missing, provided observations are 
  ###   availabe for the preceding and following time periods.
  ###
  ### Function arguments: 
  ### t=times corresponding to rows in Y, i.e., times in which observations are available
  ### Y=T*m matrix, where m is the number of observations and T is the number of time periods
  ### W=m*length(alpha)*T array; each 'slice' of W is a design matrix for detection model during a time period
  ### priors=parameters of prior distributions for lambda, alpha, and theta
  ### tune=tuning parameters for N, alpha, and theta
  ### start=starting values for N,lambda, alpha, and theta
  ###

  
  ###
  ###  Libraries and subroutines
  ###
  
  expit <- function(x) 1/(1+exp(-x))
  logit <- function(x) log(x/(1-p))
  
  ###
  ###  Setup Variables 
  ###

#     browser()
  qY <- ncol(Y)
  qW <- ifelse(is.null(ncol(W)),1,ncol(W))
  T <- nrow(Y) # Number of time periods
  m <- apply(Y,1,function(x) sum(!is.na(x))) # Number of observations per time period
  y.max <- rep(0,T)
  y.max[m>0] <- apply(Y[m>0,],1,max,na.rm=TRUE)
  
  keep <- list(N=0,alpha=0,theta=0)

  N.save <- matrix(0,n.mcmc,T)
  lambda.save <- numeric(n.mcmc)
  alpha.save <- matrix(0,n.mcmc,qW)
  theta.save <- numeric(n.mcmc)
  
  ###
  ###  Tuning and starting values 
  ###
  
  N.tune <- seq(-1*tune$N,tune$N,1)  
  alpha <- start$alpha
  N <- start$N
  lambda <- start$lambda
  theta <- start$theta
  p <- matrix(expit(W%*%alpha),T,qY,byrow=TRUE)

  ###
  ###  Begin MCMC loop
  ###
  
  for(k in 1:n.mcmc){
    if(k%%1000==0) cat(k,"");flush.console()
    
    ###
    ### Note: lines pertaining to the Poisson density in the updates for N[2:(t-1)], N[T], and theta
    ### that are 'commented' out use the identity link, i.e., N[t]~Pois(theta*N[t-1])
    ###
    
    ###
    ###  Sample alpha 
    ###
  
    alpha.star <- rnorm(qW,alpha,tune$alpha*I(qW))
    p.star <- matrix(expit(W%*%alpha.star),T,qY,byrow=TRUE)
    mh.star.alpha <- sum(dbinom(Y,N,p.star,log=TRUE),na.rm=TRUE)+sum(dnorm(alpha.star,0,priors$tau,log=TRUE))
    mh.0.alpha <- sum(dbinom(Y,N,p,log=TRUE),na.rm=TRUE)+sum(dnorm(alpha,0,priors$tau,log=TRUE))
    if(exp(mh.star.alpha-mh.0.alpha)>runif(1)){
      alpha <- alpha.star
      p <- p.star
      keep$alpha <- keep$alpha+1
    }

    
    ###
    ###  Sample N[1] (i.e., abundance for t=1)
    ###    

    N.star <- N[1] + sample(N.tune,1)
    if(N.star>y.max[1]){ # Update N[1] only if N.star>y.max[1]
      mh.star.N <- sum(dbinom(Y[1,],N.star,p[1,],log=TRUE),na.rm=TRUE)+dpois(N.star,lambda,log=TRUE)
      mh.0.N <- sum(dbinom(Y[1,],N[1],p[1,],log=TRUE),na.rm=TRUE)+dpois(N[1],lambda,log=TRUE)  
      if(exp(mh.star.N-mh.0.N)>runif(1)){
        N[1] <- N.star
        keep$N <- keep$N+1
      }    
    }      

    ###
    ###  Sample N[2:(T-1)] (i.e., abundance for t=2:(T-1))
    ###   
    
    for(t in 2:(T-1)){
      N.star <- N[t] + sample(N.tune,1)
      if(N.star>y.max[t]){ # Update N[t] only if N.star>y.max[t]
        mh.star.N <- sum(dbinom(Y[t,],N.star,p[t,],log=TRUE),na.rm=TRUE)+
          dpois(N.star,exp(log(theta)+log(N[t-1])),log=TRUE)+dpois(N[t+1],exp(log(theta)+log(N.star)),log=TRUE)
          # dpois(N.star,theta*N[t-1],log=TRUE)+dpois(N[t+1],theta*N.star,log=TRUE)
        mh.0.N <- sum(dbinom(Y[t,],N[t],p[t,],log=TRUE),na.rm=TRUE)+
          dpois(N[t],exp(log(theta)+log(N[t-1])),log=TRUE)+dpois(N[t+1],exp(log(theta)+log(N[t])),log=TRUE)
          # dpois(N[t],theta*N[t-1],log=TRUE)+dpois(N[t+1],theta*N[t],log=TRUE)
        if(exp(mh.star.N-mh.0.N)>runif(1)){
          N[t] <- N.star
          keep$N <- keep$N+1
        }    
      }
    }
    
    
    ###
    ###  Sample N[T] (i.e., abundance for t=T)
    ###    
    
    N.star <- N[T] + sample(N.tune,1)
    if(N.star>y.max[T]){ # Update N[T] only if N.star>y.max[T]
      mh.star.N <- sum(dbinom(Y[T,],N.star,p[T,],log=TRUE),na.rm=TRUE)+
        dpois(N.star,exp(log(theta)+log(N[T-1])),log=TRUE)
      mh.0.N <- sum(dbinom(Y[T,],N[T],p[T,],log=TRUE),na.rm=TRUE)+
        dpois(N[T],exp(log(theta)+log(N[T-1])),log=TRUE)  
      # mh.star.N <- sum(dbinom(Y[T,],N.star,p.tmp,log=TRUE))+dpois(N.star,theta*N[T-1],log=TRUE)
      # mh.0.N <- sum(dbinom(Y[T,],N[T],p.tmp,log=TRUE))+dpois(N[T],theta*N[T-1],log=TRUE)  
      if(exp(mh.star.N-mh.0.N)>runif(1)){
        N[T] <- N.star
        keep$N <- keep$N+1
      }    
    }      
    
    
    ###
    ###  Sample lambda  
    ###
    
    lambda <- rgamma(1,shape=N[1]+priors$r,rate=1+priors$q)
    
    
    ###
    ###  Sample theta
    ###

    theta.star <- rnorm(1,theta,tune$theta)
    mh.star.theta <- sum(dpois(N[-1],exp(log(theta.star)+log(N[-T])),log=TRUE))+
      dnorm(theta.star,1,priors$sigma,log=TRUE)
    mh.0.theta <- sum(dpois(N[-1],exp(log(theta)+log(N[-T])),log=TRUE))+
      dnorm(theta,1,priors$sigma,log=TRUE)
    # mh.star.theta <- sum(dpois(N[-1],theta.star*N[-T],log=TRUE))+dnorm(theta.star,1,priors$sigma,log=TRUE)
    # mh.0.theta <- sum(dpois(N[-1],theta*N[-T],log=TRUE))+dnorm(theta,1,priors$sigma,log=TRUE)
    if(exp(mh.star.theta-mh.0.theta)>runif(1)){
      theta <- theta.star
      keep$theta <- keep$theta+1
    }


    ###
    ###  Save Samples 
    ###
    
    alpha.save[k,] <- alpha
    N.save[k,] <- N
    lambda.save[k] <- lambda
    theta.save[k] <- theta
  }
  
  ###
  ###  Write Output 
  ###
  
  keep$N <- keep$N/(n.mcmc*T)
  keep$alpha <- keep$alpha/n.mcmc
  keep$theta <- keep$theta/n.mcmc
  cat(paste("\nN acceptance rate:",round(keep$N,2)))  
  cat(paste("\nalpha acceptance rate:",round(keep$alpha,2)))
  cat(paste("\ntheta acceptance rate:",round(keep$theta,2)))  
  list(alpha=alpha.save,N=N.save,lambda=lambda.save,theta=theta.save,keep=keep,n.mcmc=n.mcmc)
}
