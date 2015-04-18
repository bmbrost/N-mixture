base.N.mixture.MCMC <- function(Y,priors=list(a=1,b=1,r=15,q=0.1),tune=list(N=5),n.mcmc=1000){
  
  ###
  ### Brian M. Brost (17APR2015)
  ###
  ### Basic N-mixture model (Royle 2004)
  ###
  ### Arguments: 
  ### Y=m*J matrix, where m is the number of sites and J is the maximum number of
  ###   observations across all sites
  ### priors=parameters of prior distributionsfor p and lambda
  ###
  ### Model statement:
  ### Y[i,j]~Binom(N[j],p)
  ### N[j]~Pois(lambda[j])
  ### lambda[j]~Gamma(r,q)
  ### p~beta(a,b)
  ###
  
  ###
  ###  Setup Variables 
  ###

  #   browser()
  m <- nrow(Y) #Number of sites
  y <- rowSums(Y) #Total of observed counts by site
  y.min <- apply(Y,1,min)
  keep <- 0
  
  N.save <- matrix(0,m,n.mcmc)
  lambda.save <- matrix(0,m,n.mcmc)
  p.save <- matrix(0,m,n.mcmc)
  
  ###
  ###  Priors and starting values 
  ###
  
  N.tune <- seq(-1*tune$N,tune$N,1)

  N <- apply(Y,1,max)+1
  lambda <- N
  
  ###
  ###  Begin MCMC loop
  ###
  
  for(k in 1:n.mcmc){
    if(k%%1000==0) cat(k,"");flush.console()
        
    ####
    ####  Sample p 
    ####
# browser()
    p <- rbeta(m,y+priors$a,sum(N-Y)+priors$b)

    
    ####
    ####  Sample N 
    ####
    
#     N <- min(Y+rpois(m,(1-p)*lambda))
#     N <- 167
#  browser()

  N.star <- N + sample(N.tune,m,replace=TRUE)
#   N.star <- N + sample(seq(-5,5,1),1,replace=TRUE)  
  if(N.star>=0& N.star>y.min){
#       One row in Y only
      mh.star.N <- sum(dbinom(Y,N.star,p,log=TRUE))+dpois(N.star,lambda,log=TRUE)
      mh.0.N <- sum(dbinom(Y,N,p,log=TRUE))+dpois(N,lambda,log=TRUE)
      if(exp(mh.star.N-mh.0.N)>runif(1)){
        N <- N.star
        keep <- keep+1
      }
    }
  

    ####
    ####  Sample lambda 
    ####
    
    lambda <- rgamma(m,shape=N+priors$r,rate=1+priors$q)
    
    ####
    ####  Save Samples 
    ####
    
    p.save[k] <- p
    N.save[,k] <- N
    lambda.save[,k] <- lambda
    
  }
  
  ####
  ####  Write Output 
  ####
  
  list(p=p.save,N=N.save,lambda=lambda.save,keep=keep/n.mcmc,n.mcmc=n.mcmc)
}
