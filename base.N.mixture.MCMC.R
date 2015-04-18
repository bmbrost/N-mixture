base.N.mixture.MCMC <- function(Y,priors=list(p=c(1,1),lamda=c(1,0.2)),n.mcmc){
  
  ###
  ### Brian M. Brost (17APR2015)
  ###
  ### Basic N-mixture model (Royle 2004)
  ###
  ### Arguments: 
  ### Y=m*J matrix, where m is the number of sites and J is the maximum number of
  ###   observations across all sites
  ### priors=prior distribution parameters for p, the probability of detection
  ###
  ### Model statement:
  ### Y[i,j]~Binom(N[j],p)
  ### N[j]~Pois(lambda[j])
  ### lambda[j]~Gamma(r,q)
  ### p~beta(alpha,beta)
  ###
  
  ###
  ###  Setup Variables 
  ###

  #   browser()
  m <- nrow(Y) #Number of sites
  y <- rowSums(Y) #Total of observed counts by site
  
  N.save <- matrix(0,m,n.mcmc)
  # N.total.save <- rep(0,n.mcmc)
  lambda.save <- matrix(0,m,n.mcmc)
  p.save <- matrix(0,m,n.mcmc)
  
  ###
  ###  Priors and starting values 
  ###
  
  N <- y+1
  
  N <- 150
  lambda <- 150
  
  ###
  ###  Begin MCMC loop
  ###
  keep <- 0
  for(k in 1:n.mcmc){
    if(k%%100==0) cat(k,"");flush.console()
        
    ####
    ####  Sample p 
    ####
#     browser()

    p <- rbeta(m,y+priors$p[1],sum(N-Y)+priors$p[2])
        
    
    ####
    ####  Sample N 
    ####
    
#     N <- min(Y+rpois(m,(1-p)*lambda))
#     N <- 167
#  browser()
  N.star <- N + sample(seq(-5,5,1),m)    
    if(N.star>=0&sum(N.star<Y)==0){
#       One row in Y only
      mh.star.N <- sum(dbinom(Y,N.star,p,log=FALSE))+dpois(N.star,lambda,log=TRUE)
      mh.0.N <- sum(dbinom(Y,N,p,log=TRUE))+dpois(N,lambda,log=TRUE)
      if(exp(mh.star.N-mh.0.N)>runif(1)){
        N <- N.star
        keep <- keep+1
      }
    }
    N <- 167

    ####
    ####  Sample lambda 
    ####
    
    r <- priors$lambda[1]
    q <- priors$lambda[2]
    lambda <- rgamma(m,shape=N+r,rate=1+q)
#     lambda <- 150
    
    ####
    ####  Save Samples 
    ####
    
    p.save[k] <- p
    N.save[,k] <- N
    lambda.save[,k] <- lambda
    
  }
  cat("\n")
  
  ####
  ####  Write Output 
  ####
  
  list(p=p.save,N=N.save,lambda=lambda.save)
  
}
