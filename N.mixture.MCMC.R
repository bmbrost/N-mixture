base.N.mixture.MCMC <- function(Y,priors=list(a=1,b=1,r=15,q=0.1),tune=list(N=5),n.mcmc=1000){
  
  ###
  ### Brian M. Brost (17APR2015)
  ###
  ### Base N-mixture model for multiple sites (Royle 2004)
  ###
  ### Model statement:
  ### Y[i,j]~Binom(N[j],p)
  ### N[j]~Pois(lambda[j])
  ### lambda[j]~Gamma(r,q)
  ### p~beta(a,b)
  ###
  ### Function arguments: 
  ### Y=m*J matrix, where m is the number of sites and J is the maximum number of
  ###   observations across all sites
  ### priors=parameters of prior distributionsfor p and lambda
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
        
    ###
    ###  Sample p 
    ###
# browser()
  Y  
  p <- rbeta(m,y+priors$a,sapply(1:m,function(x) sum(N[x]-Y[x,]))+priors$b)
# p <- rbeta(m,y+priors$a,sum(N-Y))+priors$b)


    ###
    ###  Sample N 
    ###
    
    # Note: the full conditional for N is conjugate, but is funky; therefore, use a M-H update
    # N <- min(Y+rpois(m,(1-p)*lambda))

    N.star <- N + sample(N.tune,m,replace=TRUE)
    idx <- which(N.star>=0 & N.star>y.min)
    n.tmp <- length(idx)
#     browser()
    if(length(idx)>0){
      mh.star.N <- rowSums(matrix(dbinom(Y[idx,],N.star[idx],p[idx],log=TRUE),n.tmp,J))+dpois(N.star[idx],lambda[idx],log=TRUE)
      mh.0.N <- rowSums(matrix(dbinom(Y[idx,],N[idx],p[idx],log=TRUE),n.tmp,J))+dpois(N[idx],lambda[idx],log=TRUE)  
      idx <- idx[exp(mh.star.N-mh.0.N)>runif(n.tmp)]
      N[idx] <- N.star[idx]
      keep <- keep+length(idx)
    }

#   if(length(idx)>0){
#   mh.star.N <- sum(dbinom(Y,N.star,p,log=TRUE))+dpois(N.star,lambda,log=TRUE)
#   mh.0.N <- sum(dbinom(Y,N,p,log=TRUE))+dpois(N,lambda,log=TRUE)
#   if(exp(mh.star.N-mh.0.N)>runif(1)){
#     N <- N.star
#     keep <- keep+1
#   }
# }

    ####
    ####  Sample lambda 
    ####
    
    lambda <- rgamma(m,shape=N+priors$r,rate=1+priors$q)
    
    ####
    ####  Save Samples 
    ####
    
    p.save[,k] <- p
    N.save[,k] <- N
    lambda.save[,k] <- lambda
  }
  
  ####
  ####  Write Output 
  ####

  cat("\n")
  cat(paste("Acceptance rate:",keep/(n.mcmc*m)))  
  list(p=p.save,N=N.save,lambda=lambda.save,keep=keep/(n.mcmc*m),n.mcmc=n.mcmc)
}
