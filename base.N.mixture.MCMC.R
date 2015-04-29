base.N.mixture.MCMC <- function(Y,priors,tune,start,n.mcmc=1000){
  
  ###
  ### Brian M. Brost (17APR2015)
  ###
  ### Basic N-mixture model for multiple sites
  ###
  ### Model statement: (i indexes site, j indexes observation)
  ### Y[i,j]~Binom(N[i],p[i])
  ### N[i]~Pois(lambda[i])
  ### lambda[i]~Gamma(r,q)
  ### p[i]~beta(a,b)
  ###
  ### Function arguments: 
  ### Y=m*J matrix, where m is the number of sites and J is the maximum number of
  ###   observations across all sites
  ### priors=parameters of prior distributionsfor p and lambda
  ### tune=tuning parameter for N
  ### start=starting values for N, p, and lambda
    
  ###
  ###  Setup Variables 
  ###

#     browser()
  m <- nrow(Y) #Number of sites
  y <- rowSums(Y,na.rm=TRUE) #Total of observed counts by site
  y.max <- apply(Y,1,max,na.rm=TRUE)
  keep <- 0
  
  N.save <- matrix(0,n.mcmc,m)
  lambda.save <- matrix(0,n.mcmc,m)
  p.save <- matrix(0,n.mcmc,m)
  
  ###
  ###  Priors and starting values 
  ###
  
  N.tune <- seq(-1*tune$N,tune$N,1)
  
  N <- start$N
  lambda <- start$lambda
  p <- start$p
  
  ###
  ###  Begin MCMC loop
  ###
  
  for(k in 1:n.mcmc){
    if(k%%1000==0) cat(k,"");flush.console()
        
    ###
    ###  Sample p 
    ###
    
    p <- rbeta(m,y+priors$a,sapply(1:m,function(x) sum(N[x]-Y[x,],na.rm=TRUE))+priors$b)
    

    ###
    ###  Sample lambda 
    ###
    
    lambda <- rgamma(m,shape=N+priors$r,rate=1+priors$q)
    
    
    ###
    ###  Sample N 
    ###
    
    N.star <- N + sample(N.tune,m,replace=TRUE)
    idx <- which(N.star>y.max)
    for(i in idx){
      mh.star.N <- sum(dbinom(Y[i,],N.star[i],p[i],log=TRUE),na.rm=TRUE)+dpois(N.star[i],lambda[i],log=TRUE)
      mh.0.N <- sum(dbinom(Y[i,],N[i],p[i],log=TRUE),na.rm=TRUE)+dpois(N[i],lambda[i],log=TRUE)  
      if(exp(mh.star.N-mh.0.N)>runif(1)){
        N[i] <- N.star[i]
        keep <- keep+1
      }
    }  
    
    
    ###
    ###  Save Samples 
    ###
    
    p.save[k,] <- p
    N.save[k,] <- N
    lambda.save[k,] <- lambda
  }
  
  ###
  ###  Write Output 
  ###

  keep <- keep/(n.mcmc*m)
  cat(paste("\nN Acceptance rate:",round(keep,2)))  
  list(p=p.save,N=N.save,lambda=lambda.save,keep=keep,n.mcmc=n.mcmc)
}
