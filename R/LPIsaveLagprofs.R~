## file:LPIrun.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Copy the raw data to computing slaves and solve lag
## profiles in there. Combine the solved profiles into
## a full ACF and return it.
##
## Arguments:
##  LPIenv.name Name of the LPI environment to use
##              for the analysis.
##
## Returns:
##  ACFlist     A list that contains the solved ACF, its
##              covariance matrices, lags, etc. 
## 
## The LPI environment, which may be a named list
## as well, must be stored on the global workspace. 
##

LPIrun <- function( LPIenv.name)
  {

    # Number of lags, each full lag
    # will get its own call of LPIsolve
    nlags <- eval(LPIenv.name)[["nLags"]]
    x <- seq( nlags )

    # Number of range gates
    ngates <- eval(LPIenv.name)[['nGates']]
    maxgates <- max(ngates)

    # Are we going to calculate a full covariance matrix?
    fullcovar <- eval(LPIenv.name)[['fullCovar']]

    # Range-gate centre points
    r <- eval(LPIenv.name)[['rangeLimits']]
    rgates <- ( r[1:maxgates] + r[2:(maxgates+1)] -1 ) / 2

    # Lag-gate centre points
    l <- eval(LPIenv.name)[["lagLimits"]]
    lgates <- ( l[1:nlags] + l[2:(nlags+1)] -1 ) / 2

    # If the computing slaves do not exist, set
    # slavecl=NA and run the analysis on this process
    if(!exists('slavecl')) slavecl <- NA
  
    # If we are running on a cluster
    if( eval(LPIenv.name)[["iscluster"]]  & ( !is.na( slavecl ) ) ){
    
      # Copy the data to all computing slaves
      clusterExport( slavecl , paste( LPIenv.name ) )

      # Allocate necessary vectors on each slave
      clusterCall( slavecl , initLPIenv , LPIenv.name )

      # Run the analysis processes on the slaves
      ACFlist <-  clusterApplyLB( slavecl , x , fun=LPI:::LPIsolve , LPIenv.name=LPIenv.name )

    
    # If not running on cluster, solve the lag profiles
    # sequentially. Mimic the output list of cluster
    # calls in order to simplify further processing.
    }else{ 

      # Create a list for the lag profiles
      ACFlist <- vector(mode='list',length=nlags)

      # Allocate vectors etc. 
      initLPIenv( LPIenv.name )

      # Run the actual analysis sequentially
      for( k in 1:nlags ){
        ACFlist[[k]] <- LPI:::LPIsolve( lag=x[k] , LPIenv.name=LPIenv.name )
      }

    }
    
    # Collect the lag numbers from ACF list
    lagnums <- x
    for(k in 1:nlags ){
      lagnums[k] <- ACFlist[[k]][['lagnum']]
    }

    # Find correct order for the lag profiles
    lagorder <- x[order(lagnums)]

    # Order the ACF list
    ACFlist <- ACFlist[lagorder]

    # Make ACF and variance matrices
    ACFmat <- matrix(NA,ncol=nlags,nrow=(maxgates+1))

    # Collect the lag profiles to the ACF matrix
    for( k in 1:nlags){
      if(ngates[k]>0){
        # Copy the solved lag profile
        ACFmat[1:ngates[k],k] <- ACFlist[[k]][['lagprof']][1:ngates[k]]
        # Copy the background ACF estimate
        ACFmat[maxgates+1,k]  <- ACFlist[[k]][['lagprof']][ngates[k]+1]
      }
    }

    # If full covariance matrices were solved
    if(fullcovar){
      # allocate matrix for variances and a cube for the covariance matrices
      VARmat   <- matrix(NA,ncol=nlags,nrow=(maxgates+1))
      COVARmat <- array(NA,dim=c((maxgates+1),(maxgates+1),nlags))
      for( k in 1:nlags){
        if(ngates[k]>0){
          # Copy variances
          VARmat[1:ngates[k],k]                  <- Re(diag(ACFlist[[k]][['covariance']]))[1:ngates[k]]
          VARmat[maxgates+1,k]                   <- Re(diag(ACFlist[[k]][['covariance']]))[ngates[k]+1]
          # Copy covariance matrices
          COVARmat[1:ngates[k],1:ngates[k],k]    <- ACFlist[[k]][['covariance']][1:ngates[k],1:ngates[k]]
          COVARmat[(maxgates+1),1:ngates[k],k]   <- ACFlist[[k]][['covariance']][(ngates[k]+1),1:ngates[k]]
          COVARmat[1:ngates[k],(maxgates+1),k]   <- ACFlist[[k]][['covariance']][1:ngates[k],(ngates[k]+1)]
          COVARmat[(maxgates+1),(maxgates+1),k]  <- ACFlist[[k]][['covariance']][(ngates[k]+1),(ngates[k]+1)]
        }
      }
    # If only variances were solved
    }else{
      # Allocate a matrix for the variances,
      # set COVARmat to NULL
      VARmat   <- matrix(NA,ncol=nlags,nrow=(maxgates+1))
      COVARmat <- NULL
      for( k in 1:nlags){
        if( ngates[k] > 0 ){
          # Copy the variances
          VARmat[1:ngates[k],k]  <- Re(ACFlist[[k]][['covariance']])[1:ngates[k]]
          VARmat[(maxgates+1),k] <- Re(ACFlist[[k]][['covariance']])[ngates[k]+1]
        }
      }
    }


    # Collect the results in a list and return it.
    # A list is used because an environment
    # is much slover to transfer
    ACFreturn <- list()
    ACFreturn[["ACF"]]        <- ACFmat
    ACFreturn[["var"]]        <- VARmat
    ACFreturn[["covariance"]] <- COVARmat
    ACFreturn[["lag"]]        <- lgates
    ACFreturn[["range"]]      <- rgates
    ACFreturn[["nGates"]]     <- ngates
    
    return(ACFreturn)
    
  }
