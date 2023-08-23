## file:LPIsaveLagprofs.R
## (c) 2023- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Merge the solved lag profiles into matrices and store to files
##
## Arguments:
##  LPIparam   the LPI parameter list
##  lagprofs   a list of deconvolved lag profiles from LPIrunLagprof
##  iper       integration period numbers included in these data
##  ii         The index vector used in the cluster call to LPIrunLagprof.
##              (the actual lags and integration periods are decoded from this vector)
##  ngates     number of range gates at each lag
##  nlags      number of time lags
##
## Returns:
##  nothing, the results are written t file
## 
##

LPIsaveLagprofs <- function( LPIparam , lagprofs , iper , ii , ngates , nlags)
  {

    # Are we going to calculate a full covariance matrix?
    fullcovar <- LPIparam[['fullCovar']]

    # maximum number of range gates
      maxgates  <- max(ngates)
      
    # Range-gate centre points
    r <- LPIparam[['rangeLimits']]
    rgates <- ( r[1:maxgates] + r[2:(maxgates+1)] -1 ) / 2

    # Lag-gate centre points
    l <- LPIparam[["lagLimits"]]
    lgates <- ( l[1:nlags] + l[2:(nlags+1)] -1 ) / 2

      
    # Collect the lag numbers and index vector from ACF list
    lagnums <- ii
    iinds <- ii  
      for(k in 1:length(ii) ){
          if(!is.null(lagprofs[[k]])){
              lagnums[k] <- lagprofs[[k]][['lagnum']]
              iinds[k] <- lagprofs[[k]][['ii']]
          }else{
              lagnums[k] <- NA
              iinds[k] <- NA
          }
    }

    # the actual integration period indices
      iperinds <- ceiling( iinds / nlags )
      iperinds <- iperinds[!is.na(iperinds)]

      x <- seq(nlags)

      
    # form the ACF matrices one integration period at a time

      for(kk in unique(iperinds)){

          # lag profiles from this integration period
          curinds  <- iperinds==kk
      
          # Find correct order for the lag profiles
          lagorder <- x[order(lagnums[curinds])]

          # Order the ACF list
          ACFlist <- lagprofs[curinds][lagorder]

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
    
          
      # Store the results
      eval( as.name( LPIparam[["resultSaveFunction"]]) )( LPIparam , iper[kk] , ACFreturn )

      }
      
  }
