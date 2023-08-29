## file:LPIsolveACFfork.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

## 
## Read data for one integration period, deconvolve the lag profiles
## in a fork cluster, and write the returned
## ACF to file
## 
## Arguments: 
##  intPeriod  Integration period number, counted from
##             LPIparam[["firstTime"]] in steps of
##             LPIparam[["timeRes.s"]]
##
## Returns:
##  intPeriod  The integration period number.
## 

LPIsolveACFfork <- function( intPeriod , LPIparam )
{
    # Load packages that are needed for reading the data
    for( pn in LPIparam[["inputPackages"]] ){
        require( pn , character.only=TRUE )
    }

    # Parameter list update
    LPIparam <- eval( as.name( LPIparam[["paramUpdateFunction"]] ))( LPIparam , intPeriod )

    if( !is.null(LPIparam)){
        # Read raw data, name of the data input function
        # should be stored in a character string
        LPIdatalist.raw   <- eval(as.name(LPIparam[["dataInputFunction"]]))( LPIparam , intPeriod )

        # If data reading was successfull
        if(LPIdatalist.raw[["success"]]){

          # require that there are at least some TX and RX samples
            if( (sum(LPIdatalist.raw[["RX1"]][["idata"]]) > 0) &
                (sum(LPIdatalist.raw[["RX2"]][["idata"]]) > 0) &
                (sum(LPIdatalist.raw[["TX1"]][["idata"]]) > 0) &
                (sum(LPIdatalist.raw[["TX2"]][["idata"]]) > 0)){
                
                # Frequency mixing, filtering, etc.
                LPIdatalist.final <<- prepareLPIdata( LPIparam , LPIdatalist.raw )

                # add some missing vectors and convert into an environment in the global workspace
                initLPIenv(substitute(LPIdatalist.final))


                # Number of lags, each full lag
                # will get its own call of LPIsolve
                nlags <- LPIdatalist.final[["nLags"]]
                x <- seq( nlags )

                # Number of range gates
                ngates <- LPIdatalist.final[['nGates']]
                maxgates <- max(ngates)

                # Are we going to calculate a full covariance matrix?
                fullcovar <- LPIdatalist.final[['fullCovar']]

                # Range-gate centre points
                r <- LPIdatalist.final[['rangeLimits']]
                rgates <- ( r[1:maxgates] + r[2:(maxgates+1)] -1 ) / 2

                # Lag-gate centre points
                l <- LPIdatalist.final[["lagLimits"]]
                lgates <- ( l[1:nlags] + l[2:(nlags+1)] -1 ) / 2


                ## # fork the cluster, try 20 times if the attemps fail
                ## ncl <- parallelly::availableCores()
                ## cl <- NULL
                ## itest <- 0
                ## while(length(cl)!=ncl){
                ##     cl <- tryCatch(
                ##         parallel::makeCluster(ncl,'FORK'),
                ##         error=function(e){
                ##             Sys.sleep(.1)
                ##             return(NULL)
                ##         }
                ##     )
                ##     itest  <- itest + 1
                ##     if(itest > 20){
                ##         return(NA)
                ##     }
                ## }
                
                ## # run the actual analysis in the cluster
                ## ACFlist <- parallel::clusterApply( cl , x , fun=LPI:::LPIsolve , LPIenv.name=substitute(LPIdatalist.final) )

                ## parallel::stopCluster(cl)

                # run the actual analysis in parallel using all available cores
                ncl <- parallelly::availableCores()
                ACFlist <- parallel::mclapply( x , FUN=LPI:::LPIsolve , LPIenv.name=substitute(LPIdatalist.final) , mc.cores=ncl )

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
                

                # Collect the results in a list
                ACFreturn <- list()
                ACFreturn[["ACF"]]        <- ACFmat
                ACFreturn[["var"]]        <- VARmat
                ACFreturn[["covariance"]] <- COVARmat
                ACFreturn[["lag"]]        <- lgates
                ACFreturn[["range"]]      <- rgates
                ACFreturn[["nGates"]]     <- ngates
                
                
                # Store the results
                eval( as.name( LPIparam[["resultSaveFunction"]]) )( LPIparam , intPeriod , ACFreturn )
                
            }
        }
    }
    
    # Return the integration period
    # number to the main process
    return(intPeriod)
    
  }
