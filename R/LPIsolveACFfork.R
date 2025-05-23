## file:LPIsolveACFfork.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

## 
## Read data for one integration period, deconvolve the lag profiles
## in a fork cluster, and write the returned
## ACF to file. Repeat until end of data for every LPIparam$Ncluster'th integration period, starting from intPerFirst
## 
## Arguments: 
##  intPerFirst  Integration period number to start from, counted from
##             LPIparam[["firstTime"]] in steps of
##             LPIparam[["timeRes.s"]]
##
## Returns:
##  
## 

LPIsolveACFfork <- function( intPerFirst , LPIparam )
{
    # Load packages that are needed for reading the data
    for( pn in LPIparam[["inputPackages"]] ){
        require( pn , character.only=TRUE )
    }

    ## Parameter list update
    ## the default noUpdate will return NULL if the update is done twice for the same data
    LPIparam <- eval( as.name( LPIparam[["paramUpdateFunction"]] ))( LPIparam , intPeriod )


    if( !is.null(LPIparam)){

        ## Initialize a list for unsolved integration periods
        intPer.missing <- seq( intPerFirst , LPIparam[["lastIntPeriod"]] , by=LPIparam[['Ncluster']] )

        ## Run analysis loop until end of data
        endOfData <- FALSE

        repeat{
            
##            tt <- system.time({
            
            ## Update the last available data samples
            LPIparam[["dataEndTimes"]] <- eval( as.name( LPIparam[["dataEndTimeFunction"]] ))( LPIparam )
        
            ## Latest integration period for which data is available
            LPIparam[["maxIntPeriod"]] <- floor( ( min(unlist(LPIparam[["dataEndTimes"]])) - LPIparam[["startTime"]] ) / LPIparam[["timeRes.s"]] )

            ##  Select integration period number for the next analysis run
            ## Latest periods will be analysed first in order to simplify real-time analysis
            waitSum <- 0

            while( is.null( intPeriod <- nextIntegrationPeriods( LPIparam , 1 , intPer.missing ))){

                ## Break the loop after waiting
                ## long enough for new data
                if( waitSum > LPIparam[["maxWait.s"]] ){
                    endOfData <- TRUE
                    break
                }
        
                ## Wait 10 seconds
                Sys.sleep(10)
        
                ## Increment the wait time counter
                waitSum <- waitSum + 10
        
                ## Update the last available data samples
                LPIparam[["dataEndTimes"]] <- eval( as.name( LPIparam[["dataEndTimeFunction"]] ))( LPIparam )
 
                ## Latest integration period for which data is available
                LPIparam[["maxIntPeriod"]] <- floor( ( min(unlist(LPIparam[["dataEndTimes"]])) - LPIparam[["startTime"]] ) / LPIparam[["timeRes.s"]] )
            
            }
        
            if( endOfData ) break
      
            ## RprofFile <- paste('Rprof_',intPeriod,'.out',sep='')
            ## Rprof(filename=RprofFile,memory.profiling=TRUE,gc.profiling=TRUE,line.profiling=TRUE)
            

            ## Read raw data, name of the data input function
            ## should be stored in a character string
            LPIdatalist.raw   <- eval(as.name(LPIparam[["dataInputFunction"]]))( LPIparam , intPeriod )
            
            ## If data reading was successfull
            if(LPIdatalist.raw[["success"]]){
                
                ## require that there are at least some TX and RX samples
                if( (sum(LPIdatalist.raw[["RX1"]][["idata"]]) > 0) &
                    (sum(LPIdatalist.raw[["RX2"]][["idata"]]) > 0) &
                    (sum(LPIdatalist.raw[["TX1"]][["idata"]]) > 0) &
                    (sum(LPIdatalist.raw[["TX2"]][["idata"]]) > 0)){

                    analysisTime <- system.time({
                        
                        ## Frequency mixing, filtering, etc.
                        ## RprofFile <- paste('Rprof_',intPeriod,'.out',sep='')
                        ## Rprof(filename=RprofFile,memory.profiling=TRUE,gc.profiling=TRUE,line.profiling=TRUE)
                        
                        LPIdatalist.final <<- prepareLPIdata( LPIparam , LPIdatalist.raw )
                        
                        ## add some missing vectors and convert into an environment in the global workspace
                        if(LPIparam[["Rcomplex"]]){
                            initLPIenv(substitute(LPIdatalist.final))
                        }else{
                            initLPIenvR(substitute(LPIdatalist.final))
                        }                    
                        
                        ## Number of lags, each full lag
                        ## will get its own call of LPIsolve
                        nlags <- LPIdatalist.final[["nLags"]]
                        x <- seq( nlags )
                        
                        ## Number of range gates
                        ngates <- LPIdatalist.final[['nGates']]
                        maxgates <- max(ngates)
                        
                        ## Are we going to calculate a full covariance matrix?
                        fullcovar <- LPIdatalist.final[['fullCovar']]
                        
                        ## Range-gate centre points
                        r <- LPIdatalist.final[['rangeLimits']]
                        rgates <- ( r[1:maxgates] + r[2:(maxgates+1)] -1 ) / 2
                        
                        ## Lag-gate centre points
                        l <- LPIdatalist.final[["lagLimits"]]
                        lgates <- ( l[1:nlags] + l[2:(nlags+1)] -1 ) / 2
                        
                        
                        ## run the actual analysis in parallel using all available cores
                        if( is.null(LPIparam$nCores)){
                        ncl <- parallelly::availableCores()
                        }else{
                            ncl <- LPIparam$nCores
                        }
                        ##ACFlist <- parallel::mclapply( x , FUN=LPI:::LPIsolve , LPIenv.name=substitute(LPIdatalist.final) , mc.cores=ncl )
                                        #                    analysisTime <- system.time({
                        ACFlist <- parallel::mclapply( x , FUN=LPI:::LPIsolve , LPIenv.name=substitute(LPIdatalist.final) , intPeriod=intPeriod, mc.cores=ncl )
                                        #                    })
                        ##                    analysisTime <- NA
                        
                        ## sum of the flop counters
                        FLOP <- 0
                        ## time used for adding the theory lines to the solver
                                        #addTime <- 0
                        ## Collect the lag numbers from ACF list
                        lagnums <- x
                        for(k in 1:nlags ){
                            lagnums[k] <- ACFlist[[k]][['lagnum']]
                            FLOP <- FLOP + ACFlist[[k]][["FLOPS"]]
                                        #   addTime <- addTime + ACFlist[[k]][["addtime"]]
                        }
                        
                        ## Find correct order for the lag profiles
                        lagorder <- x[order(lagnums)]
                        
                        ## Order the ACF list
                        ACFlist <- ACFlist[lagorder]
                        
                        ## Make ACF and variance matrices
                        ACFmat <- matrix(NA,ncol=nlags,nrow=(maxgates+1))
                        
                        lagFLOP <- rep(NA,nlags)
                                        #lagAddTime <- list()
                        
                        ## Collect the lag profiles to the ACF matrix
                        for( k in 1:nlags){
                            if(ngates[k]>0){
                                ## Copy the solved lag profile
                                ACFmat[1:ngates[k],k] <- ACFlist[[k]][['lagprof']][1:ngates[k]]
                                ## Copy the background ACF estimate
                                ACFmat[maxgates+1,k]  <- ACFlist[[k]][['lagprof']][ngates[k]+1]
                                lagFLOP[k] <- ACFlist[[k]][["FLOPS"]]
                                        #lagAddTime[[k]] <- ACFlist[[k]][["addtime"]]
                            }
                        }
                        
                        ## If full covariance matrices were solved
                        if(fullcovar){
                            ## allocate matrix for variances and a cube for the covariance matrices
                            VARmat   <- matrix(NA,ncol=nlags,nrow=(maxgates+1))
                            COVARmat <- array(NA,dim=c((maxgates+1),(maxgates+1),nlags))
                            for( k in 1:nlags){
                                if(ngates[k]>0){
                                    ## Copy variances
                                    VARmat[1:ngates[k],k]                  <- Re(diag(ACFlist[[k]][['covariance']]))[1:ngates[k]]
                                    VARmat[maxgates+1,k]                   <- Re(diag(ACFlist[[k]][['covariance']]))[ngates[k]+1]
                                    ## Copy covariance matrices
                                    COVARmat[1:ngates[k],1:ngates[k],k]    <- ACFlist[[k]][['covariance']][1:ngates[k],1:ngates[k]]
                                    COVARmat[(maxgates+1),1:ngates[k],k]   <- ACFlist[[k]][['covariance']][(ngates[k]+1),1:ngates[k]]
                                    COVARmat[1:ngates[k],(maxgates+1),k]   <- ACFlist[[k]][['covariance']][1:ngates[k],(ngates[k]+1)]
                                    COVARmat[(maxgates+1),(maxgates+1),k]  <- ACFlist[[k]][['covariance']][(ngates[k]+1),(ngates[k]+1)]
                                }
                            }
                            ## If only variances were solved
                        }else{
                            ## Allocate a matrix for the variances,
                            ## set COVARmat to NULL
                            VARmat   <- matrix(NA,ncol=nlags,nrow=(maxgates+1))
                            COVARmat <- NULL
                            for( k in 1:nlags){
                                if( ngates[k] > 0 ){
                                    ## Copy the variances
                                    VARmat[1:ngates[k],k]  <- Re(ACFlist[[k]][['covariance']])[1:ngates[k]]
                                    VARmat[(maxgates+1),k] <- Re(ACFlist[[k]][['covariance']])[ngates[k]+1]
                                }
                            }
                        }
                        
                    })
                    
                    ## Collect the results in a list
                    ACFreturn <- list()
                    ACFreturn[["ACF"]]        <- ACFmat
                    ACFreturn[["var"]]        <- VARmat
                    ACFreturn[["covariance"]] <- COVARmat
                    ACFreturn[["lag"]]        <- lgates
                    ACFreturn[["range"]]      <- rgates
                    ACFreturn[["nGates"]]     <- ngates
                    ACFreturn[["FLOP"]] <- FLOP
                    ACFreturn[["analysisTime"]] <- analysisTime
                    #ACFreturn[["addTime"]] <- addTime
                    ACFreturn[["lagFLOP"]] <- lagFLOP
                    #ACFreturn[["lagAddTime"]] <- lagAddTime
                    
                    ## Store the results
                    eval( as.name( LPIparam[["resultSaveFunction"]]) )( LPIparam , intPeriod , ACFreturn )

##                    Rprof(NULL)
                    
            
                }
            }
        
            ## Remove the solved period from the list of missing ones
            intPer.missing <- setdiff( intPer.missing , intPeriod )

##        })
##    tfile <- file.path(LPIparam[["resultDir"]],sprintf("LPItimes-%05i.txt",intPerFirst))
##    cat(sprintf("%10s",names(tt)),file=tfile,append=T);cat('\n',file=tfile,append=T);cat(sprintf("%10.3f",tt),file=tfile,append=T);cat('\n',file=tfile,append=T)

            ## Stop if all integration periods are solved
            if( length(intPer.missing)==0) break

        } # repeat
        
    }

        
    
    ## Return the integration period
    ## number to the main process
    ##return(intPeriod)
    
  }
