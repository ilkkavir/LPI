## file:LPI.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.

##
## The main analysis loop of LPI
##
##
## 

LPI <- function(dataInputFunction,
                inputPackages=c(),
                startTime = 0, # 1st Jan 1970 00:00 UT
                stopTime = 4000000000, # 2nd Oct 2096 07:00 UT
                nup = LPIexpand.input( 1 ),
                filterLength = LPIexpand.input( 1 ),
                decodingFilter = "none",
                lagLimits = c(1,2),
                rangeLimits = c(1,2),
                maxRanges = Inf,
                maxClutterRange = 0,
                clutterFraction = 1,
                timeRes.s = 10,
                backgroundEstimate=TRUE,
                maxWait.s = -1,
                freqOffset = LPIexpand.input( 0 ),
                indexShifts = LPIexpand.input( list(c(0,0)) ),
                solver = "fishsr",
                nBuf = 10000,
                fullCovar = FALSE,
                rlips.options = list( type="c" , nbuf=1000 , workgroup.size=128),
                remoteRX = FALSE,
                normTX = FALSE,
                nCode = NA,
                ambInterp = FALSE,
                minNpower = 100,
                noiseSpikeThreshold = 5,
                resultDir = paste(format(Sys.time(),"%Y-%m-%d_%H:%M"),'LP',sep='_'),
                dataEndTimeFunction="currentTimes",
                resultSaveFunction = "LPIsaveACF",
                paramUpdateFunction="noUpdate",
                cl=NULL,
                nCores = NULL,
                ...
                ){
    
    # Collect all input in a list that is handy to pass forwards
    par1 <- formals()
    par1['...'] <- NULL
    par2 <- list(...)
    par1names <- names(par1)
    par1 <- lapply( names( par1 ) , FUN=function(x){ eval( as.name( x ) ) } )
    names(par1) <- par1names
    LPIparam <- c(par1,par2)

    # Expand parameters to LPI internal format and set storage modes as necessary
    LPIparam[["nup"]] <- LPIexpand.input( LPIparam[["nup"]] )
    storage.mode( LPIparam[["nup"]] ) <- "integer"
    LPIparam[["filterLength"]] <- LPIexpand.input( LPIparam[["filterLength"]] )
    storage.mode( LPIparam[["filterLength"]] ) <- "integer"
    storage.mode( LPIparam[["lagLimits"]] ) <- "integer"
    storage.mode( LPIparam[["rangeLimits"]] ) <- "integer"
    LPIparam[["maxClutterRange"]] <- LPIexpand.input( LPIparam[["maxClutterRange"]] )
    storage.mode( LPIparam[["maxClutterRange"]] ) <- "integer"
    LPIparam[["clutterFraction"]] <- LPIexpand.input( LPIparam[["clutterFraction"]] )
    LPIparam[["freqOffset"]] <- LPIexpand.input( LPIparam[["freqOffset"]] )
    if( ! is.list( LPIparam[["indexShifts"]] ) ){
      LPIparam[["indexShifts"]] <- list(LPIparam[["indexShifts"]])
    }
    LPIparam[["indexShifts"]] <- LPIexpand.input( LPIparam[["indexShifts"]] )
    for( dType in c("TX1","TX2","RX1","RX2")) storage.mode(LPIparam[["indexShifts"]][[dType]]) <- "integer"
    storage.mode( LPIparam[["nCode"]] ) <- "integer"
    storage.mode( LPIparam[["minNpower"]] ) <- "integer"


    # Print input arguments
    cat(sprintf("%20s %f (%s UT)\n","startTime:",startTime,format(as.POSIXlt(startTime,origin='1970-01-01',tz='ut'),"%Y-%m-%d %H:%M:%OS6")))
    cat(sprintf("%20s %f (%s UT)\n","stopTime:",stopTime,format(as.POSIXlt(stopTime,origin='1970-01-01',tz='ut'),"%Y-%m-%d %H:%M:%OS6")))
    cat(sprintf("%20s"," inputPackages:"))
    for(n in inputPackages){cat(n,", ")}
    cat('\n')
    cat(sprintf("%20s %s\n","dataInputFunction:",dataInputFunction))
    cat(sprintf("%20s %s\n","dataEndTimeFunction:",dataEndTimeFunction))
    ## cat(sprintf("%20s"," clusterNodes:"))
    ## if( is.list(clusterNodes )){
    ##   for(n in names(clusterNodes)){cat(sprintf("%s:",n));cat(clusterNodes[[n]],'  ')};cat('\n')
    ## }else{
    ##   cat( clusterNodes ); cat('\n')
    ## }
    cat(sprintf("%20s","nup:"));for(dType in c("RX1","RX2","TX1","TX2")){cat(' ',dType,':',LPIparam[["nup"]][[dType]],sep='')};cat('\n')
    cat(sprintf("%20s","filterLength:"));for(dType in c("RX1","RX2","TX1","TX2")){cat(' ',dType,':',LPIparam[["filterLength"]][[dType]],sep='')};cat('\n')
    cat(sprintf("%20s %s\n","decodingFilter:",decodingFilter[1]))
    cat(lagLimits,fill=70,labels=c(sprintf("%20s","lagLimits:"),rep('                    ',1000)))
    cat(rangeLimits,fill=70,labels=c(sprintf("%20s","rangeLimits:"),rep('                    ',1000)))
    cat(maxRanges,fill=70,labels=c(sprintf("%20s","maxRanges:"),rep('                    ',1000)))
    cat(sprintf("%20s %f\n","timeRes.s:",timeRes.s))
    cat(sprintf("%20s RX1:%i RX2:%i \n","maxClutterRange:",LPIparam$maxClutterRange["RX1"],LPIparam$maxClutterRange["RX2"]))
    cat(sprintf("%20s RX1:%i RX2:%i \n","clutterFraction:",LPIparam$clutterFraction["RX1"],LPIparam$clutterFraction["RX2"]))
    cat(sprintf("%20s %s\n","backgroundEstimate:",backgroundEstimate))
    cat(sprintf("%20s %f\n","maxWait.s:",maxWait.s))
    cat(sprintf("%20s RX1:%f RX2:%f TX1:%f TX2:%f\n","freqOffset:",LPIparam$freqOffset["RX1"],LPIparam$freqOffset["RX2"],LPIparam$freqOffset["TX1"],LPIparam$freqOffset["TX2"]))
    cat(sprintf("%20s","indexShifts:"));for(dType in c("RX1","RX2","TX1","TX2")){cat(' ',dType,':',sep='');cat(LPIparam$indexShifts[[dType]])};cat('\n')
    cat(sprintf("%20s %s\n","solver:",solver))
    cat(sprintf("%20s %i\n","nBuf:",nBuf))
    cat(sprintf("%20s %s\n","fullCovar:",fullCovar))
    cat(sprintf("%20s","rlips.options:"));for(n in names(rlips.options)){cat(' ',n,':',rlips.options[[n]],sep='')};cat('\n')
    cat(sprintf("%20s %s\n","remoteRX:",remoteRX))
    cat(sprintf("%20s %s\n","normTX:",normTX))
    cat(sprintf("%20s %i\n","nCode:",nCode))
    cat(sprintf("%20s %s\n","ambInterp:",ambInterp))
    cat(sprintf("%20s %s\n","resultDir:",resultDir))
    cat(sprintf("%20s %s\n","resultSaveFunction:",resultSaveFunction))
    cat(sprintf("%20s %s\n","paramUpdateFunction:",paramUpdateFunction))
#    cat(sprintf("%20s %s\n","useXDR:",useXDR))
    
    # Total number of integration periods requested
    LPIparam[["lastIntPeriod"]] <- round( ( stopTime - startTime ) / LPIparam[["timeRes.s"]] )

    # Create the result directory if a valid path was given
    if( is.character( resultDir ) ){
      if( nchar( resultDir ) > 0 ){
        dir.create( resultDir , recursive=TRUE , showWarnings=FALSE )
      }
    }


    ## check if Rcomplex or separate arrays of Re and Im data should be used
    if( any( LPIparam[["solver"]] == c("fishsr","decor") ) ){
        LPIparam[["Rcomplex"]] <- FALSE
    }else{
        LPIparam[["Rcomplex"]] <- TRUE
    }
        
      

    ## number of slave processes (one is automatically saved for the master process, we will use also that with help of future)
    if(is.null(cl)){
        Ncl <- 1
    }else{
        Ncl <- length(cl) + 1
    }
    

    ## set number of cores explicitly, since availableCores() will give an incorrect value when called within the future call
    LPIparam$nCores <- parallelly::availableCores()

    ## we cannot use the future trick with single core, set Ncl accordingly
    if(LPIparam$nCores==1){
        Ncl <- Ncl - 1
    }

    ## Ncl is assumed to be positive
    Ncl <- max(Ncl,1)

    LPIparam[["Ncluster"]] <- Ncl

    if(!is.na(LPIparam$resultDir)){
        save(LPIparam,file=file.path(resultDir,'LPIparam.Rdata'))
    }
    
    ## # find a reasonable number of parallel integration periods (Niper <= Ncl & Niper*Nlags >= Ncl)
    ## Nlags <- length(LPIparam[["lagLimits"]]) - 1

    ## let the cluster nodes do the work, except if this is not a cluster
    if (Ncl<=1){
##        print('Single core, LPIsolveACFfork')
        print( unlist( LPIsolveACFfork( 1 , LPIparam  ) ) )
##        print('Single core, LPIsolveACFfork, done')
    }else{
        ## start analysis in the current MPI node using future
        ## but this works only if there are two or more cores available.
        future::plan(multicore)
        if(LPIparam$nCores>1){
##            print('multicore, LPIsolveACFfork, future')
            futureOut <- future(LPIsolveACFfork( 1 , LPIparam  ) )
            ## start analysis in the other MPI nodes using clusterApply
##            print('multicore, LPIsolveACFfork, clusterApply')
            print(unlist(snow::clusterApply( cl , seq(2,Ncl) , fun=LPIsolveACFfork , LPIparam )))
            ## make sure that the analysis in the current MPI task is also finished
##            print('multicore, LPIsolveACFfork, clusterApply, done')
            value(futureOut)
##            print('multicore, LPIsolveACFfork, future, done')
        }else{
##            print('multicore, LPIsolveACFfork, clusterApply')
            print(unlist(snow::clusterApply( cl , seq(1,Ncl) , fun=LPIsolveACFfork , LPIparam )))
##            print('multicore, LPIsolveACFfork, clusterApply, done')
        }
    }
##    print('LPI done')

    ## this loop is now in LPIsolveACFfork
    
    ## # Initialize a list for unsolved integration periods
    ## intPer.missing <- seq( LPIparam[["lastIntPeriod"]] )

    ## # Run analysis loop until end of data
    ## endOfData <- FALSE
    ## repeat{

    ##     # Update the last available data samples
    ##     LPIparam[["dataEndTimes"]] <- eval( as.name( LPIparam[["dataEndTimeFunction"]] ))( LPIparam )
        
    ##     # Latest integration period for which data is available
    ##     LPIparam[["maxIntPeriod"]] <- floor( ( min(unlist(LPIparam[["dataEndTimes"]])) - LPIparam[["startTime"]] ) / LPIparam[["timeRes.s"]] )

    ##     # Select integration period numbers for the next analysis run
    ##     # Latest periods will be analysed first in order to simplify real-time analysis
    ##     waitSum <- 0
    ##     while( is.null( intPer.current <- nextIntegrationPeriods( LPIparam , Ncl , intPer.missing ))){

    ##         # Break the loop after waiting
    ##         # long enough for new data
    ##         if( waitSum > LPIparam[["maxWait.s"]] ){
    ##             endOfData <- TRUE
    ##             break
    ##         }
        
    ##         # Wait 10 seconds
    ##         Sys.sleep(10)
        
    ##         # Increment the wait time counter
    ##         waitSum <- waitSum + 10
        
    ##         # Update the last available data samples
    ##         LPIparam[["dataEndTimes"]] <- eval( as.name( LPIparam[["dataEndTimeFunction"]] ))( LPIparam )
 
    ##         # Latest integration period for which data is available
    ##         LPIparam[["maxIntPeriod"]] <- floor( ( min(unlist(LPIparam[["dataEndTimes"]])) - LPIparam[["startTime"]] ) / LPIparam[["timeRes.s"]] )
            
    ##     }
        
    ##     if( endOfData ) break
        
        
    ##     # run the integration periods in parallel in the MPI cluster
    ##     print(unlist(snow::clusterApply( cl , intPer.current , fun=LPIsolveACFfork , LPIparam )))
        
        
    ##     # should do the below comparison for the actual returned integration period numbers, from which we can easily exclude possibly failed ones and try them again.. 
        
    ##     # Remove the solved periods from the list of missing ones
    ##     intPer.missing <- setdiff( intPer.missing , intPer.current )

    ##     warnings()
    ##     # Stop if all integration periods are solved
    ##     if( length(intPer.missing)==0) break
    
    ## } # repeat



    
    # Shut down the cluster at end of analysis
#    if(!all(is.na(LPIparam[["clusterNodes"]]))) snow::stopCluster( cl )
#    snow::stopCluster( cl )

    # This function does not return anything,
    # results are written to files.
    invisible()

  }

