## file:LPIsaveACF.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Save resolved ACF to file
##
## Arguments:
##  LPIparam  A LPI parameter list
##  intPeriod Integration period number
##  ACF       An ACF list returned by LPIsolve
##
## Returns:
##  resFile   Result file name
##

LPIsaveACF <- function( LPIparam , intPeriod , ACF )
  {
    # Number of range gates
    ngates <- length(ACF[["range"]])

    # Number of lags
    nlags  <- length(ACF[["lag"]])

    # Seconds since 1970
    ACF[["time.s"]] <-  LPIparam[["startTime"]] + intPeriod*LPIparam[["timeRes.s"]]

    # The same time as a string, useful for debugging
    # time conversions and for plotting
    ACF[["timeString"]] <-
      format( as.POSIXct( ACF[["time.s"]] , origin='1970-01-01' , tz='UTC' ) , "%Y-%m-%d %H:%M:%OS3 UT")

    # Result file name
    resFile <- gsub(' ','0',file.path( LPIparam[["resultDir"]] , paste( sprintf( '%13.0f' , trunc( ACF[["time.s"]]  * 1000 ) ) , "LP.Rdata" , sep='') ))

    # Range
    names(ACF[["range"]]) <- paste('gate',seq(ngates),sep='')

    # Lag
    names(ACF[["lag"]]) <- paste('lag',seq(nlags),sep='')

    # Background ACF
    ACF[["backgroundACF"]] <- ACF[["ACF"]][(ngates+1),]
    ACF[["backgroundvar"]] <- ACF[["var"]][(ngates+1),]
    names(ACF[["backgroundACF"]]) <- paste('lag',seq(nlags),sep='')
    names(ACF[["backgroundvar"]]) <- paste('lag',seq(nlags),sep='')

    # ACF and variance without the background samples
    ACF[["ACF"]] <- matrix(ACF[["ACF"]][1:ngates,],ncol=nlags)
    ACF[["var"]] <- matrix(ACF[["var"]][1:ngates,],ncol=nlags)
    dimnames(ACF[["ACF"]]) <- list(paste('gate',seq(ngates),sep=''),paste('lag',seq(nlags),sep=''))
    dimnames(ACF[["var"]]) <- list(paste('gate',seq(ngates),sep=''),paste('lag',seq(nlags),sep=''))

    # Dimnames for the optional full covariance matrix
    if(LPIparam[["fullCovar"]]) dimnames(ACF[["covariance"]]) <- list( c(paste('gate',seq(ngates),sep=''),'background') , c(paste('gate',seq(ngates),sep=''),'background') , paste('lag',seq(nlags),sep=''))

    # Strip off skipped time lags
#    laginds <- apply( ACF[["ACF"]] , FUN=function(x){ any( !is.na( x ) ) } , MARGIN = 2 )
    laginds <- which( c( LPIparam[["maxRanges"]] , rep( LPIparam[["maxRanges"]][length(LPIparam[["maxRanges"]])] , nlags ))[1:nlags] >= LPIparam[["rangeLimits"]][1] )
    ACF <- stripACF( ACF , rgates = seq( ngates ) , lags=laginds , fullCovar=LPIparam[["fullCovar"]])

    # Range gate limits
    ACF[["rangeLimits"]] <- LPIparam[["rangeLimits"]]
    names(ACF[["rangeLimits"]]) <- ""

    # Lag integration limits
    ACF[["lagLimits"]] <- LPIparam[["lagLimits"]]
    names(ACF[["lagLimits"]]) <- ""

    # Maximum ranges
    ACF[["maxRanges"]] <- LPIparam[["maxRanges"]]
    names(ACF[["maxRanges"]]) <- ""

    # Write the output list to the file
    save( ACF=ACF , file=resFile )

    # Return the file name invisibly
    invisible( resFile )

  }
