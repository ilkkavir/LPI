## file:nextIntegrationPeriods.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.

##
## Indices of n latest integration periods
## that have not yet been analysed.
##
## Arguments:
##  LPIparam     A LPI parameter list
##  n            Number of new periods to search for
##  intPer.ready A list of solved period indices
##
## Returns:
##  nextIpers   Indices of the integration periods to
##              be solved next.
## 

nextIntegrationPeriods <- function( LPIparam , n , intPer.missing )
  {


    # Truly available periods
    intPer.available <- intPer.missing[ which( intPer.missing <= min( LPIparam[["maxIntPeriod"]] , LPIparam[["lastIntPeriod"]]) ) ]

    # We know that the integration periods are in order,
    # simply pick the n last ones
    nper <- length(intPer.available)
    if(nper==0) return(NULL)
    return(intPer.available[ max(1,( nper - n + 1 )) : nper ])

    
##    # A vector for the integration period numbers
##    nextIpers <- rep(0,n)
##
##    # Counter for identified new periods
##    k <- 0
##    
##    # The period from which we will start seeking backwards
##    p <- min( LPIparam[["maxIntPeriod"]] , LPIparam[["lastIntPeriod"]] )
##    
##    # If the last data sample or analysis end time is before
##    # beginning of analysis, there will be nothing to do
##    if( p < 0 ) return(NULL)
##    # Start looking backwards from the last period
##    while( k < n ){
##      # Select periods that have not yet been analysed.
##      if(!any(intPer.ready == p)){
##        k <- k+1
##        nextIpers[k] <- p
##      }
##      # Stop looking if we hit the analysis start time
##      if( p == 1) break
##      p <- p - 1
##    }
##    
##    # Return NULL if nothing was found
##    if( k== 0 ) return(NULL)
##
##    return( nextIpers[1:k] )
    
  }
