## file:currentTimes.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.

##
## Current unix time minus 5 seconds to be used for
## identifying the latest available data samples in
## real time analysis. 
##
## Arguments:
##   ...  An arbitrary list of arguments is accepted, but none
##        of them will be used.
##
## Returns:
##  curTimes A named vector ("TX1","TX2","RX1","RX2") with
##           the current unix time -5 in each element.
##
## 

currentTimes <- function( ... )
  {
    return( LPIexpand.input( as.numeric(Sys.time()-5) ) )
  }
