## file:LPIaveragePower.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Average power profiles
##
## Arguments:
##  cdata    A complex data vector
##  idatatx  A logical vector of transmitter pulse positions
##  idatarx  A logical vector of usable receiver samples
##  ndata    Number points in data vectors
##  maxrange Largest range from which the power is needed
##  nmin     Minimum number of samples to average
##
## Returns:
##  pdata   Average power profile vector
##

LPIaveragePower <- function( cdata , idatatx , idatarx , ndata , maxrange , nmin)
  {
    # Call the C function
    pow <- .Call( "average_power" , cdata , idatatx , idatarx , ndata , maxrange , nmin )

    # Check the first element, .1 means that number of
    # summed power values is 10 in average.
    # The first element will be NA if no pulses were found,
    # then it does not really matter  what we do..
    if( is.na( pow[1] ) ){
      pow[] <-  mean( abs( cdata[idatarx])**2 )
    }else if( pow[1] > 1/nmin ){
      pow[] <-  mean( abs( cdata[idatarx])**2 )
    }

    return(pow)
  }
