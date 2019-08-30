## file:noUpdate.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## LPI parameter list update function
## Return the list itself in first call,
## NULL in the second call with the same list
##
## Arguments:
##  LPIparam  A LPI parameter list
##  intPeriod Integration period number
##
## Returns:
##  LPIparam An exact copy  of the input LPIparam in
##           first call, NULL in second call with the
##           same list
##

noUpdate <-  function( LPIparam , intPeriod )
    {

        if(is.null(LPIparam[["callN"]])){
            LPIparam[["callN"]] <- 1
            return(LPIparam)
        }

        return(NULL)

    }
        
