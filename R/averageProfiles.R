## file:averageProfiles.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Lag-profile pre-averaging before the actual inversion.
## Provides significant speed-up but may lead to somewhat
## reduced estimation accuracy
##
## This routine is intended to be used in real-time
## analysis with limited computing resources when speed
## gain with reduced accuaracy and flexibility is accepable.
##
##
## Arguments:
##  LPIenv A LPI environment
##  l      Lag number
##
## Returns:
##  success TRUE if both lagged products and range ambiguity
##           functions were successfully averaged.
##
## The averaged profiles are overwritten to
## LPIenv[["cprod"]] and LPIenv[["camb"]]
##

averageProfiles <- function( LPIenv , l )
  {

    s1 <- .Call( "average_profile" , LPIenv[["cprod"]] , LPIenv[["TX1"]][["idata"]] , as.integer( LPIenv[["nData"]] - l ) , as.integer( LPIenv[["nCode"]] ) )

    s2 <- .Call( "average_profile" , LPIenv[["camb"]]  , LPIenv[["TX1"]][["idata"]] , as.integer( LPIenv[["nData"]] - l ) , as.integer( LPIenv[["nCode"]] ) )
    
    invisible( ( s1 & s2 ) )
    
  }
