## file:ffts.init.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## FFT deconvolution.
## Initialization function.
##
## Arguments:
##  rrange Extreme ranges to be solved c(rmin,rmax)
##  itx    A logical vector of transmitter pulse positions.
##
## Returns:
##  s     A ffts solver environment
##

ffts.init <- function( rrange , itx )
  {
    # Minimum range
    rmin        <- min( rrange )

    # Maximum range
    rmax        <- max( rrange )

    # longest inter-pulse period
    ippmax <- max( diff( which( diff( itx > 0 ) == 1 ) ) , showWarnings=FALSE )

    # Select the FFT length
    n <- max( nextn( ippmax ) , nextn( rmax*2 ) )

    # Allocate vectors
    fy          <- rep( 0+0i , n )
    amb.tmp     <- famb.tmp         <- rep( 0+0i , n )
    meas.tmp    <- rep( 0+0i , n )
    sqfamb      <- rep( 0 , n )
    varsum      <- 0
    nmeas       <- 0

    # Set storage modes
    storage.mode( rmin )     <- "integer"
    storage.mode( rmax )     <- "integer"
    storage.mode( n )        <- "integer"
    storage.mode(nmeas)      <- "integer"
    storage.mode( fy )       <- "complex"
    storage.mode( amb.tmp )  <- "complex"
    storage.mode( famb.tmp ) <- "complex"
    storage.mode( meas.tmp ) <- "complex"
    storage.mode( sqfamb )   <- "double"
    storage.mode(varsum)     <- "double"

    # Create a new environment and assign everything to it
    s <- new.env()
    assign( 'n'         , n         , s )
    assign( 'rmin'      , rmin      , s )
    assign( 'rmax'      , rmax      , s )
    assign( 'fy'        , fy        , s )
    assign( 'sqfamb'    , sqfamb    , s )
    assign( 'amb.tmp'   , amb.tmp   , s )
    assign( 'famb.tmp'  , famb.tmp  , s )
    assign( 'meas.tmp'  , meas.tmp  , s )
    assign( 'nmeas'     , nmeas     , s )
    assign( 'varsum'    , varsum    , s )

    # return the environment
    return( s )

  }
