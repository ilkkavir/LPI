## file:fftws.init.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## FFT deconvolution with optimized fft length and the fast fftw library.
## Initialization function.
##
## Arguments:
##  rrange Extreme ranges to be solved c(rmin,rmax)
##  itx    A logical vector of transmitter pulse positions.
##
## Returns:
##  s     A ffts solver environment
##

fftws.init <- function( rrange , itx , nData )
{
    require(fftw)
    
    ## Minimum range
    rmin        <- min( rrange )
    
    ## Maximum range
    rmax        <- max( rrange )

    ## profile length
    lprof <- rmax - rmin
    
    ## longest pulse
    ps <- which(diff(itx>0)==1) 
    pe <- which(diff(itx>0)==-1)
    ps <- ps[ps<nData]
    pe <- ps[pe<=nData]
    ps <- ps[ps<max(pe)]
    pe <- pe[pe>min(ps)]
    plenmax <- max(pe-ps)

    ## Select the FFT length
    n <- nextn( lprof*2 + plenmax*4 )

    ## the fft plan, try with a reasonable effort (should do this earlier and only once if more effort is used)
    FFTplan <- planFFT(n,effort=1)
    
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
    assign( 'FFTplan'   , FFTplan   , s )

    # return the environment
    return( s )

  }
