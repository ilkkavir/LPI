## file:fftws.add.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## FFT deconvolution with optimized fft length and the fast fftw library.
## Data accumulation function.
##
## Arguments:
##  e       An ffts solver environemnt
##  M.data  Measurement vector
##  M.ambig Range ambiguity function
##  I.ambig Indices of non-zero ambiguity values
##  I.prod  Indices of usable lagged products
##  E.data  Measurement variance vector
##  nData   Number of points in data vectors
## 
## Returns:
##  success TRUE if the data was successfully added
##

fftws.add <- function( e , M.data , M.ambig , I.ambig , I.prod , E.data , nData )
{
    ##
    ## FFT deconvolution. Data accumulation function.
    ##
    ## I. Virtanen 2012, 2025
    ##

    ## Remove possibly remaining non-zero values
    ## from points with unset index vector
    M.data[  !I.prod  ] <- 0+0i
#    E.data[  !I.prod  ] <- 0 # no need for this, because there is an error for every single point
    M.ambig[ !I.ambig ] <- 0+0i

    dambig <- diff(I.ambig > 0)
    
    ## Locate pulse start positions 
    ps <- which( dambig == 1 )
    
    ## The first point should be adjusted to pulse start,
    ## so it is safe to use if the index is set
    if( I.ambig[1] ) ps <- c( 1 , ps )
    npulse <- length( ps )
    
    ## Locate pulse end positions
    pe <- which( dambig == -1 )
    
    ## pe and ps should be of the same length,
    ## but check anyway...
    npulse <- min( length(pe) , length(ps) )

    ## return if no pulses found
    if( npulse < 1 ) return() 
    
    ## Add data from one IPP at a time
    for( k in seq( npulse ) ){
        
        ## Set temporary vectors to zero
        e[["amb.tmp"]][] <- e[["meas.tmp"]][] <- 0.+0.i
    
        ## Pulse end or data end (should always be pulse end,
        ## but check anyway)
        pe1              <- min( nData , pe[k] )
        
        ## Copy one pulse
        e[["amb.tmp"]][ 1 : ( pe1 - ps[k] + 1 ) ]   <- M.ambig[ ps[k] : pe1 ]
        
        ## Take fft
        e[["famb.tmp"]][] <- fftw::FFT( e[["amb.tmp"]] , plan=e[["FFTplan"]])
        
        ## first echo sample to use
        s1 <- min( nData , ps[k] + e[["rmin"]] )

        ## last echo sample to use
        s2 <- min( nData , e[["rmax"]] + pe[k] )

        ## Copy data
        e[["meas.tmp"]][ 1 : ( s2 - s1 + 1 ) ] <- M.data[ s1 : s2 ]
    
        ## Actual addition to the solver
        e[["fy"]][]     <- e[["fy"]]     + Conj( e[["famb.tmp"]] ) * fftw::FFT( e[["meas.tmp"]] , plan=e[["FFTplan"]] )
        e[["sqfamb"]][] <- e[["sqfamb"]] + abs( e[["famb.tmp"]] )**2

        ## Variances
        e[["varsum"]] <- e[["varsum"]] + sum( E.data[ s1 : s2 ] )
        e[["nmeas"]]  <- e[["nmeas"]]  + sum( ( I.prod[ s1 : s2 ] > 0 ) )

    }
    
    
    invisible()
    
  }
