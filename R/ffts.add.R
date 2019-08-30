## file:ffts.add.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## FFT deconvolution. 
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

ffts.add <- function( e , M.data , M.ambig , I.ambig , I.prod , E.data , nData )
  {
    #
    # FFT deconvolution. Data accumulation function.
    #
    # I. Virtanen 2012 
    #

    # Return immediately if the ambiguity
    # function is zero at all points
    if( ! any( I.ambig[1:nData] ) ) return()
  
    # Remove possibly remaining non-zero values
    # from points with unset index vector
    M.data[  which(!I.prod)  ] <- 0+0i
    E.data[  which(!I.prod)  ] <- 0
    M.ambig[ which(!I.ambig) ] <- 0+0i
    
    # Locate pulse start positions 
    ps <- which( diff( I.ambig[1:nData] > 0 ) == 1 )
    
    # The first point should be adjusted to pulse start,
    # so it is safe to use if the index is set
    if( I.ambig[1] ) ps <- c( 1 , ps )
    npulse <- length( ps )
    
    # Locate pulse end positions
    pe <- which( diff( I.ambig[1:nData] > 0 ) == -1 )
    
    # pe and ps should be of the same length,
    # but check anyway...
    npulse <- min( length(pe) , length(ps) )

    # Add data from one IPP at a time
    for( k in seq( npulse ) ){

      # Set temporary vectors to zero
      e[["amb.tmp"]][] <- e[["meas.tmp"]][] <- 0.+0.i
    
      # Pulse end or data end (should always be pulse end,
      # but check anyway)
      pe1              <- min( nData , pe[k] )
      
      # max range or data end
      pe2              <- min( nData , ( ps[k] + e[["n"]] - 1 ) )
    
      # Copy one pulse
      e[["amb.tmp"]][ 1 : ( pe1 - ps[k] + 1 ) ]   <- M.ambig[ ps[k] : pe1 ]
      
      # Take fft
      e[["famb.tmp"]][] <- fft( e[["amb.tmp"]] )
      
      # Copy data
      e[["meas.tmp"]][ 1 : ( pe2 - ps[k] + 1 ) ] <- M.data[ ps[k] : pe2 ]
    
      # Actual addition to the solver
      e[["fy"]][]     <- e[["fy"]]     + Conj( e[["famb.tmp"]] ) * fft( e[["meas.tmp"]] )
      e[["sqfamb"]][] <- e[["sqfamb"]] + abs( e[["famb.tmp"]] )**2
      
    }

    # Variances
    e[["varsum"]] <- e[["varsum"]] + sum( E.data[ 1 : nData ] )
    e[["nmeas"]]  <- e[["nmeas"]]  + sum( ( I.prod[ 1 : nData ] > 0 ) )
    
    invisible()
    
  }
