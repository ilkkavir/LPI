## file:decoFilter.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Voltage level decoding, either matched or inverse
## filtering, using measured transmitter samples.
##
## Arguments:
##  cdata       A complex receiver data vector
##  cenv        A complex transmitter data vector
##  idata       A logical vector of transmitter data indices
##  filterType  Decoding filter. Either a complex vector of filter taps, 'matched' or 'inverse'
##
## Returns:
##  cdata The complex receiver data vector after decoding
##

decoFilter.cdata <- function( cdata , cenv , idata , filterType='inverse')
  {
    # Pulse start positions and number of pulses
    txstarts <- which( diff(idata>0) == 1 )
    if(idata[1]) txstarts <- c(0,txstarts)
    ntx <- length(txstarts)
    txstarts <- c(txstarts,length(cdata))


    # If there are no transmission pulses, then simply return
    if(ntx<1) return(cdata)

    # Set the data points before the first pulse to zero
    if( txstarts[1] > 0 ) cdata[1:txstarts[1]] <- 0+0i


    # Set transmitter data to zero at points that are not transmitter samples
    cenv[!idata] <- 0+0i

    # Filtering with user-defined coefficients
    if( is.numeric( filterType ) ){
        nfilter <- length(filterType)
        for( k in seq(  ntx ) ){
            cenv[ (txstarts[k]+1) : (txstarts[k+1]) ] <- 0+0i
            cenv[ (txstarts[k]+1) :  (txstarts[k]+nfilter)] <- filterType
            cdata[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                fft(
                    fft( cdata[ (txstarts[k]+1) : (txstarts[k+1]) ] ) /
                    fft(  cenv[ (txstarts[k]+1) : (txstarts[k+1]) ] )
                    , inverse=TRUE ) /
                        (txstarts[k+1]-txstarts[k]) * sqrt(sum(abs(cenv[ (txstarts[k]+1) : (txstarts[k+1]) ])**2))
        }        
    }else if( is.character( filterType ) ){
        # Inverse filtering
        if(filterType=="inverse"){
            for( k in seq(  ntx ) ){
                cdata[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                    fft(
                        fft( cdata[ (txstarts[k]+1) : (txstarts[k+1]) ] ) /
                        fft(  cenv[ (txstarts[k]+1) : (txstarts[k+1]) ] )
                        , inverse=TRUE ) /
                            (txstarts[k+1]-txstarts[k]) * sqrt(sum(abs(cenv[ (txstarts[k]+1) : (txstarts[k+1]) ])**2))
            }
        # Matched filtering
        }else if(filterType=="matched"){
            for( k in seq( ntx ) ){
                cdata[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                    fft(
                        fft( cdata[ (txstarts[k]+1) : (txstarts[k+1]) ] ) *
                        Conj( fft( cenv[ (txstarts[k]+1) : (txstarts[k+1]) ] ) )
                        , inverse=TRUE ) /
                            (txstarts[k+1]-txstarts[k]) / sqrt(sum(abs(cenv[ (txstarts[k]+1) : (txstarts[k+1]) ])**2))
            }
        # Other filters are not supported at the moment
        }else{
            stop("Unknown decoding filter")
        }
    }else{
        stop("Unknown decoding filter")
    }

  return(cdata)

}

##
## Index corrections for decoded receiver data
## 
##
## Arguments:
##  idata A logical vector of transmitter data indices
##
## Returns:
##  idata A corrected index vector with only first index
##        of each pulse set.
##

decoFilter.idata <- function( idata )
  {

    # Pulse start positions
    txstarts <- which( diff(idata>0) == 1 )
    if(idata[1]) txstarts <- c(0,txstarts)
    ntx <- length(txstarts)
    txstarts <- c(txstarts,length(idata))
    
    # Each pulse should have been compressed into
    # a single sample in the decoding
    for( k in seq( ntx ) ){
      idata[(txstarts[k]+2):txstarts[k+1]] <- FALSE
    }
    
    return(idata)
    
}
