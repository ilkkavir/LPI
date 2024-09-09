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
##  cdataT        A complex transmitter data vector
##  idataT      A logical vector of transmitter data indices
##  idataR      A logical vector of accepted RX data indices
##  filterType  Decoding filter. Either a complex vector of filter taps, 'matched' or 'inverse'
##
## Returns:
##  cdata The complex receiver data vector after decoding
##

##
## this is not working, do not use! Will replace with something better... 
##


decoFilter2 <- function( cdataT , idataT , cdataR , idataR , ndata , filterType='inverse'){
    ## Pulse start positions and number of pulses
    txstarts <- which( diff(idataT>0) == 1 )
    if(idataT[1]) txstarts <- c(0,txstarts)
    txstarts <- txstarts[txstarts<ndata]
    ntx <- length(txstarts)    
    txstarts <- c(txstarts,ndata)

    
    ## If there are no transmission pulses, then simply return
    if(ntx<1) return(cdataR)
    
    ## Set the data points before the first pulse to zero
    if( txstarts[1] > 0 ) cdataR[1:txstarts[1]] <- 0+0i
    
    
    ## Set transmitter data to zero at points that are not transmitter samples
    cdataT[!idataT] <- 0+0i

    ## the same for the samples of the received signal
    cdataR[!idataR] <- 0+0i

    ## make sure that there are only zeros and ones in idataR
    idataR = idataR!=0
    
    ## Filtering with user-defined coefficients
    if( is.numeric( filterType ) ){
        nfilter <- length(filterType)
        for( k in seq(  ntx ) ){
            cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] <- 0+0i
            cdataT[ (txstarts[k]+1) :  (txstarts[k]+nfilter)] <- filterType
            cpow <- abs(cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ])**2
            dscale <- sqrt(fft( fft( cpow ) * Conj( fft( idataR[ (txstarts[k]+1) : (txstarts[k+1]) ] ) ) , inverse=TRUE ) / (txstarts[k+1]-txstarts[k]) )
            cdataR[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                fft(
                    fft( cdataR[ (txstarts[k]+1) : (txstarts[k+1]) ] ) *
                    Conj( fft(  cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) ) , inverse=TRUE ) / ((txstarts[k+1]-txstarts[k]) * dscale)
            cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                fft(
                    fft( cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) *
                    Conj( fft(  cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) ) , inverse=TRUE ) / ((txstarts[k+1]-txstarts[k]) * sqrt(sum(abs(cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ])**2)) )
            idataR[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                abs( fft(
                    fft( idataR[ (txstarts[k]+1) : (txstarts[k+1]) ] ) *
                    Conj( fft(  idataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) ) , inverse=TRUE ) / ((txstarts[k+1]-txstarts[k]) ) 
                    ) > .1
        }        
    }else if( is.character( filterType ) ){
        ## Inverse filtering
        if(filterType=="inverse"){
            for( k in seq(  ntx ) ){
                cpow <- abs(cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ])**2
                dscale <- sqrt(fft( fft( cpow ) * Conj( fft( idataR[ (txstarts[k]+1) : (txstarts[k+1]) ] ) ) , inverse=TRUE ) / (txstarts[k+1]-txstarts[k]) )
                cdataR[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                    fft(
                        fft( cdataR[ (txstarts[k]+1) : (txstarts[k+1]) ] ) /
                        fft(  cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) , inverse=TRUE ) / ( (txstarts[k+1]-txstarts[k]) * dscale )
                cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                    fft(
                        fft( cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) /
                        fft(  cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) , inverse=TRUE ) / ( (txstarts[k+1]-txstarts[k]) * sqrt(sum(abs(cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ])**2)) )
                idataR[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                    abs( fft(
                        fft( idataR[ (txstarts[k]+1) : (txstarts[k+1]) ] ) *
                        Conj( fft(  idataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) ) , inverse=TRUE ) / ((txstarts[k+1]-txstarts[k]) ) 
                        ) > .1
                ## each pulse should has been compressed into a single short pulse in the decoding
                idataT[(txstarts[k]+2):txstarts[k+1]] <- FALSE
            }
        # Matched filtering
        }else if(filterType=="matched"){
            for( k in seq( ntx ) ){

                cpow <- abs(cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ])**2
                dscale <- abs(sqrt(fft( fft( cpow ) * Conj( fft( idataR[ (txstarts[k]+1) : (txstarts[k+1]) ] ) ) , inverse=TRUE ) / (txstarts[k+1]-txstarts[k]) ))
                cdataR[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                    fft(
                        fft( cdataR[ (txstarts[k]+1) : (txstarts[k+1]) ] ) *
                        Conj( fft( cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) ) , inverse=TRUE ) / ( (txstarts[k+1]-txstarts[k]) * dscale )
                cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                    fft(
                        fft( cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) *
                        Conj( fft( cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) ) , inverse=TRUE ) / ( (txstarts[k+1]-txstarts[k]) * sqrt(sum(abs(cdataT[ (txstarts[k]+1) : (txstarts[k+1]) ])**2)) )
                idataR[ (txstarts[k]+1) : (txstarts[k+1]) ] <-
                    abs( fft(
                        fft( idataR[ (txstarts[k]+1) : (txstarts[k+1]) ] ) *
                        Conj( fft( idataT[ (txstarts[k]+1) : (txstarts[k+1]) ] ) ) , inverse=TRUE ) / ((txstarts[k+1]-txstarts[k]) ) 
                        ) > .1
                ## each pulse should has been compressed into a single short pulse in the decoding
                idataT[(txstarts[k]+2):txstarts[k+1]] <- FALSE

            }
            ## Other filters are not supported at the moment
        }else{
            stop("Unknown decoding filter")
        }
    }else{
        stop("Unknown decoding filter")
    }
    
    return(list(cdataT=cdataT , idataT=idataT , cdataR=cdataR , idataR=idataR))

    
}

