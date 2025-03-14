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
## only matched filtering with recorded waveforms is implemented at the moment
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
    if( txstarts[1] > 0 ){
        cdataR[1:txstarts[1]] <- 0+0i
        cdataT[1:txstarts[1]] <- 0+0i
    }
    
    
    ## Set transmitter data to zero at points that are not transmitter samples
    cdataT[!idataT] <- 0+0i

    ## the same for the samples of the received signal
    cdataR[!idataR] <- 0+0i

    ## make sure that there are only zeros and ones in idataR
    idataR = idataR!=0

    ## Filtering with user-defined coefficients
    if( is.numeric( filterType ) ){
        stop("the filter is not implemented in this version")
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
            stop("the inverse filter is not implemented in this version")
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
                
                ## use IPP + pulse length samples in decoding to avoid data loss in bistatic measurements
                ii1 <- txstarts[k]+1 # first sample
                ii2 <- txstarts[k+1] # last sample to decode
                nn2 <- ii2-ii1+1 # number of samples to decode
                plen <- sum(idataT[ ii1 : ii2 ])## pulse length
                ii3 <- txstarts[k] + nextn(txstarts[k+1] - txstarts[k] + plen) # last sample in fft
                ii3 <- min(ii3,ndata) ## will this cause an error after the last pulse?
                nii <- ii3-ii1+1 # number of samples in fft

                ## zero-padded transmission envelope
                cdataTtmp <- c( cdataT[ ii1 : ii2 ]  , rep(0+0i , ii3-ii2) )
                idataTtmp <- c( idataT[ ii1 : ii2 ]  , rep( FALSE , ii3-ii2) )

                ## power of the complex transmission envelope
                cpow <- abs(cdataTtmp)**2

                ## complex conjugate of fft of the transmission envelope
                cdataTfftConj <- Conj( fft(cdataTtmp) ) 

                ## data scaling factors
                dscale <- sqrt( abs( fft( Conj( fft( cpow) ) * fft( idataR[ ii1 : ii3 ] ) , inverse=TRUE ) / nii ))[ 1:nn2 ]

                ## decode the complex received signal
                cdataR[ ii1 : ii2] <- fft( fft( cdataR[ ii1 : ii3 ] ) * cdataTfftConj , inverse=TRUE )[ 1:nn2 ] / ( nii * dscale )

                ## decode the complex transmitted signal
                cdataT[ ii1 : ii2] <- fft( abs(cdataTfftConj)**2 , inverse=TRUE )[ 1:nn2 ] / ( nii * sqrt(sum( cpow )) )

                ## fix the RX indices
#                idataR[ ii1 : ii2 ]  <-  abs( fft( fft(idataR[ ii1 : ii3 ] ) * Conj( fft( idataTtmp) ) , inverse=TRUE)[ 1:nn2 ]  / nii ) > .1
                idataR[ ii1 : ii2 ]  <-  dscale/sqrt(sum(cpow)) > .1

                ## each pulse should have been compressed into a single short pulse in the decoding
                idataT[ (ii1 + 1 ):ii2] <- FALSE

            }
            ## Other filters are not supported at the moment
        }else{
            stop("Unknown decoding filter")
        }
    }else{
        stop("Unknown decoding filter")
    }
    
    return(list(cdataT=cdataT , idataT=idataT , cdataR=cdataR , idataR=idataR ))

    
}

