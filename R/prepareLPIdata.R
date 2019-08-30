## file:prepareLPIdata.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Create the LPI environment that is passed from local
## control nodes to remote nodes. A list is created instead
## of the final environment because it is faster to transfer.
##
## Arguments:
##  LPIparam        An LPI parameter list
##  LPIdatalist.raw A raw data list returned by a data input
##                  function. (See e.g. readLPIdata.gdf)
##
## Returns:
##  LPIdatalist.final The final data list that is transferred
##                     to the solver nodes.
##

prepareLPIdata <- function( LPIparam , LPIdatalist.raw )
  {
    # Internally used data vectors
    dTypes <- c( "RX1" , "RX2" , "TX1" , "TX2" )

    # An empty list for the output data
    LPIdatalist.final <- vector(mode="list",length=4)
    names(LPIdatalist.final) <- dTypes


    # A list for TX1  pulse start positions in all data
    # vectors (these will be different if sample rates
    # are different).
    # Initialise with zeros to handle data vectors without
    # pulses (they will also go through the whole system
    # and NA results will be written). The pulseStarts will
    # be passed to c-routines as such, and 0 is thus the
    # firs index.
    pulseStarts <- list( TX1 = c(0) , TX2 = c(0) , RX1 = c(0) , RX2 = c(0) )

    # A list for first sample to use in decimation
    # in each data vector
    firstSample <- c( TX1 = 0 , TX2 = 0 , RX1 = 0 , RX2 = 0 )

    # Pulse start positions in TX1 ( >0 used because
    # c-routines may have put values larger than one
    #  to the idata vector)
    pulseStarts[["TX1"]] <- which( diff( LPIdatalist.raw[["TX1"]][["idata"]][1:LPIdatalist.raw[["TX1"]][["ndata"]]] > 0 ) == 1 )

    # Calculate the corresponding pulse
    # start positions in other data vectors
    for( XXN in dTypes ){
      pulseStarts[[XXN]] <- round( as.numeric(pulseStarts[["TX1"]]) / LPIparam[["filterLength"]][["TX1"]] * LPIparam[["nup"]][["TX1"]] * LPIparam[["filterLength"]][[XXN]] / LPIparam[["nup"]][[XXN]] )
      firstSample[[XXN]] <- pulseStarts[[XXN]][1]
    }

    # The below fix does not work if 'nup' are not common for all data vectors.
    # Disable in this case.

    if(all(LPIparam[["nup"]]==LPIparam[["nup"]]["TX1"])){
      # Strip off samples to make each
      # IPP a multiple of filter length
      for( XXN in dTypes ){
        # New pulse start positions that
        # are even multiples of the filter length
        pstarts2 <- pulseStarts[[XXN]] - round( ( pulseStarts[[XXN]] - firstSample[[XXN]] ) %% ( LPIparam[["filterLength"]][[XXN]] / LPIparam[["nup"]][[XXN]] ) )

        # Do something only if the pulse positions
        # really need to be modified
        if( any( pstarts2 != pulseStarts[[XXN]] ) ){

          # Amount of shift needed in original data
          ncut <- pulseStarts[[XXN]] - pstarts2
          ntx <- length(ncut)

          # Because we are cutting off data samples,
          # the start point k-1 will already be adjusted
          # when handling point k. We will thus need to
          # subtract the number of points cut in point
          # k-1 from the original ncut[k]. Then take
          # modulus to make sure that no points will be
          # cut  unless really necessary and that number
          # of points to cut is not negative
          ncut[2:ntx] <- ncut[2:ntx] - ncut[1:(ntx-1)]
          ncut <- ncut %% round( LPIparam[["filterLength"]][[XXN]] / LPIparam[["nup"]][[XXN]] )
          ind <- rep( TRUE , LPIdatalist.raw[[XXN]][["ndata"]])
          for( k in seq(length(pstarts2)) ){
            if( ( ncut[k] > 0 ) & (pulseStarts[[XXN]][k]<LPIdatalist.raw[[XXN]][["ndata"]]) ) ind[(pulseStarts[[XXN]][k]-ncut[k]+1):pulseStarts[[XXN]][k]] <- FALSE
          }
          # Number of data points must have changed
          # as samples were cut off, update the values
          LPIdatalist.raw[[XXN]][["ndata"]] <- min( LPIdatalist.raw[[XXN]][["ndata"]] , sum(ind) )
          LPIdatalist.raw[[XXN]][["cdata"]] <- LPIdatalist.raw[[XXN]][["cdata"]][ind][1:LPIdatalist.raw[[XXN]][["ndata"]]]
          LPIdatalist.raw[[XXN]][["idata"]] <- LPIdatalist.raw[[XXN]][["idata"]][ind][1:LPIdatalist.raw[[XXN]][["ndata"]]]
        }
      }

    }

    # The idata vectors will be modified according
    # to LPIparam$indexShift before decimation.
    # Take this into account in firstSamples.
    # Again keep 0 as the first index, because
    # the indices will be passed to c-routines as such
    firstSample[["TX1"]] <- firstSample[["TX1"]] + LPIparam[["indexShifts"]][["TX1"]][1]
#    while( firstSample[["TX1"]] < 0 ){
#      firstSample[["TX1"]]  <- firstSample[["TX1"]] - LPIparam[["filterLength"]][["TX1"]] / LPIparam[["nup"]][["TX1"]]
    while( firstSample[["TX1"]] < 0 ){
      firstSample[["TX1"]]  <- firstSample[["TX1"]] + LPIparam[["filterLength"]][["TX1"]] / LPIparam[["nup"]][["TX1"]]
    }

    firstFraction <- c( TX1 = 0 , TX2 = 0 , RX1 = 0 , RX2 = 0 )
    for( XXN in dTypes ){
      firstSampleF <- firstSample[["TX1"]] * LPIparam[["filterLength"]][[XXN]] / LPIparam[["filterLength"]][["TX1"]] / LPIparam[["nup"]][[XXN]] * LPIparam[["nup"]][["TX1"]] 
      firstSample[[XXN]] <- round( firstSampleF )
      firstFraction[[XXN]] <- round( ( firstSample[[XXN]] - firstSampleF ) * LPIparam[["nup"]][[XXN]] )
    }


    # Conversions to integer mode
    storage.mode( firstSample ) <- "integer"
    storage.mode( LPIparam[["filterLength"]] ) <- "integer"
    storage.mode( firstFraction ) <- "integer"

    # Index corrections, frequency mixing,
    # and filtering in C routines
    for( XXN in dTypes ){

      storage.mode( LPIparam[["indexShifts"]][[XXN]] ) <- "integer"

      LPIdatalist.final[[XXN]] <-
        .Call( "prepare_data"                   ,
              LPIdatalist.raw[[XXN]][["cdata"]] ,
              LPIdatalist.raw[[XXN]][["idata"]] ,
              LPIdatalist.raw[[XXN]][["ndata"]] ,
              LPIparam[["freqOffset"]][XXN]     ,
              LPIparam[["indexShifts"]][[XXN]]  ,
              LPIparam[["nup"]][XXN]            ,
              LPIparam[["filterLength"]][XXN]   ,
              firstSample[[XXN]]                ,
              firstFraction[[XXN]]              ,
              TRUE
              )

    }

    # Use length of the shortest data vector
    LPIdatalist.final[["nData"]] <-
      min(
          LPIdatalist.final[["RX1"]][["ndata"]],
          LPIdatalist.final[["RX2"]][["ndata"]],
          LPIdatalist.final[["TX1"]][["ndata"]],
          LPIdatalist.final[["TX2"]][["ndata"]]
          )


    # Optional TX amplitude normalisation
    if( LPIparam[["normTX"]] ){
      itx1 <- which(LPIdatalist.final[["TX1"]][["idata"]][1:LPIdatalist.final[["nData"]]])
      itx2 <- which(LPIdatalist.final[["TX2"]][["idata"]][1:LPIdatalist.final[["nData"]]])
      txamp1 <- mean(abs(LPIdatalist.final[["TX1"]][["cdata"]][itx1]))
      txamp2 <- mean(abs(LPIdatalist.final[["TX2"]][["cdata"]][itx2]))
      LPIdatalist.final[["TX1"]][["cdata"]][itx1] <- exp(1i*Arg(LPIdatalist.final[["TX1"]][["cdata"]][itx1])) * txamp1
      LPIdatalist.final[["TX2"]][["cdata"]][itx2] <- exp(1i*Arg(LPIdatalist.final[["TX2"]][["cdata"]][itx2])) * txamp2
    }

    # Optional ground clutter suppression
    if( ( LPIparam[["maxClutterRange"]]["RX1"] > 0 ) & ( LPIparam[["clutterFraction"]][["RX1"]] > 0 )){
      clutterSuppress( LPIdatalist.final[["TX1"]] , LPIdatalist.final[["RX1"]] , LPIparam[["rangeLimits"]][1] , LPIparam[["maxClutterRange"]]["RX1"] , LPIdatalist.final[["nData"]] , LPIparam[["clutterFraction"]][["RX1"]] )
    }
    if( ( LPIparam[["maxClutterRange"]]["RX2"] > 0 ) & ( LPIparam[["clutterFraction"]][["RX2"]] > 0 )){
      clutterSuppress( LPIdatalist.final[["TX2"]] , LPIdatalist.final[["RX2"]] , LPIparam[["rangeLimits"]][1] , LPIparam[["maxClutterRange"]]["RX2"] , LPIdatalist.final[["nData"]] , LPIparam[["clutterFraction"]][["RX2"]] )
    }


    # Optional voltage level decoding
    if( is.numeric( LPIparam[["decodingFilter"]] ) ){

        LPIdatalist.final[["RX1"]][["cdata"]][!LPIdatalist.final[["RX1"]][["idata"]]] <- 0+0i
        LPIdatalist.final[["RX2"]][["cdata"]][!LPIdatalist.final[["RX2"]][["idata"]]] <- 0+0i
        LPIdatalist.final[["TX1"]][["cdata"]][!LPIdatalist.final[["TX1"]][["idata"]]] <- 0+0i
        LPIdatalist.final[["TX2"]][["cdata"]][!LPIdatalist.final[["TX2"]][["idata"]]] <- 0+0i

        nd <- LPIdatalist.final[["nData"]]


        LPIdatalist.final[["RX1"]][["cdata"]] <- LPI:::decoFilter.cdata( LPIdatalist.final[["RX1"]][["cdata"]][1:nd] , LPIdatalist.final[["TX1"]][["cdata"]][1:nd] , LPIdatalist.final[["TX1"]][["idata"]][1:nd] , LPIparam[["decodingFilter"]][1] )

        LPIdatalist.final[["TX1"]][["cdata"]] <- LPI:::decoFilter.cdata( LPIdatalist.final[["TX1"]][["cdata"]][1:nd] , LPIdatalist.final[["TX1"]][["cdata"]][1:nd] , LPIdatalist.final[["TX1"]][["idata"]][1:nd] , LPIparam[["decodingFilter"]][1] )

        LPIdatalist.final[["RX2"]][["cdata"]] <- LPI:::decoFilter.cdata( LPIdatalist.final[["RX2"]][["cdata"]][1:nd] , LPIdatalist.final[["TX2"]][["cdata"]][1:nd] , LPIdatalist.final[["TX2"]][["idata"]][1:nd] , LPIparam[["decodingFilter"]][1] )

        LPIdatalist.final[["TX2"]][["cdata"]] <- LPI:::decoFilter.cdata( LPIdatalist.final[["TX2"]][["cdata"]][1:nd] , LPIdatalist.final[["TX2"]][["cdata"]][1:nd] , LPIdatalist.final[["TX2"]][["idata"]][1:nd] , LPIparam[["decodingFilter"]][1] )

    }else if( is.character( LPIparam[["decodingFilter"]] )){

        if( any( LPIparam[["decodingFilter"]][1] == c("matched","inverse") ) ){

            LPIdatalist.final[["RX1"]][["cdata"]][!LPIdatalist.final[["RX1"]][["idata"]]] <- 0+0i
            LPIdatalist.final[["RX2"]][["cdata"]][!LPIdatalist.final[["RX2"]][["idata"]]] <- 0+0i
            LPIdatalist.final[["TX1"]][["cdata"]][!LPIdatalist.final[["TX1"]][["idata"]]] <- 0+0i
            LPIdatalist.final[["TX2"]][["cdata"]][!LPIdatalist.final[["TX2"]][["idata"]]] <- 0+0i

            nd <- LPIdatalist.final[["nData"]]


            LPIdatalist.final[["RX1"]][["cdata"]] <- LPI:::decoFilter.cdata( LPIdatalist.final[["RX1"]][["cdata"]][1:nd] , LPIdatalist.final[["TX1"]][["cdata"]][1:nd] , LPIdatalist.final[["TX1"]][["idata"]][1:nd] , LPIparam[["decodingFilter"]][1] )

            LPIdatalist.final[["TX1"]][["cdata"]] <- LPI:::decoFilter.cdata( LPIdatalist.final[["TX1"]][["cdata"]][1:nd] , LPIdatalist.final[["TX1"]][["cdata"]][1:nd] , LPIdatalist.final[["TX1"]][["idata"]][1:nd] , LPIparam[["decodingFilter"]][1] )

            LPIdatalist.final[["RX2"]][["cdata"]] <- LPI:::decoFilter.cdata( LPIdatalist.final[["RX2"]][["cdata"]][1:nd] , LPIdatalist.final[["TX2"]][["cdata"]][1:nd] , LPIdatalist.final[["TX2"]][["idata"]][1:nd] , LPIparam[["decodingFilter"]][1] )

            LPIdatalist.final[["TX2"]][["cdata"]] <- LPI:::decoFilter.cdata( LPIdatalist.final[["TX2"]][["cdata"]][1:nd] , LPIdatalist.final[["TX2"]][["cdata"]][1:nd] , LPIdatalist.final[["TX2"]][["idata"]][1:nd] , LPIparam[["decodingFilter"]][1] )

            LPIdatalist.final[["TX1"]][["idata"]] <- LPI:::decoFilter.idata( LPIdatalist.final[["TX1"]][["idata"]][1:nd] )

            LPIdatalist.final[["TX2"]][["idata"]] <- LPI:::decoFilter.idata( LPIdatalist.final[["TX2"]][["idata"]][1:nd] )
        }
    }


    # Largest range in rangeLimits
    maxr <- as.integer(max(LPIparam[["rangeLimits"]]))

    # Average signal powers, loop three times in order to make simple noise spike detection as well
    for(niter in seq(3)){

      # Average power in signal vector RX1
      LPIdatalist.final[["RX1"]][["power"]] <- LPIaveragePower( LPIdatalist.final[["RX1"]][["cdata"]] , LPIdatalist.final[["TX1"]][["idata"]] , LPIdatalist.final[["RX1"]][["idata"]] , LPIdatalist.final[["nData"]] , maxr )

      # Average power in signal vector RX2
      LPIdatalist.final[["RX2"]][["power"]] <- LPIaveragePower( LPIdatalist.final[["RX2"]][["cdata"]] , LPIdatalist.final[["TX2"]][["idata"]] , LPIdatalist.final[["RX2"]][["idata"]] , LPIdatalist.final[["nData"]] , maxr )

      # Flag data points whose power is more than four times the average at a given height,
      # but only if there were reasonably many samples in the averages
      if(LPIdatalist.final[["RX1"]][["power"]][1] < .05 ){
          itx1 <- which( abs(LPIdatalist.final[["RX1"]][["cdata"]][1:LPIdatalist.final[["nData"]]]) > (sqrt(LPIdatalist.final[["RX1"]][["power"]])*4) )
          LPIdatalist.final[["RX1"]][["idata"]][itx1] <- FALSE
      }
      if(LPIdatalist.final[["RX2"]][["power"]][1] < .05 ){
          itx2 <- which( abs(LPIdatalist.final[["RX2"]][["cdata"]][1:LPIdatalist.final[["nData"]]]) > (sqrt(LPIdatalist.final[["RX2"]][["power"]])*4) )
          LPIdatalist.final[["RX2"]][["idata"]][itx1] <- FALSE
      }
    }

    # maxr points in the beginning will not have
    # a reasonable  power estimate, flag these points as well
    LPIdatalist.final[["RX1"]][["idata"]][1:maxr] <- FALSE
    LPIdatalist.final[["RX2"]][["idata"]][1:maxr] <- FALSE

    ######################################
    ## Copy parameters from LPIparam to ##
    ## the final data list as necessary ##
    ######################################

    # Lag values
    LPIdatalist.final[["lagLimits"]] <- LPIparam[["lagLimits"]]
    LPIdatalist.final[["nLags"]]     <- length(LPIdatalist.final[["lagLimits"]]) - 1


    # Maximum ranges, repeat the last value as necessary
    LPIdatalist.final[["maxRanges"]] <-  LPIparam[["maxRanges"]]
    nmaxr <- length(LPIdatalist.final[["maxRanges"]])
    if( nmaxr < LPIdatalist.final[["nLags"]] ){
      LPIdatalist.final[["maxRanges"]] <- c( LPIdatalist.final[["maxRanges"]] , rep(LPIdatalist.final[["maxRanges"]][nmaxr],(LPIdatalist.final[["nLags"]]-nmaxr)))
    }




    # Range gate limits
    LPIdatalist.final[["rangeLimits"]] <- LPIparam[["rangeLimits"]]
    LPIdatalist.final[["nGates"]] <- rep( length(LPIparam[["rangeLimits"]]) - 1 , LPIdatalist.final[["nLags"]] )
    for( k in seq(LPIdatalist.final[["nLags"]]) ){
      LPIdatalist.final[["nGates"]][k] <- length( LPIparam[["rangeLimits"]][ LPIparam[["rangeLimits"]] < LPIdatalist.final[["maxRanges"]][k] ] ) - 1
    }

    # The TX vectors are always decimated
    # in the present version
    LPIdatalist.final[["nDecimTX"]] <- 1

    # Number of theory matrix rows to buffer
    LPIdatalist.final[["nBuf"]] <- LPIparam[["nBuf"]]

    # Inverse problem solver
    LPIdatalist.final[["solver"]] <- LPIparam[["solver"]]

    # Options to rlips
    LPIdatalist.final[["rlips.options"]] <- LPIparam[["rlips.options"]]

    # Do we calculate background ACF estimates
    LPIdatalist.final[["backgroundEstimate"]] <- LPIparam[["backgroundEstimate"]]

    # Should full covariance matrix or only its
    # diagonal be calculated
    LPIdatalist.final[["fullCovar"]] <- LPIparam[["fullCovar"]]

    # Are we running in a cluster or locally
    LPIdatalist.final[["iscluster"]] <- LPIparam[["iscluster"]]

    # Is the rx data from a remote site?
    LPIdatalist.final[["remoteRX"]] <- LPIparam[["remoteRX"]]

    # Number of codes if pre-averaging is being used
    LPIdatalist.final[["nCode"]] <- LPIparam[["nCode"]]

    # Should interpolation be used when calculating
    # the range ambiguity functions
    LPIdatalist.final[["ambInterp"]] <- LPIparam[["ambInterp"]]

    # Make sure that the storage modes are correct
    storage.mode(LPIdatalist.final[["TX1"]][["cdata"]])  <- "complex"
    storage.mode(LPIdatalist.final[["TX2"]][["cdata"]])  <- "complex"
    storage.mode(LPIdatalist.final[["TX1"]][["idata"]])  <- "logical"
    storage.mode(LPIdatalist.final[["TX2"]][["idata"]])  <- "logical"
    storage.mode(LPIdatalist.final[["RX1"]][["cdata"]])  <- "complex"
    storage.mode(LPIdatalist.final[["RX2"]][["cdata"]])  <- "complex"
    storage.mode(LPIdatalist.final[["RX1"]][["idata"]])  <- "logical"
    storage.mode(LPIdatalist.final[["RX2"]][["idata"]])  <- "logical"
    storage.mode(LPIdatalist.final[["RX1"]][["power"]])  <- "double"
    storage.mode(LPIdatalist.final[["RX2"]][["power"]])  <- "double"
    storage.mode(LPIdatalist.final[["lagLimits"]])       <- "integer"
    storage.mode(LPIdatalist.final[["rangeLimits"]])     <- "integer"
    storage.mode(LPIdatalist.final[["nDecimTx"]])        <- "integer"
    storage.mode(LPIdatalist.final[["nBuf"]])            <- "integer"
    storage.mode(LPIdatalist.final[["nData"]])           <- "integer"
    storage.mode(LPIdatalist.final[["nGates"]])          <- "integer"
    storage.mode(LPIdatalist.final[["nLags"]])           <- "integer"
    storage.mode(LPIdatalist.final[["nCode"]])           <- "integer"
    storage.mode(LPIdatalist.final[["ambInterp"]])       <- "logical"
    storage.mode(LPIdatalist.final[["backgroundEstimate"]]) <- "logical"

    return( LPIdatalist.final )

  }
