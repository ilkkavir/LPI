## file:rangeAmbiguity.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Calculation of range ambiguity functions.
##
## Arguments:
##  LPIenv  A LPI environment
##  lag     Lag number
##
##
## Returns:
##  success  TRUE if at least one point was successfully
##           calculated, FALSE otherwise.
##           The range ambiguity function is
##           (over)written to LPIenv$camb.
##
##
##
##
##

rangeAmbiguity <- function( LPIenv , lag )
  {
  
    # True oversampling is not supported.
    if( LPIenv[['nDecimTX']] != 1) stop("True transmitter signal oversampling is not supported.")

    # Make sure that lag is an integer
    storage.mode(lag) <- "integer"

    # Simulate oversampling by means of interpolation.
    # This works well if the pulses have
    # sharp edges and constant amplitude.
    if( LPIenv[["ambInterp"]] ){
      return( .Call( "range_ambiguity"          ,
                    LPIenv[["TX1"]][["cdata"]] ,
                    LPIenv[["TX2"]][["cdata"]] ,
                    LPIenv[["TX1"]][["idata"]] ,
                    LPIenv[["TX2"]][["idata"]] ,
                    LPIenv[["camb"]]           ,
                    LPIenv[["iamb"]]           ,
                    LPIenv[["nData"]]          ,
                    LPIenv[["nData"]]          ,
                    lag
                    )
             ) 
    }
    
    # Simple lagged products of decimated data,
    # works with strong codes.
    return( .Call( "lagged_products"          ,
                  LPIenv[["TX1"]][["cdata"]] ,
                  LPIenv[["TX2"]][["cdata"]] ,
                  LPIenv[["TX1"]][["idata"]] ,
                  LPIenv[["TX2"]][["idata"]] ,
                  LPIenv[["camb"]]           ,
                  LPIenv[["iamb"]]           ,
                  LPIenv[["nData"]]          ,
                  LPIenv[["nData"]]          ,
                  lag
                  )
           )
    
  }
