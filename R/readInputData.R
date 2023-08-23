## file:LPIsolve.acf.R
## (c) 2023- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

## 
## Read raw data from one integration period, and apply filtering and decimation
## 
## 
## 
## Arguments: 
##  intPeriod  Integration period number, counted from
##             LPIparam[["firstTime"]] in steps of
##             LPIparam[["timeRes.s"]]
##
## Returns:
##  LPIenv     An "LPI environment" that contains the data vectors
## 

readInputData <- function( intPeriod , LPIparam )
  {
    # Load packages that are needed for reading the data
    for( pn in LPIparam[["inputPackages"]] ){
      require( pn , character.only=TRUE )
    }

    # Parameter list update (I do not thing this is needed, but will not cause any harm either...)
    LPIparam <- eval( as.name( LPIparam[["paramUpdateFunction"]] ))( LPIparam , intPeriod )

    if( !is.null(LPIparam)){
        # Read raw data, name of the data input function
        # should be stored in a character string
        LPIdatalist.raw   <- eval(as.name(LPIparam[["dataInputFunction"]]))( LPIparam , intPeriod )

        # If data reading was successfull
        if(LPIdatalist.raw[["success"]]){

          # require that there are at least some TX and RX samples
          if( (sum(LPIdatalist.raw[["RX1"]][["idata"]]) > 0) &
              (sum(LPIdatalist.raw[["RX2"]][["idata"]]) > 0) &
              (sum(LPIdatalist.raw[["TX1"]][["idata"]]) > 0) &
              (sum(LPIdatalist.raw[["TX2"]][["idata"]]) > 0)){
    
            # Frequency mixing, filtering, etc.
            LPIdatalist.final <- prepareLPIdata( LPIparam , LPIdatalist.raw )


          }
        }
    }


    # 
    # Return the data enviroment
    return(LPIdatalist.final)
    
  }
