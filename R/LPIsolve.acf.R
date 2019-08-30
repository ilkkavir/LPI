## file:LPIsolve.acf.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

## This function is run in local control nodes. 
## Read data for one integration period, send it for
## analysis in a remote computer, and write the returned
## ACF to file
## 
## Arguments: 
##  intPeriod  Integration period number, counted from
##             LPIparam[["firstTime"]] in steps of
##             LPIparam[["timeRes.s"]]
##
## Returns:
##  intPeriod  The integration period number.
## 

LPIsolve.acf <- function( intPeriod , LPIparam )
  {
    # Load packages that are needed for reading the data
    for( pn in LPIparam[["inputPackages"]] ){
      require( pn , character.only=TRUE )
    }

    # Parameter list update
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
    
            # Frequency mixing, filtering, etc., the output is
            # collected in a list and stored on the user workspace
            LPIdatalist.final <<- prepareLPIdata( LPIparam , LPIdatalist.raw )

            # Call the function that will send the data to
            # proper place and run the actual analysis
            ACF <- LPI:::LPIrun.remote( substitute(LPIdatalist.final) )

            # Store the results
            eval( as.name( LPIparam[["resultSaveFunction"]]) )( LPIparam , intPeriod , ACF )

          }
        }
    }
    
    # Return the integration period
    # number to the main process
    return(intPeriod)
    
  }
