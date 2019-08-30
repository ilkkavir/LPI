## file:lagprodVar.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Variances of lagged products. Calculated
## as lagged products of average power values. 
##
## Arguments:
##   LPIenv A LPI environment
##   lag    Lag number
##
## Returns:
##   success  TRUE if a variance estimate was successfully
##            calculated for at least one data point,
##            FALSE otherwise.
## The variances are (over)written to LPIenv[["var"]]
##

lagprodVar <- function( LPIenv , lag )
  {

    # Make sure that lag is an integer
    storage.mode(lag) <- "integer"
  
    # Call the C function
    return( .Call( "lagged_products_r"        ,
                  LPIenv[["RX1"]][["power"]] ,
                  LPIenv[["RX2"]][["power"]] ,
                  LPIenv[["var"]]            ,
                  LPIenv[["nData"]]          ,
                  LPIenv[["nData"]]          ,
                  lag
                  )
           ) 
  }
