## file:laggedProducts.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Calculation of lagged products
##
## Arguments:
##   LPIenv An LPI environment
##   lag    Lag number
##
##
## Returns:
##   success  TRUE if at least one lagged product was
##            successfully calculated, FALSE otherwise.
##
## The lagged products are (over)written to
## the vector LPIenv[["cprod."]]
##

laggedProducts <- function( LPIenv , lag )
  {

    # Make sure that the lag number is an integer
    storage.mode(lag) <- "integer"
  
    # Call the c function
    return( .Call( "lagged_products" ,
                  LPIenv[["RX1"]][["cdata"]] ,
                  LPIenv[["RX2"]][["cdata"]] ,
                  LPIenv[["RX1"]][["idata"]] ,
                  LPIenv[["RX2"]][["idata"]] ,
                  LPIenv[["cprod"]]          ,
                  LPIenv[["iprod"]]          ,
                  LPIenv[["nData"]]          ,
                  LPIenv[["nData"]]          ,
                  lag
                  )
           )    
  }
