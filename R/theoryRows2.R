## file:theoryRows.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Form theory matrix rows for lag profile inversion
##
## Arguments:
##   LPIenv    A LPI environment
##   lag       Lag number
##
## Returns:
##   success  TRUE if at least one theory matrix row
##            was successfully produces, FALSE otherwise.
##
## The rows are written to LPIenv[["arows"]],
## the correspoding measurements to LPIenv[["meas"]],
## variance to LPIen[["mvar"]], and number of rows
## generated to LPIenv[["nrows"]]
##           
##

theoryRows2 <- function( LPIenv , lag )
  {

    # Call the C routine
    return( .Call( "theory_rows2" ,
                  LPIenv[['camb']] ,
                  LPIenv[['iamb']] ,
                  LPIenv[['cprod']],
                  LPIenv[['iprod']],
                  LPIenv[['var']] ,
                  LPIenv[['nData']] ,
                  LPIenv[['nCur']] ,
                  as.integer(LPIenv[['nCur']]+LPIenv[['nBuf']]) ,
                  LPIenv[['rangeLimits']] ,
                  LPIenv[['nGates']][lag] ,
                  LPIenv[['arowsR']] ,
                  LPIenv[['arowsI']] ,
                  LPIenv[['irows']] ,
                  LPIenv[['measR']] ,
                  LPIenv[['measI']] ,
                  LPIenv[['mvar']],
                  LPIenv[['nrows']],
                  LPIenv[["backgroundEstimate"]],
                  LPIenv[["remoteRX"]]
                  )
           )
  }
