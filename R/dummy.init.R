## file:dummy.init.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Dummy inverse problem solver that calculates
## simple averages.
## Initialization function.
##
## Arguments:
##  rrange extreme ranges to be solved c(rmin,rmax)
##
## Returns:
##  s     A dummy solver environment
##

dummy.init <- function( rrange )
  {
    
    # A new environment for the solver
    s <- new.env()

    # Number of ranges (this is different
    # from number of final range gates)
    nr <- abs(diff(rrange))

    # A vector for sum of weighted measurements
    msum <- rep(0+0i,nr)

    # A vector for sum of information
    vsum <- rep(0,nr)

    # Minimum range
    rmin <- min(rrange)

    # Maximum range
    rmax <- max(rrange)

    # Make sure that storage modes are correct
    storage.mode(msum) <- "complex"
    storage.mode(vsum) <- "double"
    storage.mode(rmin) <- "integer"
    storage.mode(rmax) <- "integer"

    # Assign the variables to the environment
    assign( 'msum' , msum , s )
    assign( 'vsum' , vsum , s )
    assign( 'rmin' , rmin , s )
    assign( 'rmax' , rmax , s )

    # Return the environment
    
    return(s)
    
  }
