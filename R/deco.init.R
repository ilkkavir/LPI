## file:deco.init.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Matched filter decoder. Initialization function.
##
## Arguments:
##  ncols Number of unknowns (theory matrix columns)
##  ...   Additional arguments are allowed by not used
##        in order to make the solver more compatible
##        with others.
##
## Returns:
##  e    A deco solver environment
##

deco.init <- function( ncols , ... )
  {

    # A new environment for the solver
    s <- new.env()

    # Number of columns in theory matrix
    assign( 'ncol' , ncols , s )

    # Diagonal of the Fisher information matrix
    assign( 'Qvec' , rep(0,ncols) , s )

    # Scaled measurements
    assign( 'y'    , rep(0,ncols) , s )

    # Make sure that the storage modes are
    # correct for later c function calls
    storage.mode(s$Qvec) <- storage.mode(s$y) <- "complex"
    storage.mode(s$ncol) <- "integer"

    # return the environment
    return(s)

  }
