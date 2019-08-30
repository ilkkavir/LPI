## file:fishs.init.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Linear inverse problem solution by means of direct calculation
## of Fisher information matrix. Initialization function.
##
## Arguments:
##  ncols Number of unknowns (theory matrix columns)
##
## Returns:
##  s     A fishs solver environment
##

fishs.init <- function( ncols , ... )
  {
    # New environment for the solver
    s <- new.env()

    # Number of columns in the theory matrix
    assign( 'ncol' , ncols , s )
    
    # A vector for upper triangular part of
    # the Fisher information matrix
    assign( 'Qvec' , rep(0+0i,(ncols*(ncols+1)/2)) , s )
    
    # A vector for weighted measurements
    assign( 'y'    , rep(0+0i,ncols) , s )

    # Make sure that the storage modes are
    # correct for later c function calls
    storage.mode(s$Qvec) <- storage.mode(s$y) <- "complex"
    storage.mode(s$ncol) <- "integer"

    return(s)

  }
