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

decor.init <- function( ncols , ... )
{
    
    ## A new environment for the solver
    s <- new.env()
    
    ## Number of columns in theory matrix
    assign( 'ncol' , ncols , s )
    
    ## Diagonal of the Fisher information matrix, only real part needed in this decoder
    assign( 'QvecR' , rep(0,ncols) , s )
    
    ## Scaled measurements
    assign( 'yR'    , rep(0,ncols) , s )
    assign( 'yI'    , rep(0,ncols) , s )
    
    assign( 'FLOPS' , 0 , s )
    ## Make sure that the storage modes are
    ## correct for later c function calls
    storage.mode(s$QvecR) <- storage.mode(s$yR) <- storage.mode(s$yI) <- "double"
    storage.mode(s$ncol) <- "integer"
    storage.mode(s$FLOPS) <- 'double'
    
    ## return the environment
    return(s)
    
}
