## file:deco.add.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Matched filter decoder. Data accumulation function.
##
## Arguments:
##  e      A deco solver environemnt
##  A.data Theory matrix rows as a vector (row-by-row)
##  M.data Measurement vector
##  E.data Measurement variance vector
##
## Returns:
##  success TRUE if the rows were successfully added.
##


deco.add <- function( e , A.data ,  M.data ,  E.data=1 )
  {
    # Number of theory rows
    nrow <- as.integer(length(M.data))

    # Measurement variance vector
    E.data <- rep(E.data,length.out=nrow)

    # Set storage modes
    storage.mode(A.data) <- "complex"
    storage.mode(M.data) <- "complex"
    storage.mode(E.data) <- "double"
    storage.mode(nrow)   <- "integer"

    # Call the c routine
    return( .Call( "deco_add" , e$Qvec , e$y , A.data , M.data , E.data , e$ncol , nrow ))

  }
