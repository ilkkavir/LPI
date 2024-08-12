## file:fishs.add.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Linear inverse problem solution by means of direct
## calculationof Fisher information matrix.
## Data accumulation function.
##
## Arguments:
##  e      A fishs solver environemnt
##  A.data Theory matrix rows as a vector (row-by-row)
##  I.data Indices of non-zero theory matrix elements
##  M.data Measurement vector
##  E.data Measurement variance vector
##
## Returns:
##  success TRUE if the rows were successfully added.
##

fishs.add <- function( e , A.data , I.data ,  M.data ,  E.data=1 )
  {

    # Number of theory rows to add
    nrow <- as.integer(length(M.data))

    # Variance vector
    E.data <- rep(E.data,length.out=nrow)

    # Check storage modes before calling the c function
    storage.mode(A.data) <- "complex"
    storage.mode(I.data) <- "logical"
    storage.mode(M.data) <- "complex"
    storage.mode(E.data) <- "double"
    storage.mode(nrow)   <- "integer"

    # Call the c function
    return( .Call( "fishs_add" , e[["Qvec"]] , e[["y"]] , A.data , I.data , M.data , E.data , e[["ncol"]] , nrow ))

  }
