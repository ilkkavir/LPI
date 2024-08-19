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

fishsr.add <- function( e , A.Rdata , A.Idata , I.data ,  M.Rdata , M.Idata ,  E.data , nrow )
{


    # Call the c function
    return( .Call( "fishsr_add" , e[["QvecR"]] , e[["QvecI"]] , e[["yR"]] , e[["yI"]] , A.Rdata , A.Idata , I.data , M.Rdata , M.Idata , E.data , e[["ncol"]] , nrow ))

  }
