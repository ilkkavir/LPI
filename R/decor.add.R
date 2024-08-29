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
##  I.data Indices of non-zero theory matrix elements
##  M.data Measurement vector
##  E.data Measurement variance vector
##
## Returns:
##  success TRUE if the rows were successfully added.
##


decor.add <- function( e , A.Rdata , A.Idata ,  I.data , M.Rdata , M.Idata ,  E.data , nrow )
  {

    # Call the c routine
    return( .Call( "decor_add" , e[["QvecR"]] , e[["yR"]] , e[["yI"]] , A.Rdata , A.Idata , I.data , M.Rdata , M.Idata , E.data , e[["ncol"]] , nrow , e[["FLOPS"]] ))

  }
