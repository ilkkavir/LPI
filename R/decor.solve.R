## file:deco.solve.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Matched filter decoder. Final solver function.
##
## Arguments:
##  e               A deco solver environment
##  full.covariance Logical, full covariance matrix is
##                  calculated if TRUE, otherwise only
##                  variances are returned.
##
## Returns:
##  Nothing, the solution is assigned to
##  the solver environment
##

decor.solve <- function( e , ... )
{
    ## Diagonal of the Fisher information matrix
    ## (Matched filter decoding is equivalent with assuming
    ## that the nondiagonal elements are zeros)
    Qdiag <- e[["QvecR"]]
    
    ## The points at which Qdiag is zero were not measured
    ## at all, flag these points
    nainds <- Qdiag == 0
    
    ## Put unit values to the unmeasured points. This does
    ## not affect the other points as they cannot
    ## correlated with the unmeasured ones.
    Qdiag[nainds] <- 1
    
    ## Variance is simply the inverse of the diagonal
    ## of the Fisher information
    variance  <- 1 / Qdiag
          
    ## Assign the solution to the solver environment
    y <- e[["yR"]] + 1i*e[["yI"]]
    assign( 'solution'   , variance * y , e )
    
    ## Set NAs to the unmeasured points
    e[["solution"]][nainds] <- NA
    
    ## Same for the variances
    assign( 'covariance' , variance , e )
    e[["covariance"]][nainds] <- NA

    invisible()
    
  }
