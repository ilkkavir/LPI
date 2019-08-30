## file:dummy.solve.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Dummy inverse problem solver that
## calculates simple averages.
## Final solver function.
##
## Arguments:
##  e      A dummy solver environment
##  rlims  Range gate limits
##
## Returns:
##  Nothing, the solution is assigned to
##  the solver environment
##

dummy.solve <- function( e , rlims )
  {
    #
    # 
    # Final solver function.
    # 
    # I. Virtanen 2012
    #

    # Number of range gates
    nr <- length(rlims) - 1

    # Vectors for the solution and variance
    solution <- rep(0+0i,nr)
    covariance <- rep(0,nr)

    # Range integration for the data points that have
    # the best possible resolution at this point.
    for( r in seq(nr) ){

      # Lower limit of this range gates
      r1 <- rlims[r] - rlims[1] + 1

      # Upper limit of this range gate
      r2 <- rlims[r+1] - rlims[1]
      
      # The vector e$msum contains variance weighted sum,
      # we can simply sum its elements.
      solution[r] <- sum(e[["msum"]][r1:r2])

      # The vector e$vsum contains informations, sum them.
      covariance[r] <- sum(e[["vsum"]][r1:r2])

    }

    # Variance is inverse of the information
    covariance <- c( 1/covariance , NA )

    # Multiply the solution with the final variances
    solution <- c( solution , NA ) * covariance

    # Vectors solution and covariance will now contain
    # variance-weighted averages of the lag profiles
    # and their variances. Assign to the solver environment
    assign( 'solution'   , solution   , e )
    assign( 'covariance' , covariance , e )

    invisible()
  
  }
