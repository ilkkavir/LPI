## file:rlips.solve2.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

## Call rlips.solve after regularization for
## unknowns that were not measured at all
## Set the corresponding values to NA before returning
##
## Arguments:
##  e              An rlips solver environment
##  fullCovariance Logical, if TRUE full covariance matrix
##                 is calculated, otherwise only the
##                 variances.
## Returns:
##  Nothing, the solution is assigned to the
##  solver environment.
##

rlips.solve2 <- function( e , full.covariance = TRUE )
  {
    # Read data from gpu memory
    rlips.get.data( e )

    # Select non-measured points
    nainds <- which( Re( diag( e$R.mat ) ) == 0 )

    # Add regularizing imaginary measurements
    regrow <- rep(0+0i,e$ncols)
    for( n in nainds ){
      regrow[]  <- 0+0i
      regrow[n] <- 1+0i
      rlips.add( e , A.data = regrow , M.data = 1.0+0.0i )
    }

    # Solve the problem
    rlips.solve( e , calculate.covariance = TRUE , full.covariance = full.covariance )

    # Set NAs to appropriate points in the solution
    sol <- e$solution
    sol[nainds] <- NA
    assign( 'solution' , sol , e )

    # Set the unmeasured points to NA
    # in the covariance matrix as well.
    covar <- e$covariance
    if( full.covariance ){
      covar[ , nainds ] <- NA
      covar[ nainds , ] <- NA
    }else{
      covar[nainds] <- NA
    }

    # Assign the covariance matrix to the solver environment
    assign( 'covariance' , covar , e )
    
    invisible()
    
  }
