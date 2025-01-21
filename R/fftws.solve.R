## file:fftws.solve.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## FFT deconvolution with optimized fft length and the fast fftw library.
## Final solver function.
##
## Arguments:
##  e      A ffts solver environment
##  rlims  Range gate limits
##
## Returns:
##  Nothing, the solution is assigned to the solver environment
##
fftws.solve <- function( e , rlims )
  {
    #
    # FFT deconvolution. Final solver function.
    #
    # I. Virtanen 2012
    #

    # Solve the lag profile by means of FFT
    sol <- fftw:IFFT( e[["fy"]]   / e[["sqfamb"]] , plan=e[["FFTplan"]] ) / e[["n"]]

    # Variance, the same value will be repeated at all ranges
    var <- e[["varsum"]] / as.double(e[["nmeas"]])  * mean( 1/ e[["sqfamb"]] )

    # Number of range gates
    nr  <- length(rlims) - 1

    # Final solution and variance vectors
    solution   <- rep(0+0i,nr)
    covariance <- rep(0,nr)
    
    for( r in seq(nr) ){

      # Lower limit of range gate
      r1            <- rlims[r]  + 1 - e[["rmin"]]

      # Upper limit of range gate
      r2            <- rlims[r+1] - e[["rmin"]]

      # All points have equal variances, calculate simple average
      solution[r]   <- mean( sol[r1:r2] , na.rm=TRUE )

      # Scale the variance
      covariance[r] <- var/(r2-r1+1)
    }

    # The background ACF cannot be measured with this technique, set it to NA.
    covariance <- c( covariance , NA )
    solution   <- c( solution , NA )

    # Assign the results to the solver environment.
    assign( 'solution'   , solution   , e )
    assign( 'covariance' , covariance , e )
    
    invisible()
    
  }
