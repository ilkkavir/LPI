## file:fishs.solve.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Linear inverse problem solution by means of direct
## calculation of Fisher information matrix.
## Final solver function.
##
## Arguments:
##  e               A fishs solver environment
##  full.covariance Logical, full covariance matrix is calculated
##                  if TRUE, otherwise only variances are returned.
##
## Returns:
##  Nothing, the solution is assigned to
##  the solver environment
##

fishs.solve <- function( e , full.covariance = TRUE , ... )
  {

    # Allocate a matrix for the full
    # Fisher information matrix
    Q <- matrix( 0 , ncol=e[["ncol"]] , nrow=e[["ncol"]] )

    # Copy the upper triangular part form e$Qvec
    i <- 1
    for( k in seq( e$ncol ) ){
      Q[ k , k : e[["ncol"]] ] <- e[["Qvec"]][ i : ( i + ( e[["ncol"]] - k ) ) ]
      i <- i + e[["ncol"]] - k + 1
    }

    # The lower triangular part is
    # complex conjugate of the upper one
    Q         <- Q + Conj( t( Q ) )

    # The above row multiplies the diagonal
    # with 2, divide accordingly
    diag( Q)  <- diag( Q ) / 2

    # Select points at which the diagonal of Q is zero,
    # these points have not been measured at all and
    # need to be regularized before inverting the matrix
    nainds    <- Re( diag( Q ) ) == 0

    # Set unit values on the diagonal at unmeasured points.
    # This will not affect the other unknowns because
    # they cannot correlate with this one
    diag( Q )[ nainds ] <- 1


    # normalize diagonal of Q to 1 for better numerical stability
    dQsqrt <- sqrt(diag(Q))
    Qscale <- outer(dQsqrt,dQsqrt)
    Qnorm  <- Q/Qscale

      
    # Covariance matrix is inverse matrix of
    # the Fisher information matrix
    # Even if there were measurements the matrix might not be invertible
    # return NA matrix in this case
    covariance          <- tryCatch( solve( Qnorm ) , error=function(e){Qnorm*NA})

    # back to unnormaized units  
    covariance <- covariance/Qscale
      
    # Multiply the covariance matrix with e$y from right.
    # For some reason the direct matrix multiplication
    # with %*% does not work properly in some machines.
    solution <- rep(0+0i,e[["ncol"]])
    for(k in seq(e[["ncol"]])) solution[k] <- sum( covariance[k,] * e[["y"]] )

    # Set NAs to points that were not actually measured
    solution[ nainds ]  <- NA

    # Assign the solution to the solver environment e
    assign( 'solution'   , solution , e )

    # The full covariance matrix was already calculated, pick
    # the diagonal if that is enough.
    # Put NA to unmeasured points. 
    if( full.covariance ){
      covariance[ nainds ,        ] <- NA
      covariance[        , nainds ] <- NA
    }else{
      covariance                    <- diag( covariance )
      covariance[ nainds ]          <- NA
    }

    # Assign the covariance to the solver environment e
    assign( 'covariance' , covariance , e )
    
    invisible()
    
  }
