## file:clutterSuppress.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Ground clutter suppression as follows:
##
## 1. Scattered signal in ranges between rmin and rmax is 
##    solved by means of voltage-level inversion.
## 2. The solved profile is convolved with the transmission
##    envelope and the convolution is subtracted from the
##    receiver samples.
##
## Arguments:
##  txdata  A transmitter data list that contains named
##          vectors 'cdata' and 'idata'
##  rxdata  A receiver data list that cntains named
##          vectors 'cdata' and 'idata'
##  rmin    Smallest range from which clutter should
##          be suppressed
##  rmax    Largest range from which clutter should
##          be suppressed
##  ndata   Number of points in data vectors
##  clutterFraction Fraction of the full integration
##          period used for the clutter profile estimation
##          A float from the interval (0,1]
##
## Returns:
##  solution The solved clutter profile
##
## Clutter-suppressed receiver data is written to the
## vector rxdata[["cdata"]]
## 

clutterSuppress <- function( txdata , rxdata , rmin , rmax , ndata , clutterFraction )
  {

    # If rmin > rmax there will be nothing to subtract
    if( rmin > rmax ) return()

    # No reason to continue if ndata is not positive
    if( ndata <= 0 ) return()

    # Set negative ranges to zero
    rmin <- max( rmin , 0 )
    rmax <- max( rmax , 0 )
    

    # Number of range gates to solve
    nr <- rmax - rmin + 1

    # Initialize a fishs object
    e <- fishs.init( ncols = nr )

    # number of points used in clutter profile estimation
    nclutter <- round( ndata * min( clutterFraction , 1 ) )
    
    # Set correct storage modes
    storage.mode( ndata ) <- "integer"
    storage.mode( nclutter ) <- "integer"
    storage.mode( rmin ) <- "integer"
    storage.mode( rmax ) <- "integer"
    
    # Add data to the inverse problem
    nrow <- .Call( "clutter_meas",
                  txdata[["cdata"]],
                  txdata[["idata"]],
                  rxdata[["cdata"]],
                  rxdata[["idata"]],
                  ndata,
                  rmin,
                  rmax,
                  e[["Qvec"]],
                  e[["y"]]
                  )

    # Do not subtract if the number of measurement rows
    # is smaller than number of unknowns
    if( nrow < nr ){
      warning("Not enough data points for clutter suppression.")
      invisible( NULL )
    }

    # Otherwise solve the inverse problem
    fishs.solve(e)

    # The unmeasured points should be zero instead of NA
    e[["solution"]][is.na(e[["solution"]])] <- 0+0i

    # Do the actual subtraction
    ncor <- .Call( "clutter_subtract",
                  txdata[["cdata"]],
                  txdata[["idata"]],
                  rxdata[["cdata"]],
                  rxdata[["idata"]],
                  ndata,
                  rmin,
                  rmax,
                  e[["solution"]]
                  )

    invisible(e$solution)

  }
