## file:LPIrun.remote.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

## 
## Send data to a remote computer and run analysis in there
## 
## Arguments:
##  LPIenv.name Name of the LPI environment to use, the
##              environment was copied on global workspace
##              by LPIsolve.ACF
##
## Returns:
##  ACFlist     A list that contains the solved ACF, its
##              covariance matrices, lags, etc. 
## 

LPIrun.remote <- function( LPIenv.name )
  {

    # Check if we are running in a cluster or not
    if( eval(LPIenv.name)[["iscluster"]]  ){

      # Check that this is a local control node
      if(!is.na(remcl)){
      
        # Send the data environment to the remote node
        clusterExport( remcl , paste(LPIenv.name)  )

        # Run the remote analysis
        ACF <-  clusterCall( remcl , LPI:::LPIrun , LPIenv.name )

        # Return the ACF
        return( ACF[[1]] )
      }
      
    }

    # If the analysis is run in a single process,
    # or if this is a remote control node,
    # just run the analysis in this process
    return( LPI:::LPIrun( LPIenv.name ) )

   }
