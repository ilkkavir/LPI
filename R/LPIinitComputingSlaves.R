## file:LPIinitComputingSlaves.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.

##
## Init the actual worker processes,
## "computing slaves", of a remote control node.
##
## Arguments:
##  slaveNodes Cluster node definition, either an integer
##             number of cluster nodes or a string vector
##             of host names.
##
## Returns:
##  slavecl   An object of class SOCKcluster. The same
##            object is also stored on the global workspace
##

LPIinitComputingSlaves <- function( slaveNodes , useXDR )
  {

    # If only one slave, do not allocate it but 
    # run analysis in the control process
    if( slaveNodes == 1 ){
      return( slavecl <<- NA )
    }
    
    # Create the cluster of computing slaves
    slavecl <<- makeCluster( slaveNodes , useXDR=useXDR )

    # Load LPI package to each of the nodes
    clusterEvalQ( slavecl , library(LPI) )
    
    return(slavecl)
    
  }
