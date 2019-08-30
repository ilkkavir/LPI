## file:LPIinitRemoteNode.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.

##
## Establish a scoket connection between a local
## control process and a remote control process
## and initialise the computing slaves at the remote
##
## Arguments:
##   remnode  a list of the form
##            remname = list(compslave1,compslave2,...)
##            where remname is the remote analysis 
##            computer and the list constains either
##            the number of computing slaves
##            on that particular computer, or a vector
##            of computer names at which to create the
##            computing slaves
##

LPIinitRemoteNode <- function( remNode , useXDR )
  {

    # Node name
    nodeName <- names(remNode)

    # If nodeName is localhost, do not start the remote
    # control process but make direct connections
    # to the slaves instead.
    if( nodeName == "localhost"){
      remcl <<- NA
      LPI:::LPIinitComputingSlaves( remNode[[1]] , useXDR )
      return(remcl)
    }

    # Establish the connection to the remote control nodes
    remcl   <<- makeCluster( names(remNode) , useXDR=useXDR )

    # Load package LPI
    clusterEvalQ( remcl , library(LPI) )

    # Initialise the computing slaves
    clusterCall( remcl , LPIinitComputingSlaves , remNode[[1]] , useXDR )

    return(remcl)

  }
