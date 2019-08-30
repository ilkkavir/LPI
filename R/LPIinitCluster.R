## file:LPIinitCluster.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.

## 
## Initialise the analysis cluster, which consists of:
##   - a master process (which calls this function)
##   - length(nodes) control processes running on the same
##     computer with the master process (if the input list
##     does not contain the entry localControl=FALSE)
##   - length(nodes) control processes running on any
##     computer on the cluster (names of these computers
##     are given in names(nodes) )
##   - length(nodes[[i]]) computing slaves running on each
##     remote computer nodes[i] (if nodes[[i]] is a
##     character vector, these nodes may be also on several
##     different computers)
##
##   The computing slaves do most of the  actual work,
##   the other processs are mainly for data transfer
##   (the control processes closest to master do all
##   disk I/O)
##
## Each remote control process will be given one
## integration period of raw data, whose ACF will be
## calculated in parallel by the computing slaves of the
## process. Thus, there are length(nodes) integration
## periods in parallel, and each of them has
## length(nodes[[i]]) parallel lag profiles inversions
## running.
##
## Arguments:
##  nodes  A list of remote control machine names and
##         definitions of computing slaves for each of them
##         e.g. nodes = list( tesla1=8, tesla2=8,
##                            tesla3=8, tesla4=8, tesla5=8)
##              or
##              nodes = list(
##            tesla1=c( rep('tesla1',8), rep('tesla3',8) ),
##            tesla4=c( rep('tesla4',8), rep('tesla5',8) )
##                          )
##
##         The former example starts a control process
##         on each computer of the tesla cluster, and
##         allocates one computing slave per core
##         (each of the computers has 8 cores).
##         Thus, five integration periods are analysed
##         in parallel with 8 lag profiles in parallel
##         in each of them. The latter one runs only two
##         integration periods at a time, but each of
##         them has 16 lag profiles in parallel.
##
##         Notice that the latter option leads to
##         signicantly larger amount of network traffic,
##         as the remote control nodes transfer the data
##         to each computing slave separately
## 
##         Alternatively, one can give an integer number,
##         which will start the given number of parallel
##         processes, running one integration period each
##         on localhost. Any combination of the above
##         inputs are also accepted.
##
##         The list nodes is treated as follows:
##           1. Put NAs to values <= 0
##           2. If only NA's were left from 1.,, do not
##              start a cluster
##              (analysis sequentially in the main process)
##           3. Unnamed entries are replaced with equal
##              number of entries localhost=1
##
##         More examples:
##
##           Start 5 parallel integration periods on
##           localhost, and another 5 on "remotecomputer":
##             nodes=list(5,remotecomputer=1,
##                      remotecomputer=1,remotecomputer=1,
##                      remotecomputer=1,remotecomputer=1)
##
##           Start 1 integration period with five parallel
##           lags in localhost and another similar one on 
##           "remotecomputer":
##             nodes=list(localhost=5,remotecomputer=5)
##
##           Start 4 parallel integration periods on both
##           remotecomputer1 and remotecomputer2, but do
##           not create the local control processes. This
##           requires that both computers have the input
##           and output data directories mounted on same
##           paths.
##             nodes=list(remotecomputer1=4,
##                        remotecomputer2=4,
##                        localControl=FALSE)
##
##           Do not use parallelism, solve everything
##           sequentally in the main process
##             nodes=NA
##
##             
##
##
## Returns:
##   ctrlcl  A list of class cluster of
##           the local control nodes
##
##           The corresponding lists of remote control
##           clusters and computing slaves
##           are stored on the cluster nodes
##

LPIinitCluster <- function( nodes , useXDR=FALSE )
  {

    # Check if nodes has an entry "localControl",
    # if not, use default (TRUE)
    localControl <- TRUE
    if(is.list(nodes)){
      if(is.logical(nodes[["localControl"]])) localControl <- nodes[["localControl"]]
    }
    
    # Replace negative values with NAs
    for(k in seq(length(nodes))){
      if(is.numeric(nodes[[k]])){
        if(nodes[[k]]<=0)  nodes[[k]] <- NA
      }
    }

    # If only NA values, we will run locally
    if(all(is.na(nodes))) return(NA)

    # Strip off all NAs (original NAs
    # and those from non-positive values)
    nodes[is.na(nodes)] <- NULL

    # Named nodes are left as such, unnamed are
    # assumed to denote the number of local
    # parallel integration periods
    nnames <- names(nodes)
    if(is.null(nnames)){
      ncnames <- rep(0,length(nodes))
    }else{
      ncnames <- nchar(names(nodes))
    }
    nodes2 <- c( nodes[ncnames>0] , rep( list( localhost=1 ) , sum( unlist( nodes[ ncnames==0 ] ) ) ) )

    # Remove the localControl entry
    nodes2[["localControl"]] <- NULL
    
    # Create the (optional) local control nodes
    if(localControl){
      # Create the local control nodes
      ctrlcl <- makeCluster( length( nodes2 ) , useXDR=useXDR )

      # Load packaages LPI and parallel to each of the local control nodes
      clusterEvalQ( ctrlcl , library(LPI) )

      # Create the remote computer control processes
      for(k in 1:length(nodes2)){
        # Run initialisation at the local control process to
        # create the remote control process and its slaves
        clusterCall( ctrlcl[k] , LPIinitRemoteNode , nodes2[k] , useXDR )
      }

    # Otherwise the remote nodes will
    # act as control nodes as well
    }else{

      # Create the remote control nodes directly
      ctrlcl <- makeCluster( names(nodes2) , useXDR=useXDR )
      
      # Load packaages LPI and parallel to
      # each of the remote control nodes
      clusterEvalQ( ctrlcl , library(LPI) )

      # Set remcl=NA on each node to notify that
      # the additional control step does not exist
      remcl <<- NA
      clusterExport( ctrlcl , 'remcl')

      # Initialise the computing slaves
      for(k in 1:length(nodes2)){
        clusterCall( ctrlcl[k] , LPI:::LPIinitComputingSlaves , nodes2[[k]] , useXDR )
      }
    }
    
    return(ctrlcl)

  }
