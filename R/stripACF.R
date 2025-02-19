## file:stripACF.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Return an ACF list with only selected ranges and lags
##
## Arguments:
##   ACFlist   An ACF list returned by runLPI or stored by LPIsaveACF
##   rgates    Range gate indices
##   lags      Lag indices
##   fulCovar  TRUE if the ACFlist contains the full covariance matrices
##
## Returns:
##   ACFlist   A modified ACF list
##

stripACF <- function( ACFlist , rgates , lags , fullCovar=FALSE)
  {

    # An empty list for the output
    ACFlist2 <- list()

    # If rgates and lags are logical vectors
    # convert them into indices
    if(is.logical(rgates)) rgates <- which(rgates)
    if(is.logical(lags)) lags <- which(lags)

    # Pick the ACF and variance values
    ACFlist2[["ACF"]] <- ACFlist[["ACF"]][rgates,lags]
    ACFlist2[["var"]] <- ACFlist[["var"]][rgates,lags]
    
    # Make sure that ACF, var, and covariance are still arrays
    dim(ACFlist2[["ACF"]]) <- c( length(rgates) , length(lags) )
    dim(ACFlist2[["var"]]) <- c( length(rgates) , length(lags) )
    if(fullCovar){
      covdims <- dim(ACFlist[["covariance"]])
      ACFlist2[["covariance"]] <- ACFlist[["covariance"]][c(rgates,covdims[1]),c(rgates,covdims[2]),lags]
      dim(ACFlist2[["covariance"]]) <- c( (length(rgates)+1) , (length(rgates)+1) , length(lags) )
    }

    
    ACFlist2[["lag"]] <- ACFlist[["lag"]][lags]
    ACFlist2[["range"]] <- ACFlist[["range"]][rgates]
    ACFlist2[["nGates"]] <- pmin(rep(length(rgates),length(lags)),ACFlist[["nGates"]][lags])
    ACFlist2[["backgroundACF"]] <- ACFlist[["backgroundACF"]][lags]
    ACFlist2[["backgroundvar"]] <- ACFlist[["backgroundvar"]][lags]
    ACFlist2[["timeString"]] <- ACFlist[["timeString"]]
    ACFlist2[["time.s"]] <- ACFlist[["time.s"]]

    # Udpate names to match with the new indexing
    nlags <- length(lags)
    ngates <- length(rgates)

    names(ACFlist2[["range"]]) <- paste('gate',seq(ngates),sep='')
    names(ACFlist2[["lag"]]) <- paste('lag',seq(nlags),sep='')
    names(ACFlist2[["backgroundACF"]]) <- paste('lag',seq(nlags),sep='')
    names(ACFlist2[["backgroundvar"]]) <- paste('lag',seq(nlags),sep='')
    dimnames(ACFlist2[["ACF"]]) <- list(paste('gate',seq(ngates),sep=''),paste('lag',seq(nlags),sep=''))
    dimnames(ACFlist2[["var"]]) <- list(paste('gate',seq(ngates),sep=''),paste('lag',seq(nlags),sep=''))

      ACFlist2[["analysisTime"]] <- ACFlist[["analysisTime"]]
      ACFlist2[["FLOP"]] <- ACFlist[["FLOP"]]
#      ACFlist2[["addTime"]] <- ACFlist[["addTime"]]
      ACFlist2[["lagFLOP"]] <- ACFlist[["lagFLOP"]]
#      ACFlist2[["lagAddTime"]] <- ACFlist[["lagAddTime"]]

    return(ACFlist2)

  }
