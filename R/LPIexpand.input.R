## file:LPIexpand.input.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.

##
## Expand input argument list or vector
## into the internally used format
##
## Arguments:
##  parvec A vector (or list)
##
## Returns:
##  outvec A named vector or list with elements
## "RX1", "RX2", "TX1", and "TX2".
##

LPIexpand.input <- function( parvec )
  {

    # Names of the input list / vector
    namevec <- names(parvec)
    
    # If the input does not have names attributes, assume
    # that the elements are in order RX1 , RX2 , TX1 , TX2
    # and repeat as necessary.
    if(is.null(namevec)){
      # Repeat the input
      outvec        <- rep(parvec,length.out=4)
      # Set names
      names(outvec) <- c( "RX1" , "RX2" , "TX1" , "TX2" )
      # Return the named vector / list
      return(outvec)
    }

    # If the input had names(s), start inspecting them

    # A vector for the output
    outvec <- rep(NA,4)
    names(outvec) <- c( "RX1" , "RX2" , "TX1" , "TX2" )

    # First look if any of the internally used
    # names is used in the input
    if( any(namevec=="RX1")) outvec[1] <- parvec["RX1"]
    if( any(namevec=="RX2")) outvec[2] <- parvec["RX2"]
    if( any(namevec=="TX1")) outvec[3] <- parvec["TX1"]
    if( any(namevec=="TX2")) outvec[4] <- parvec["TX2"]

    # If the vector had elements "RX1" , "RX2" , "TX1" ,
    # and "TX2", return them in correct order
    if( !any(is.na(outvec))) return(outvec)

    # If there are still missing values,
    # look for elements "RX" and "TX"
    if( is.na(outvec[1])){
      if(any(namevec=="RX")) outvec[1] <- parvec["RX"]
    }
    if( is.na(outvec[2])){
      if(any(namevec=="RX")) outvec[2] <- parvec["RX"]
    }
    if( is.na(outvec[3])){
      if(any(namevec=="TX")) outvec[3] <- parvec["TX"]
    }
    if( is.na(outvec[4])){
      if(any(namevec=="TX")) outvec[4] <- parvec["TX"]
    }

    # If the vector is now properly filled, return it
    if( !any(is.na(outvec))) return(outvec)

    # Now look for elements "TR1" and "TR2"
    if( is.na(outvec[1])){
      if(any(namevec=="TR1")) outvec[1] <- parvec["TR1"]
    }
    if( is.na(outvec[2])){
      if(any(namevec=="TR2")) outvec[2] <- parvec["TR2"]
    }
    if( is.na(outvec[3])){
      if(any(namevec=="TR1")) outvec[3] <- parvec["TR1"]
    }
    if( is.na(outvec[4])){
      if(any(namevec=="TR2")) outvec[4] <- parvec["TR2"]
    }

    # If the vector is now properly filled, return it
    if( !any(is.na(outvec))) return(outvec)

    # Finally remove the named elements from parvec and
    # try to fill the output vector
    parvec <- parvec[ nchar(namevec) == 0 ]
    if( length(parvec) > 0 ) outvec[is.na(outvec)] <- rep(parvec,length.out=sum(is.na(outvec)))

    # If the output is now full, return it
    if( !any(is.na(outvec))) return(outvec)

    # If still unsuccesfull, stop the whole analysis
    stop("Cannot parse the input vector ",paste(substitute(parvec)))
  }
