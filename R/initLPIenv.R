## file:initLPIenv.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Allocate and initialise necessary vectors and variables
## for the actual lag profile inversion. This function is
## called once per integration period in each computing slave
##
## Arguments:
##  LPIenv.name Name of the LPI environment used for
##              the analysis.
##
## Returns:
##    Nothing, the udpated environment is stored on
##    the global workspace.
## 

initLPIenv <- function( LPIenv.name )
  {

    # Get the LPI environment (transferred as a list,
    # convert into an environment first)
    LPIenv <- as.environment( eval( LPIenv.name ) )

    # Allocate vector for the range ambiguity function
    assign( 'camb' , vector(mode='complex',length=(LPIenv[["nData"]]*LPIenv[["nDecimTX"]]))     , LPIenv )
    
    # Range ambiguity indices
    assign( 'iamb' , vector(mode='logical',length=(LPIenv[["nData"]]*LPIenv[["nDecimTX"]]))     , LPIenv )
    
    # Laged products
    assign( 'cprod', vector(mode='complex',length=LPIenv[["nData"]])                            , LPIenv )
    
    # Lagged product indices
    assign( 'iprod', vector(mode='logical',length=LPIenv[["nData"]])                            , LPIenv )
    
    # Lagged product variances
    assign( 'var'  , vector(mode='numeric',length=LPIenv[["nData"]])                            , LPIenv )
    
    # Theory matrix rows, one extra row because
    # theory_rows needs a temp vector
    assign( 'arows', vector(mode='complex',length=((max(LPIenv[["nGates"]])+1)*(LPIenv[["nBuf"]]+1))), LPIenv )
    
    # Indices for theory matrix rows, one extra row because
    # theory_rows needs a temp vector
    assign( 'irows', vector(mode='logical',length=((max(LPIenv[["nGates"]])+1)*(LPIenv[["nBuf"]]+1))), LPIenv )
    
    # Measurement vector
    assign( 'meas' , vector(mode='complex',length=LPIenv[["nBuf"]])                             , LPIenv )
    
    # Measurement variances
    assign( 'mvar' , vector(mode='numeric',length=LPIenv[["nBuf"]])                             , LPIenv )
    
    # Buffer row counter
    assign( 'nrows', as.integer(0)                                                              , LPIenv )

    # Copy the modified environment back
    # to the user workspace
    assign( paste(LPIenv.name) , LPIenv , envir=.GlobalEnv)

    return()

  }
