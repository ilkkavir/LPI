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

initLPIenvR <- function( LPIenv.name )
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
#    assign( 'arows', vector(mode='complex',length=((max(LPIenv[["nGates"]])+1)*(LPIenv[["nBuf"]]+1))), LPIenv )
    
    # Indices for theory matrix rows, one extra row because
    # theory_rows needs a temp vector
    assign( 'irows', vector(mode='logical',length=((max(LPIenv[["nGates"]])+1)*(LPIenv[["nBuf"]]+1))), LPIenv )
    
    # Measurement vector
#    assign( 'meas' , vector(mode='complex',length=LPIenv[["nBuf"]])                             , LPIenv )
    
    # Measurement variances
    assign( 'mvar' , vector(mode='numeric',length=LPIenv[["nBuf"]])                             , LPIenv )
    
    # Buffer row counter
    assign( 'nrows', as.integer(0)                                                              , LPIenv )

    ## this version uses separate double arrays for Re and Im
    assign( 'Rcomplex' , FALSE , LPIenv )
      
    ## make sure that the values are stored in correct format
    storage.mode( LPIenv$camb ) <- 'complex'
    storage.mode( LPIenv$iamb ) <- 'logical'
    storage.mode( LPIenv$cprod ) <- 'complex'
    storage.mode( LPIenv$iprod ) <- 'logical'
    storage.mode( LPIenv$var ) <- 'double'
#    storage.mode( LPIenv$arows ) <- 'complex'
    storage.mode( LPIenv$irows ) <- 'logical'
#    storage.mode( LPIenv$meas ) <- 'complex'
    storage.mode( LPIenv$mvar ) <- 'double'
    storage.mode( LPIenv$nrows ) <- 'integer'

      ## real and imaginary parts separately...
      assign( 'arowsR' , vector( mode='double' , length=((max(LPIenv[["nGates"]])+1)*(LPIenv[["nBuf"]]+1))), LPIenv )
      assign( 'arowsI' , vector( mode='double' , length=((max(LPIenv[["nGates"]])+1)*(LPIenv[["nBuf"]]+1))), LPIenv )
      assign( 'measR'  , vector( mode='double' , length=LPIenv[["nBuf"]]) , LPIenv )
      assign( 'measI'  , vector( mode='double' , length=LPIenv[["nBuf"]]) , LPIenv )

      storage.mode( LPIenv$arowsR ) <- 'double'
      storage.mode( LPIenv$arowsI ) <- 'double'
      storage.mode( LPIenv$measR ) <- 'double'
      storage.mode( LPIenv$measI ) <- 'double'

      
    # Copy the modified environment back
    # to the user workspace
    assign( paste(LPIenv.name) , LPIenv , envir=.GlobalEnv)

    return()

  }
