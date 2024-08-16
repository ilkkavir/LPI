## file:LPIsolve.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

## 
## Solve the MAP estimate of a lag profile,
## starting from raw voltage samples
## 
## Arguments:
##   LPIenv     A lag profile inversion environment
##   lag        Lag number, all fractional lags from
##              LPIenv[["lagLimits"]][lag] to 
##              LPIenv[["lagLimits"]][lag+1]-1
##              are integrated in the same profile
## 
## Returns:
##  lagprof     A named list  containing the MAP estimate
##              of the lag profile together with
##              its (co)variance.
## 

LPIsolve <- function( lag , LPIenv.name , intPeriod=0)
{
    
#    if(lag==1){
#        EnvFile <-  paste('Env_',intPeriod,'.Rdata',sep='')
#        RprofFile <- paste('Rprof_',intPeriod,'.out',sep='')
#        FlopsFile <- paste('FLOPS_',intPeriod,'.Rdata',sep='')
#        Rprof(filename=RprofFile,memory.profiling=TRUE,gc.profiling=TRUE,line.profiling=TRUE)
#    }
    
    
    ## Get the LPI environment from the global workspace
    LPIenv <- eval(LPIenv.name)
    
    ## Return immediately if number of gates is <= 0
    if( LPIenv[["nGates"]][lag] <= 0 ) return(list(lagnum=lag))
    
    ## If rlisp is used, make sure it has been loaded.
    ## rlips is not required in startup in order to
    ## allow analysis without installing it. Other
    ## solvers are included in the LPI package.
    ## Switch quietly to fishs if rlips is not available.
    if(LPIenv$solver=="rlips"){
        require(rlips) -> rres
        if( !rres ) assign( 'solver' , 'fishs' , LPIenv )
    }
    
    ## Initialise the inverse problem solver
    if(LPIenv$solver=="rlips"){
        solver.env <- rlips.init( ncols = LPIenv$nGates[lag] + 1 , nrhs = 1 , type = LPIenv$rlips.options[["type"]] , nbuf = LPIenv$rlips.options[["nbuf"]] , workgroup.size = LPIenv$rlips.options[["workgroup.size"]] )
    }else if ( LPIenv$solver=="fishs" ){
        solver.env <- fishs.init2( LPIenv[["nGates"]][lag] + 1 )
    }else if ( LPIenv[["solver"]]=="deco" ){
        solver.env <- deco.init( LPIenv[["nGates"]][lag] + 1 )
    }else if ( LPIenv[["solver"]]=="dummy" ){
        solver.env <- dummy.init( range( LPIenv[["rangeLimits"]][ 1 : (LPIenv[["nGates"]][lag]+1) ]) )
    }else if ( LPIenv[["solver"]]=="ffts" ){
        solver.env <- ffts.init( LPIenv[["nGates"]][lag] , LPIenv[["TX1"]][["idata"]][1:LPIenv[["nData"]]])
    }
    
    ## Copy of LPIenv[["nData"]]
    ndcpy <- LPIenv[["nData"]]
    
    ## Walk through all fractional time-lags
    for( l in seq( LPIenv[["lagLimits"]][lag] , ( LPIenv[["lagLimits"]][lag+1] - 1 ) )){
        
        ## If the lag is longer than the data vector
        ## it cannot be calculated
        if( l >= LPIenv[["nData"]]) break
        
        ## Current position in data vector, we will skip the first nGates samples
        assign( "nCur" , as.integer(LPIenv[["rangeLimits"]][LPIenv[["nGates"]][lag]+1]+1) , LPIenv)

        ## Calculate the lagged products
        laggedProducts( LPIenv , l )

        ## Variances of lagged products
        lagprodVar( LPIenv , l )

        ## Calculate range ambiguity function
        rangeAmbiguity( LPIenv , l )
        
        ## Optional pre-averaging of lag-profiles
        if( !is.null( LPIenv[["nCode"]] )){
            if( !is.na( LPIenv[["nCode"]] )){
                if( LPIenv[["nCode"]] > 0 ){
                    averageProfiles( LPIenv , l )
                    nd <- min( LPIenv[["nData"]] , which( diff( LPIenv[["TX1"]][["idata"]] ) == 1 )[ LPIenv[["nCode"]] + 1 ] )
                    LPIenv[["nData"]] <- ifelse( is.na(nd) , LPIenv[["nData"]] , nd )
                    ## Approximate the variance.
                    ## This is not exactly accurate!
                    if(!is.na(nd))  LPIenv[["var"]] <- LPIenv[["var"]]  / ( sum(diff(LPIenv[["TX1"]][["idata"]])==1) /  LPIenv[["nCode"]] )
                }
            }
        }
        
        ## Solvers "dummy" and "ffts" operate
        ## directly with the product vectors
        if( LPIenv[["solver"]]=="dummy" ){
            
            dummy.add( e      = solver.env        ,
                      M.data  = LPIenv[["cprod"]] ,
                      M.ambig = LPIenv[["camb"]]  ,
                      I.ambig = LPIenv[["iamb"]]  ,
                      I.prod  = LPIenv[["iprod"]] ,
                      E.data  = LPIenv[["var"]] , nData = as.integer( LPIenv[["nData"]] - l ) )
            
        }else if( LPIenv[["solver"]]=="ffts"){
            
            ffts.add( e       = solver.env        ,
                     M.data  = LPIenv[["cprod"]] ,
                     M.ambig = LPIenv[["camb"]]  ,
                     I.ambig = LPIenv[["iamb"]]  ,
                     I.prod  = LPIenv[["iprod"]] ,
                     E.data  = LPIenv[["var"]]   ,
                     nData   = as.integer(LPIenv[["nData"]] - l)
                     )
            
            ## Other solvers need theory matrix rows
        }else{
            
            ## Produce theory matrix rows in
            ## (small) sets and add them to the solver
            FLOPS <- 0
            NROWS <- 0
            addtime <- system.time({
            while( newrows <- theoryRows2( LPIenv , lag ) ){
                NROWS <- NROWS + LPIenv[["nrows"]]
                ## If new rows were produced
                if( LPIenv[["nrows"]]>0){
                    
                    ## select the correct solver
                    if(LPIenv$solver=="rlips"){
                        
                        rlips.add( e = solver.env ,
                                  A.data = LPIenv[["arows"]][1:(LPIenv[["nrows"]]*(LPIenv[["nGates"]][lag]+1))] ,
                                  M.data = LPIenv[["meas"]][1:LPIenv[["nrows"]]] ,
                                  E.data = LPIenv[["mvar"]][1:LPIenv[["nrows"]]]
                                  )
                        
                    }else if(LPIenv$solver=='fishs'){
                        
                        ## FLOPS <- FLOPS + as.double(fishs.add( e = solver.env ,
                        ##           A.data = LPIenv[["arows"]] ,
                        ##           I.data = LPIenv[["irows"]] ,
                        ##           M.data = LPIenv[["meas"]] ,
                        ##           E.data = LPIenv[["mvar"]] ,
                        ##           nrow = LPIenv[['nrows']]
                        ##           ))
                        FLOPS <- FLOPS + as.double(fishs.add2( e = solver.env ,
                                  A.Rdata = LPIenv[["arowsR"]],
                                  A.Idata = LPIenv[["arowsI"]] ,
                                  I.data = LPIenv[["irows"]] ,
                                  M.Rdata = LPIenv[["measR"]] ,
                                  M.Idata = LPIenv[["measI"]] ,
                                  E.data = LPIenv[["mvar"]],
                                  nrow = LPIenv[["nrows"]]
                                  ))
                        
                    }else if(LPIenv[["solver"]] == "deco" ){
                        
                        deco.add( e = solver.env ,
                                 A.data = LPIenv[["arows"]][1:(LPIenv[["nrows"]]*(LPIenv[["nGates"]][lag]+1))] ,
                                 I.data = LPIenv[["irows"]][1:(LPIenv[["nrows"]]*(LPIenv[["nGates"]][lag]+1))] ,
                                 M.data = LPIenv[["meas"]][1:LPIenv[["nrows"]]] ,
                                 E.data = LPIenv[["mvar"]][1:LPIenv[["nrows"]]]
                                 )
                        
                    }
                }
            }
            })
        }
        
        ## Make sure that the original value is
        ## stored in LPIenv[["nData"]]
        LPIenv[["nData"]] <- as.integer(ndcpy)
        
        
    }
    
    ## Solve the inverse problem
    if(LPIenv$solver=="rlips"){
        rlips.solve2( e = solver.env ,full.covariance = LPIenv[["fullCovar"]])
    }else if(LPIenv$solver=="fishs"){
        fishs.solve2( e = solver.env , full.covariance = LPIenv[["fullCovar"]] )
    }else if(LPIenv[["solver"]]=="deco"){
        deco.solve( e = solver.env )
    }else if(LPIenv[["solver"]]=="dummy"){
        dummy.solve( e = solver.env , LPIenv[["rangeLimits"]][1:(LPIenv[["nGates"]][lag]+1)])
    }else if( LPIenv[["solver"]]=="ffts"){
        ffts.solve( e = solver.env , LPIenv[["rangeLimits"]][1:(LPIenv[["nGates"]][lag]+1)])
    }
    
    ## Create the return environment
    lagprof <- new.env()

    ## Assign the solution to the new environment
    assign( "lagprof" , solver.env[["solution"]] , lagprof )
    assign( "covariance" , solver.env[["covariance"]] , lagprof )
    assign( "lagnum" , lag , lagprof )
    assign( "FLOPS" , FLOPS , lagprof )
    assign( "addtime" , addtime , lagprof)
    assign( "NROWS" , NROWS , lagprof )
    
    ## Kill the solver object
    if(LPIenv$solver=="rlips") rlips.dispose(solver.env)
    
#    if(lag==1){
#        Rprof(NULL)
#        save(FLOPS=FLOPS,NROWS=NROWS,addtime=addtime,file=FlopsFile)
#        save(LPIenv=LPIenv,file=EnvFile)
#    }
    
    ## Conversion to list because it is faster to transfer
    return(as.list(lagprof))
    
}
