## file:LPIrunLagprof.R
## (c) 2023- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## 
## Solve the MAP estimate of a lag profile,
## starting from raw voltage samples. This function
## first selects the correct LPIenv. 
## 
## Arguments:
##   LPIenvNames    Name of a list of lag profile inversion environments
##   lagind     An index, from which the lag number and integration period are decoded
## 
## Returns:
##  lagprof     A named list  containing the MAP estimate
##              of the lag profile together with
##              its (co)variance.
## 

LPIrunLagprof <- function( lagind , LPIenvNames , nlags )
  {

      # read the data lists from the workspace
      LPIenvs <- eval(LPIenvNames)

      # integration period within the data chunk
#      iper <- ceiling( lagind / LPIenvs[[1]][["nLags"]] )
      iper <- ceiling( lagind / nlags )

      # lag number
#      lag <- (lagind-1)%%LPIenvs[[1]][["nLags"]] + 1
      lag <- (lagind-1)%%nlags + 1

      # assign the correct part of the data to the workspace (because the solver reads it from there...)
      if (!is.null( LPIenv <<- LPIenvs[[iper]] )){
          
          # initialization
          initLPIenv(substitute(LPIenv))
          
          # the actual inversion
          lagprof <- LPIsolve( lag , substitute(LPIenv))
          
          # add the index to the profile to allow one to identify it later
          lagprof[["ii"]] <- lagind
          
          return(lagprof)
          
      }else{
          return(NULL)
      }
              
  } 
