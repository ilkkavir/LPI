## file:dummy.add.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Dummy inverse problem solver that
## calculates simple averages.
## Data accumulation function.
##
## Arguments:
##  e       A dummy solver environemnt
##  M.data  Measurement vector
##  M.ambig Range ambiguity function
##  I.ambig Indices of non-zero ambiguity values
##  I.prod  Indices of usable lagged products
##  E.data  Measurement variance vector
##  nData   Number of points in data vectors
## 
## Returns:
##  success TRUE if the data was successfully added
##

dummy.add <- function( e , M.data , M.ambig , I.ambig , I.prod , E.data , nData )
  {

    # Call the C routine
    return( .Call( "dummy_add" ,
                  e[["msum"]] ,
                  e[["vsum"]] ,
                  e[["rmin"]] ,
                  e[["rmax"]] ,
                  M.data ,
                  M.ambig ,
                  I.ambig ,
                  I.prod ,
                  E.data ,
                  nData)
           )
    
  }
