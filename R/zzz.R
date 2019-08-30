## file:zzz.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Initialization when the package is loaded
##
## Arguments:
##   libname (see ?.onLoad)
##   pkgname (see ?.onLoad)
##
##

.onLoad <- function(libname,pkgname)
  {
    ctrlcl  <<- NA
    slavecl <<- NA
    remcl   <<- NA
    LPIdatalist.final <<- NULL
  }
