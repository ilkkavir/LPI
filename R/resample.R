resample <- function( cdata ,  idata ,  ndata ,  nup , nfilter , nfirst , ipartial)
  {

    if( nup > nfilter ) stop("Upsampling is not currently supported, select nup <= nfilter.")
    return(
           .Call( "resample_R" , cdata , idata , ndata , nup , nfilter , nfirst ,ipartial )
           )

  }
