indexAdjust <- function( idata , ndata , shifts ){

  shifts2 <- as.integer(shifts)

  return( .Call( "index_adjust_R" , idata , ndata , shifts2 ) )

}
