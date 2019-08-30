// file:index_adjust.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/*
  Adjust tx / rx indices. The rising edges are shifted 
  shifts[0] samples and the falling edges shifts[1] 
  samples towards larger indices. Also negative
  shifts are allowed.

  This function allocates new data vectors

  Arguments:
   idata   ndata integer vector of TX pulse / RX positions
   ndata   Number of data points in idata
   shifts  2-vector of shifts 
           (shifts at rising and falling edges)
  
  Returns:
   ans       A list with elements
              idata   Index vector after adjustements
	      success Logical, set if all processing
	              was successful
*/

SEXP index_adjust_R( SEXP idata , SEXP ndata , SEXP shifts )
{
  SEXP ans;
  SEXP idata_new;
  SEXP s;
  SEXP names;
  char *cnames[2] = {"idata","success"};
  int *inew;
  int *iold;
  register uint64_t k;


  // Output list ans[[1]] = idata , ans[[2]] = success
  PROTECT( ans = allocVector( VECSXP , 2 ) );

  // Allocate the new logical vector
  PROTECT( idata_new = allocVector( LGLSXP , *(INTEGER(ndata)) ) );

  // A pointer to the new data vector
  inew = LOGICAL( idata_new );

  // A pointer to the old data vector
  iold = LOGICAL( idata );

  // Copy data from old to new
  for( k = 0 ; k < *(INTEGER(ndata)) ; ++k ){
    inew[k] = iold[k];
  }

  // The success logical
  PROTECT( s = allocVector( LGLSXP , 1 ) );

  // The actual work
  s = index_adjust( idata_new , ndata , shifts );

  // Collect the data into the return list
  SET_VECTOR_ELT( ans , 0 , idata_new );
  SET_VECTOR_ELT( ans , 1 , s );

  // Set the name attributes
  PROTECT( names = allocVector( STRSXP , 2 ));
  SET_STRING_ELT( names , 0 , mkChar( cnames[0] ) );
  SET_STRING_ELT( names , 1 , mkChar( cnames[1] ) );
  setAttrib( ans , R_NamesSymbol , names);

  UNPROTECT(4);

  return(ans);

}

/*
  Adjust TX / RX indices. The rising edges are shifted
  shifts[0] samples and the falling edges shifts[1] 
  samples towards larger indices.
  Also negative shifts are allowed.

  This function overwrites the idata vector

  Arguments:
   idata   ndata integer vector of TX pulse / RX positions
   ndata   Number of data points in idata
   shifts  2-vector of shifts
           (shifts at rising and falling edges)
  
  Returns:
   success 1 if all processing was successful, 0 otherwise

*/

SEXP index_adjust( SEXP idata , SEXP ndata , SEXP shifts)
{
  int *id = INTEGER(idata);
  int *nd = INTEGER(ndata);
  int *sh = INTEGER(shifts);
  // temporary variables
  int sh1;
  register int64_t k;
  int lasttrue;
  int ncut;
  int nadd;
  // for the return value
  SEXP success;
  int *isuccess;

  // Allocate the return value and initialise it
  PROTECT(success = allocVector(LGLSXP,1));
  isuccess = LOGICAL(success);
  *isuccess = 1;

  // The shift on rising edges is done by 
  //shifting the whole index vector

  // Find the last true index in the whole vector,
  // it will be needed later
  lasttrue = 0;
  for( k = ( *nd - 1 ) ; k >= 0 ; --k ){
    if( id[k] ){
      lasttrue = k;
      break;
    }
  }

  // If sh[0] < 0, shift towards smaller indices
  if( sh[0] < 0 ){
    for( k = 0 ; k < ( *nd + sh[0] ) ; ++k ){
      id[k] = id[ k - sh[0] ];
    }
    // The last value is repeated in the remaining points
    for( k = ( *nd + sh[0]) ; k < *nd ; ++k){
      id[k] = id[( *nd - 1 )];
    }
  }

  // If sh[0] > 0, shift towards larger indices
  if( sh[0] > 0 ){
    for( k = ( *nd - 1 ) ; k >= sh[0] ; --k ){
      id[k] = id[ k - sh[0] ];
    }
    // The first value is repeated in the first sh[0] points
    for( k = ( sh[0] - 1 ) ; k > 0 ; --k ){
      id[k] = id[0];
    }
  }

  // Add the shift that was already done to sh[1]
  sh1 = sh[1] - sh[0];

  // If sh1 < 0  we are supposed to shift 
  // the falling edges towards smaller indices
  if( sh1 < 0 ){
    ncut = 0;
    for( k = ( *nd - 1 ) ; k >= 0 ; --k ){
      if( id[ k ] == 0 ){
        ncut = 0;
      }else{
        --ncut;
      }
      if( ncut >= sh1 ) id[k] = 0;
    }
  }
  // If sh1 > 0 we are supposed to shift
  // the falling edges towards larger indices
  if( sh1 > 0 ){
    nadd = 0;
    for( k = 0 ; k < *nd ; ++k ){
      if( id[ k ] == 0 ){
        ++nadd;
      }else{
        nadd = 0;
      }
      if( nadd <= sh1 ) id[k] = 1;
    }

  }

  // Now there may be errors in the very end of the index
  //  vector, correct using the stored index lasttrue
  for( k = ( lasttrue + sh[1] + 1 ) ; k < *nd ; ++k ){
    id[k] = 0;
  }

  // Remove protection from the return value
  UNPROTECT(1);

  // Return the variable success only, the data is stored
  // in the R vectors 'cdata', 'idatar', and 'idatai'
  return(success);

}
