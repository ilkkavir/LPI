// file:prepare_data.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/*
  Frequency mixing, index adjustments,
  and filtering in a single function

  Arguments:
   cdata     Complex voltage data vector
   idata     Integer vector of usable data indices
   ndata     Data vector length
   frequency Frequency offset
   shifts    Corrections to idata
   nup       Upsaling factor in resampling
   nfilter   Filter length (for upsampled data, final 
             filter length is nfilter / nup)
   nfirst    Decimation start index
   ipartial  Logical, are partial matches of
             idata with the filter accepted?

  Returns:
   ans       A list with elements
              cdata   Final complex data vector
              idata   Final index vector
              ndata   Final data vector length
	      success Logical, set if all processing
                      was successfull

 */


SEXP prepare_data( SEXP cdata , SEXP idata , SEXP ndata , SEXP frequency, SEXP shifts , SEXP nup , SEXP nfilter , SEXP nfirst , SEXP nfirstfrac , SEXP ipartial )		  
{
  SEXP ans;
  SEXP cdata_new;
  SEXP idata_new;
  SEXP ndata_new;
  SEXP s;
  SEXP names;
  char *cnames[4] = {"cdata","idata","ndata","success"};
  Rcomplex * restrict cnew;
  Rcomplex * restrict cold;
  int * restrict inew;
  int * restrict iold;
  uint64_t k;
  PROTECT_INDEX cpind=0;
  PROTECT_INDEX ipind=0;



  // Output list ans[[1]] = cdata ans[[2]] = pdata ,
  // ans[[3]] = idata , ans[[4]] = ndata , ans[[5]] = success
  PROTECT( ans = allocVector( VECSXP , 5 ) );

  // Allocate the new complex vector
  PROTECT_WITH_INDEX( cdata_new = allocVector( CPLXSXP , *(INTEGER(ndata)) ) , &cpind );

  // Allocate the new logical vector
  PROTECT_WITH_INDEX( idata_new = allocVector( LGLSXP , *(INTEGER(ndata)) ) , &ipind );

  // Allocate the new ndata variable
  PROTECT( ndata_new = allocVector( INTSXP , 1 ) );

  // A pointer to the new cdata vector
  cnew = COMPLEX( cdata_new );

  // A pointer to the old cdata vector
  cold = COMPLEX( cdata );

  // A pointer to the new idata vector
  inew = LOGICAL( idata_new );

  // A pointer to the old idata vector
  iold = LOGICAL( idata );

  // Copy data from old cdata to new cdata
  for( k = 0 ; k < *(INTEGER(ndata)) ; ++k ){
    cnew[k].r = cold[k].r;
    cnew[k].i = cold[k].i;
  }

  // Copy data from old idata to new idata
  for( k = 0 ; k < *(INTEGER(ndata)) ; ++k ){
    inew[k] = iold[k];
  }

  // Use the same pointers to copy old ndata to new ndata
  inew = INTEGER( ndata_new);
  iold = INTEGER( ndata );
  *inew = *iold;

  // The  success logical
  PROTECT( s = allocVector( LGLSXP , 1 ) );

  // Frequency mixing
  s = mix_frequency( cdata_new , ndata_new , frequency );

  // Index adjustments
  s = index_adjust( idata_new , ndata_new , shifts );

  // Filtering
  s = resample( cdata_new , idata_new , ndata_new , nup , nfilter , nfirst , nfirstfrac , ipartial );

  // Set cdata_new to zero at all points where idata_new==0
  inew = LOGICAL( idata_new );
  for( k = 0 ; k < *INTEGER(ndata_new) ; ++k ){
    if( inew[k] == 0 ){
      cnew[k].r = .0;
      cnew[k].i = .0;
    }
  }

  // Reallocate the vectors to match with the new data length
  SET_LENGTH( cdata_new , *INTEGER(ndata_new) );
  REPROTECT( cdata_new , cpind );
  SET_LENGTH( idata_new , *INTEGER(ndata_new) );
  REPROTECT( idata_new , ipind );

  // Collect the data into the return list
  SET_VECTOR_ELT( ans , 0 , cdata_new );
  SET_VECTOR_ELT( ans , 1 , idata_new );
  SET_VECTOR_ELT( ans , 2 , ndata_new );
  SET_VECTOR_ELT( ans , 3 , s );

  // Set the name attributes
  PROTECT( names = allocVector( STRSXP , 4 ));
  SET_STRING_ELT( names , 0 , mkChar( cnames[0] ) );
  SET_STRING_ELT( names , 1 , mkChar( cnames[1] ) );
  SET_STRING_ELT( names , 2 , mkChar( cnames[2] ) );
  SET_STRING_ELT( names , 3 , mkChar( cnames[3] ) );
  setAttrib( ans , R_NamesSymbol , names);

  UNPROTECT(6);

  return(ans);

}
