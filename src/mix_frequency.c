// file:mix_frequency.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/*
  Frequency mixing for IQ data

  This function allocates new vectors

  Argumnets:
   cdata      ndata complex vector of data samples
   ndata      Number of samples in cdata
   frequency  The mixing frequency

  Returns:
   ans       A list with elements
              cdata   Complex data samples after
                      frequency mixing
	      success Logical, set if all processing
                      was successful
*/


SEXP mix_frequency_R( SEXP cdata , SEXP ndata , SEXP frequency )
{
  SEXP ans;
  SEXP cdata_new;
  SEXP s;
  SEXP names;
  char *cnames[2] = {"cdata","success"};
  Rcomplex *cnew;
  Rcomplex *cold;
  register uint64_t k;


  // Output list ans[[1]] = cdata , ans[[2]] = success
  PROTECT( ans = allocVector( VECSXP , 2 ) );

  // Allocate the new complex vector
  PROTECT( cdata_new = allocVector( CPLXSXP , *(INTEGER(ndata)) ) );

  // A pointer to the new data vector
  cnew = COMPLEX( cdata_new );

  // A pointer to the old data vector
  cold = COMPLEX( cdata );

  // Copy data from old to new
  for( k = 0 ; k < *(INTEGER(ndata)) ; ++k ){
    cnew[k].r = cold[k].r;
    cnew[k].i = cold[k].i;
  }

  // The success logical
  PROTECT( s = allocVector( LGLSXP , 1 ) );

  // The actual frequency mixing
  s = mix_frequency( cdata_new , ndata , frequency );

  // Collect the data into the return list
  SET_VECTOR_ELT( ans , 0 , cdata_new );
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
  Frequency mixing for IQ data

  This function overwrites the cdata vector

  Argumnets:
   cdata      ndata complex vector of data samples
   ndata      Number of samples in cdata
   frequency  The mixing frequency

  Returns:
   success    1 if all processing was successful, 0 otherwise
*/
SEXP mix_frequency( SEXP cdata , SEXP ndata , SEXP frequency)
{
  // Pointers to the R variables
  Rcomplex *cd = COMPLEX(cdata);
  int *nd = INTEGER(ndata);
  double *fr = REAL(frequency);
  register uint64_t k, nc;
  double arg;
  Rcomplex ctmp;
  // Temporary variables
  int ncycle;
  double tmpprod;
  double idiff;
  double *coefr;
  double *coefi;
  // For the return value
  SEXP success;
  int *isuccess;

  // Allocate the return value and initialise it
  PROTECT(success = allocVector(LGLSXP,1));
  isuccess = LOGICAL(success);
  *isuccess = 1;

  // The multiplicand will be cyclic, find the cycle length
  ncycle = *nd;
  for( k = 1 ; k < *nd ; ++k){
    tmpprod = *fr * (double)(k);
    idiff = tmpprod - (double)((int)(tmpprod));
    if( fabs(idiff) <= FLT_MIN ){
      ncycle = k;
      break;
    }
  }

  // If the cycle length is one, the mixing would not change anything
  if( ncycle == 1 ){
    UNPROTECT(1);
    return(success);
  }

  // Tabulate the cyclic coefficients.
  // This usually saves time as radar engineers tend to
  // select nice numerical values for the frequencies
  coefr = (double*) R_Calloc( ncycle , double );
  coefi = (double*) R_Calloc( ncycle , double );
  for( k = 0 ; k < ncycle ; ++k ){
    arg      = 2.0 * M_PI * *fr * (double)(k);
    coefr[k] = cos(arg);
    coefi[k] = sin(arg);
  }

  // Actual mixing
  nc = 0;
  for( k = 0 ; k < *nd ; ++k ){
    ctmp.r = cd[k].r;
    ctmp.i = cd[k].i;
    cd[k].r = ctmp.r * coefr[nc] - ctmp.i * coefi[nc];
    cd[k].i = ctmp.i * coefr[nc] + ctmp.r * coefi[nc];
    ++nc;
    if( nc == ncycle ) nc = 0;
  }

  // Free the memory allocated for the coefficient tables
  Free(coefr);
  Free(coefi);

  // Remove protection from the return value
  UNPROTECT(1);

  // Return the variable success only, the data is stored in
  // the R vectors 'cdata', 'idatar', and 'idatai'
  return(success);
  
}

