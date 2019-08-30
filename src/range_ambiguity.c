// file:range_ambiguity.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/*
  Range ambiguity function with linear 
  interpolation of TX data

  Arguments:
   cdata1  First complex transmitter samples
   cdata2  Second complex transmitter samples
   idata1  First transmitter sample indices
   idata2  Seconds transmitter sample indices
   cdatap  Complex range ambiguity function
   idatap  Range ambiguity index vector
   ndata1  Length of vectors cdata1 and idata1
   ndata2  Length of vectors cdata2 and idata2
   lag     Lag

  Returns:
   success 1 if all processing was successful, 0 otherwise

*/

SEXP range_ambiguity( SEXP cdata1 , SEXP cdata2 , SEXP idata1 , SEXP idata2 , SEXP cdatap , SEXP idatap , SEXP ndata1 , SEXP ndata2 , SEXP lag )
{
  Rcomplex *cd1 = COMPLEX(cdata1);
  Rcomplex *cd2 = COMPLEX(cdata2);
  int *id1 = LOGICAL(idata1);
  int *id2 = LOGICAL(idata2);
  Rcomplex *cdp = COMPLEX(cdatap);
  int *idp =  LOGICAL(idatap);
  int nd1 = *INTEGER(ndata1);
  int nd2 = *INTEGER(ndata2);
  int l = *INTEGER(lag);
  SEXP success;
  int *isuccess;
  int k = 0;
  int npr;
  int ninterp = AMB_N_INTERP;
  int i;
  double * tmpr1;
  double * tmpi1;
  double * tmpr2;
  double * tmpi2;

  // Allocate temporary vectors for interpolated data
  tmpr1 = (double*) R_Calloc( 2*ninterp , double );
  tmpi1 = (double*) R_Calloc( 2*ninterp , double );
  tmpr2 = (double*) R_Calloc( 2*ninterp , double );
  tmpi2 = (double*) R_Calloc( 2*ninterp , double );

  // Output data length will be minimum of the 
  // two input data lengths, minus the lag
  npr = nd1 - l;
  if( nd1 > nd2 ) npr = nd2 - l;

  // Allocate the success return value
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // A local pointer to the success value
  isuccess = LOGICAL( success );
  *isuccess = 1;

  // The actual lagged product calculation
  for( k = 0 ; k < npr ; ++k ){
    // The index vector
    idp[k] = (id1[k] * id2[k+ l]);
    // Multiply data values only if the index vector was set
    if(idp[k]){
      // Initialize the temporary vectors to zero
      for( i = 0 ; i < ( 2 * ninterp ) ; ++i ){
	tmpr1[i] = .0;
	tmpi1[i] = .0;
	tmpr2[i] = .0;
	tmpi2[i] = .0;
      }
      // Linear interpolation towards the previous data point
      if( k > 1 ){
	for( i = 0 ; i < ninterp ; ++i ){
	  tmpr1[i] = cd1[k-1].r + ( cd1[k].r - cd1[k-1].r ) * ( 1. - (double)i / (double)( 2 * ninterp ) );
	  tmpi1[i] = cd1[k-1].i + ( cd1[k].i - cd1[k-1].i ) * ( 1. - (double)i / (double)( 2 * ninterp ) );
	  tmpr2[i] = cd2[k-1+l].r + ( cd2[k+l].r - cd2[k-1+l].r ) * ( 1. - (double)i / (double)( 2 * ninterp ) );
	  tmpi2[i] = cd2[k-1+l].i + ( cd2[k+l].i - cd2[k-1+l].i ) * ( 1. - (double)i / (double)( 2 * ninterp ) );
	}
      }
      // Linear interpolation towards the next data point
      if( k < npr ){
        for( i = 0 ; i < ninterp ; ++i ){
          tmpr1[i+ninterp] = cd1[k].r + ( cd1[k+1].r - cd1[k].r ) * ( (double)i / (double)( 2 * ninterp ) );
          tmpi1[i+ninterp] = cd1[k].i + ( cd1[k+1].i - cd1[k].i ) * ( (double)i / (double)( 2 * ninterp ) );
          tmpr2[i+ninterp] = cd2[k+l].r + ( cd2[k+1+l].r - cd2[k+l].r ) * ( (double)i / (double)( 2 * ninterp ) );
          tmpi2[i+ninterp] = cd2[k+l].i + ( cd2[k+1+l].i - cd2[k+l].i ) * ( (double)i / (double)( 2 * ninterp ) );
        }
      }
      // Initialize the final data value to zero
      cdp[k].r = .0;
      cdp[k].i = .0;
      // Add products of the interpolated data
      for( i = 0 ; i < ( 2 * ninterp ) ; ++i ){
        cdp[k].r += tmpr1[i] * tmpr2[i] + tmpi1[i] * tmpi2[i];
        cdp[k].i += tmpr1[i] * tmpi2[i] - tmpi1[i] * tmpr2[i];
      }
      // Divide with number of summed values
      cdp[k].r /= (double)(2*ninterp);
      cdp[k].i /= (double)(2*ninterp);
    }
  }

  // Set l index values from the beginning to false
  for( k = 0 ; k < l ; ++k ){
    idp[npr+k] = 0;
  }

  // Free the temporary vectors
  Free(tmpr1);
  Free(tmpi1);
  Free(tmpr2);
  Free(tmpi2);

  UNPROTECT(1);

  return(success);

}


