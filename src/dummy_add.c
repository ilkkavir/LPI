// file:dummy_add.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/* 
   Simple variance- and power-weighted average lag profile.
   Works only below one IPP range.

   Arguments:
    msum  Sum of normalised measurements
    vsum  sum of normalised inverse variances
    rmin  Lower edge of the measurement
    rmax  Upper edge
    mdata Complex measurement vector (lag profile)
    mamb  Complex range ambiguity function
    iamb  Range ambiguity function indices
    iprod Lagged product indices
    edata Measurement variances
    ndata Data vector length

   Returns:
    success 1 if the processing was succesful, 0 otherwise

*/

SEXP dummy_add( SEXP msum , SEXP vsum , SEXP rmin , SEXP rmax , SEXP mdata , SEXP mamb , SEXP iamb , SEXP iprod , SEXP edata , SEXP ndata )
{
  Rcomplex *ms = COMPLEX(msum);
  double *vs = REAL(vsum);
  int r1 = *INTEGER(rmin);
  int r2 = *INTEGER(rmax);
  Rcomplex *cd = COMPLEX(mdata);
  Rcomplex *ad = COMPLEX(mamb);
  int *ia = LOGICAL(iamb);
  int *ip = LOGICAL(iprod);
  double *vd = REAL(edata);
  int nd = *INTEGER(ndata);

  int i, j, r, r0;

  SEXP                success;
  int      * restrict i_success;

  // success output
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // local pointer to the success output
  i_success = LOGICAL( success );

  // set the success output
  *i_success = 1;

  // Skip first r2 points, their range ambiguity function
  // is not known
  r = r2+1;
  r0 = 0;

  // Walk through the data vector
  for( i = 0 ;  i < nd ; ++i ){
    // Check that we are above r1
    if( r >= r1 ){
      // Check that we are below r2
      if( r < r2 ){
	// Check that the point is flagged as usable
	if(ip[i]){
	  // The average vector starts from range r1
	  j = r-r1;
	  // Divide the lagged product with its variance and 
	  // multiply with TX power
	  ms[j].r += cd[i].r  / vd[i] * ad[r0].r;
	  ms[j].i += cd[i].i  / vd[i] * ad[r0].r;
	  // Inverse of variance scaled accordingly
	  vs[j] += ad[r0].r * ad[r0].r / vd[i];
	}
      }
    }

    // If a new pulse is transmitted set range to zero, 
    // otherwise increment the range counter. 
    if( ia[i] ){
      r = 0;
      r0 = i;
    }else{
      ++r;
    }
  }

  UNPROTECT(1);

  return(success);

}
