// file:deco_add.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/* 
   Matched filter decoding, modified from fishs_add.

   Arguments:
    Qvec  Diagonal of the Fisher information matrix
    yvec  Modified measurement vector
    arows Theory matrix rows
    irows Indices of non-zero theory matrix elements    
    meas  Measurements
    var   Measurement variances
    nx    Number of unknowns
    nrow  Number of theory rows in arows

   Returns:
    success 1 if the processing was succesful, 0 otherwise

*/

SEXP deco_add(  SEXP Qvec , SEXP yvec , const SEXP arows , const SEXP irows , const SEXP meas  , const SEXP var  , const SEXP nx   , const SEXP nrow  )
{
  Rcomplex *q = COMPLEX(Qvec);
  Rcomplex * restrict qtmp;

  Rcomplex *y = COMPLEX(yvec);
  Rcomplex * restrict ytmp;

  Rcomplex * restrict acpy = COMPLEX(arows);

  int * restrict icpy = LOGICAL(irows);  

  Rcomplex * restrict mcpy = COMPLEX(meas);

  double   * restrict vcpy = REAL(var);

  int n  = *INTEGER(nx);

  int nr = *INTEGER(nrow);

  int i  = 0;
  int j  = 0;
  int l  = 0;

  SEXP success;
  int * restrict i_success;

  // Success output
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // Local pointer to the success output
  i_success = LOGICAL( success );

  // Set the success output
  *i_success = 1;

  // Go through all theory matrix rows
  for( l = 0 ; l < nr ; ++l ){

    // Pointers to y-vector and Fisher information matrix diagonal
    ytmp = y;
    qtmp = q;

    // Go through all range gates
    for( i = 0 ; i < n ; ++i ){

      if ( *icpy ) {
	
	// Add information (only diaonal)
	qtmp->r += ( acpy->r * acpy->r + acpy->i * acpy->i ) / *vcpy;
	qtmp->i += ( acpy->r * acpy->i - acpy->i * acpy->r ) / *vcpy;

      }
      
      // Increment information matrix counter (only diagonal)
      ++qtmp;

      // Add the corresponding measurement to the y-vector
      ytmp->r += ( mcpy->r * acpy->r + mcpy->i * acpy->i ) / *vcpy;
      ytmp->i += ( mcpy->i * acpy->r - mcpy->r * acpy->i ) / *vcpy;

      // Increment the y-vector counter
      ++ytmp;

      // Increment the theory matrix counter
      ++acpy;
      ++icpy;

    }

    // Increment the variance and measurement vector counters
    ++mcpy;
    ++vcpy;

  }

  UNPROTECT(1);

  return(success);

}


