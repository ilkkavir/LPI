// file:fishs_add.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/* 
   Inverse problem solver using direct calculation of the 
   Fisher information matrix. Data accumulation.

   Arguments:
    Qvec  Upper triangular part of the Fisher
          information matrix as a vector
    yvec  Modified measurement vector
    arows Theory matrix rows
    irows Indices of non-zero theory matrix elements
    meas  Measurements
    var   Measurement variances
    nx    Number of unknowns
    nrow  Number of theory rows in arows

   Returns:
    success 1 if the processing was successful, 0 otherwise

*/

SEXP fishs_add( const SEXP Qvec , const SEXP yvec , const SEXP arows , const SEXP irows , const SEXP meas  , const SEXP var  , const SEXP nx   , const SEXP nrow  )               
{
  Rcomplex *q = COMPLEX(Qvec);
  Rcomplex *y = COMPLEX(yvec);
  int n  = *INTEGER(nx);
  int nr = *INTEGER(nrow);
  int i  = 0;
  int j  = 0;
  int l  = 0;

  Rcomplex * restrict qtmp;
  Rcomplex * restrict acpy = COMPLEX(arows);
  Rcomplex * restrict atmp;
  int * restrict icpy = LOGICAL(irows);  
  int * restrict itmp;  
  Rcomplex * restrict ytmp;
  Rcomplex * restrict mcpy = COMPLEX(meas);
  double   * restrict vcpy = REAL(var);

  SEXP                success;
  int      * restrict i_success;

  // success output
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // local pointer to the success output
  i_success = LOGICAL( success );

  // set the success output
  *i_success = 1;

  // Go through all theory matrix rows
  for( l = 0 ; l < nr ; ++l ){

    // Pointers to y-vector and Fisher information matrix
    ytmp = y;
    qtmp = q;

    // Go through all range gates
    for( i = 0 ; i < n ; ++i ){

      // Second pointer to the theory matrix
      atmp = acpy;
      itmp = icpy;

      if ( *icpy ) {
	// Go through all columns in the upper triangular part
	for( j = 0 ; j < ( n - i ) ; ++j ){
	  
	  // Add information

	  if( *itmp ){
	    qtmp->r += ( acpy->r * atmp->r + acpy->i * atmp->i ) / *vcpy;
	    qtmp->i += ( acpy->r * atmp->i - acpy->i * atmp->r ) / *vcpy;
	  }
	  
	  // Increment the second theory matrix counter
	  ++atmp;
	  ++itmp;
	  
	  // Increment the information matrix counter
	  ++qtmp;
	  
	}
	
	
	// Add the corresponding measurement to the y-vector
	ytmp->r += ( mcpy->r * acpy->r + mcpy->i * acpy->i ) / *vcpy;
	ytmp->i += ( mcpy->i * acpy->r - mcpy->r * acpy->i ) / *vcpy;
	
	// Increment the y-vector counter
	++ytmp;
	
      }else{
	// Jump to the next diagonal element in q
	qtmp += n-i;
	++ytmp;
      }

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


