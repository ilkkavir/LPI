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

SEXP fishs_add2( SEXP QvecR , SEXP QvecI , SEXP yvecR , SEXP yvecI , const SEXP arowsR , const SEXP arowsI , const SEXP irows , const SEXP measR , const SEXP measI  , const SEXP var  , const SEXP nx   , const SEXP nrow  )               
{
  double *qR = REAL(QvecR);
  double *qI = REAL(QvecI);
  double * restrict qtmpR;
  double * restrict qtmpI;

  double *yR = REAL(yvecR);
  double *yI = REAL(yvecI);
  double * restrict ytmpR;
  double * restrict ytmpI;

  double *acpyR = REAL(arowsR);
  double *acpyI = REAL(arowsI);
  double *atmpR;
  double *atmpI;

  int *icpy = LOGICAL(irows);  
  int *itmp;  

  double * restrict mcpyR = REAL(measR);
  double * restrict mcpyI = REAL(measI);

  double * restrict vcpy = REAL(var);

  int n  = *INTEGER(nx);

  int nr = *INTEGER(nrow);

  int i = 0;
  int j = 0;
  int l = 0;
  int k = 0;
  int addlines = 0;
  int naddlines  = 0;

  SEXP success;
  int * restrict i_success;

  // success output
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // local pointer to the success output
  i_success = LOGICAL( success );

  // set the success output
  *i_success = 1;

  // Go through all theory matrix rows
  for( l = 0 ; l < nr ; ++l ){

    // Pointers to y-vector and Fisher information matrix
    ytmpR = yR;
    ytmpI = yI;
    qtmpR = qR;
    qtmpI = qI;
    // Go through all range gates
    for( i = 0 ; i < n ; ++i ){

      // Second pointer to the theory matrix
      atmpR = acpyR;
      atmpI = acpyI;
      itmp = icpy;

      if ( *icpy ) {

        // Go through all columns in the upper triangular part
	// Check in blocks of 8 and skip those that contain only zeros.
        for( j = 0 ; j < ( n - i ) ; j+=8 ){
	  
	  addlines = 0;
	  naddlines = 0;
	  for ( k = 0 ; ( k < 8 ) & ((k+j) < ( n - i )) ; ++k ){
	    addlines += *itmp;
	    ++itmp;
	    ++naddlines;
	  }
	  
	  if (addlines){
#pragma GCC ivdep
	    for (k = 0 ; k<naddlines ; ++k ){
	      // Add information
	      
	      //if( *itmp ){
	      *qtmpR += ( *acpyR * *atmpR + *acpyI * *atmpI ) / *vcpy;
	      *qtmpI += ( *acpyR * *atmpI - *acpyI * *atmpR ) / *vcpy;
	      /* *qtmpR += ( *acpyR * *atmpR + *acpyI * *atmpI ) / *vcpy; */
	      /* *qtmpI += ( *acpyR * *atmpI - *acpyI * *atmpR ) / *vcpy; */
	      
	      // Use the return value as a flop counter for testing. Will overflow in many cases...
	      *i_success += 10;
	      //}
	      
	      
	      // Increment the second theory matrix counter
	      ++atmpR;
	      ++atmpI;
	      
	      // Increment the information matrix counter
	      ++qtmpR;
	      ++qtmpI;
	    }
	  }else{
	    atmpR += naddlines;
	    atmpI += naddlines;
	    qtmpR += naddlines;
	    qtmpI += naddlines;
	  }
	  
        }
	
	
        // Add the corresponding measurement to the y-vector
	*ytmpR += ( *mcpyR * *acpyR + *mcpyI * *acpyI ) / *vcpy;
        *ytmpI += ( *mcpyI * *acpyR - *mcpyR * *acpyI ) / *vcpy;
	
        // Use the return value as a flop counter for testing. Will overflow in many cases...
        *i_success += 10;
        
        // Increment the y-vector counter
        ++ytmpR;
        ++ytmpI;
	
      }else{
        // Jump to the next diagonal element in q
        qtmpR += n-i;
        qtmpI += n-i;
        ++ytmpR;
        ++ytmpI;
      }

      // Increment the theory matrix counter
      ++acpyR;
      ++acpyI;
      ++icpy;

    }

    // Increment the variance and measurement vector counters
    ++mcpyR;
    ++mcpyI;
    ++vcpy;

  }

  UNPROTECT(1);

  return(success);

}
