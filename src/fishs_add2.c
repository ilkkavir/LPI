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


/* A FASTER VERSION BELOW. THIS ONE MINIMIZES FLOPS, BUT APPARENTLY THE VARIABLE BLOCK SIZE SLOWS DOWN THE COMPUTATIONS */
/* 	naddlines = 0; */
/* 	j = 0; */
/* 	// check all elements in this row */
/* 	while ( j < ( n - i ) ){ */
/* 	  // the non-zero data are in continuous blocks due to the pulsed transmissions. Find length of the current block. */
/* 	  if(*itmp){ */
/* 	    ++naddlines; */
/* 	  }else{ */
/* 	    // add information from this block (pulse) */
/* 	    if (naddlines){ */
/* #pragma GCC ivdep */
/* 	      for (k = 0 ; k < naddlines ; ++k ){ */
/* 		*qtmpR += ( *acpyR * *atmpR + *acpyI * *atmpI ) / *vcpy; */
/* 		*qtmpI += ( *acpyR * *atmpI - *acpyI * *atmpR ) / *vcpy; */
/* 		++atmpR; */
/* 		++atmpI; */
/* 		++qtmpR; */
/* 		++qtmpI; */
/* 		*i_success += 10; */
/* 	      } */
/* 	      // the lines have been added, set naddlines to 0 */
/* 	      naddlines = 0; */
/* 	      // just increment the counters when zero-data are found.  */
/* 	    }else{ */
/* 	      ++atmpR; */
/* 	      ++atmpI; */
/* 	      ++qtmpR; */
/* 	      ++qtmpI; */
/* 	    } */
/* 	  } */
/* 	  ++j; */
/* 	  ++itmp; */
/* 	} */
/* 	// add information from pulses at the edge */
/* 	if (naddlines){ */
/* #pragma GCC ivdep */
/* 	  for (k = 0 ; k < naddlines ; ++k ){ */
/* 	    *qtmpR += ( *acpyR * *atmpR + *acpyI * *atmpI ) / *vcpy; */
/* 	    *qtmpI += ( *acpyR * *atmpI - *acpyI * *atmpR ) / *vcpy; */
/* 	    ++atmpR; */
/* 	    ++atmpI; */
/* 	    ++qtmpR; */
/* 	    ++qtmpI; */
/* 	    *i_success += 10; */
/* 	  } */
/* 	  // the lines have been added, set naddlines to 0 */
/* 	  naddlines = 0; */
/* 	} */


	// THE FASTER VERSION WITH CONSTANT BLOCK SIZE (8).

	
        // Go through all columns in the upper triangular part
	// Check in blocks of 8 and skip those that contain only zeros.
        for( j = 0 ; j < ( n - i ) ; j+=8 ){

	  // check if there are any non-zero data in the next 8 elements
	  // naddlines is needed to avoid overflow at end of the vector
	  addlines = 0;
	  naddlines = 0;
	  for ( k = 0 ; ( k < 8 ) & ((k+j) < ( n - i )) ; ++k ){
	    addlines += *itmp;
	    ++itmp;
	    ++naddlines;
	  }

	  // if there is something to add
	  if (addlines){
	    // add the information to real and imaginary parts of matrix Q
#pragma GCC ivdep
	    for (k = 0 ; k<naddlines ; ++k ){

	      // Add information
	      *qtmpR += ( *acpyR * *atmpR + *acpyI * *atmpI ) / *vcpy;
	      *qtmpI += ( *acpyR * *atmpI - *acpyI * *atmpR ) / *vcpy;
	      
	      // Use the return value as a flop counter for testing. Will overflow in many cases...
	      *i_success += 10;
	      
	      
	      // Increment the second theory matrix counter
	      ++atmpR;
	      ++atmpI;
	      
	      // Increment the information matrix counter
	      ++qtmpR;
	      ++qtmpI;
	    }
	  }else{
	    // move forward if only zeros were found
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
