// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/* 
   Matched filter decoder. With re and im in separate arrays.

   Arguments:
    QvecR  Diagonal of the precision matrix, real part (imaginary is always zero)
    yvecR  Modified measurement vector, real part
    yvecI  Modified measurement vector, imaginary part
    arowsR Theory matrix rows, real part
    arowsI Theory matrix rows, imaginary part
    irows  Indices of non-zero theory matrix elements
    measR  Measurements, real part
    measI  Measurements, imaginary part
    var    Measurement variances
    nx     Number of unknowns
    nrow   Number of theory rows in arows

   Returns:
    success 1 if the processing was successful, 0 otherwise

*/

SEXP decor_add( SEXP QvecR , SEXP yvecR , SEXP yvecI , const SEXP arowsR , const SEXP arowsI , SEXP irows , const SEXP measR , const SEXP measI  , const SEXP var  , const SEXP nx   , const SEXP nrow , SEXP flops )               
{
  double *qR = REAL(QvecR);
  double * restrict qtmpR;

  double *yR = REAL(yvecR);
  double *yI = REAL(yvecI);
  double * restrict ytmpR;
  double * restrict ytmpI;

  double *acpyR = REAL(arowsR);
  double *acpyI = REAL(arowsI);

  int *icpy = LOGICAL(irows);  

  double * restrict mcpyR = REAL(measR);
  double * restrict mcpyI = REAL(measI);

  double * restrict vcpy = REAL(var);

  int n  = *INTEGER(nx);

  int nr = *INTEGER(nrow);

  double *flop_count = REAL(flops);
  

  int i = 0;
  int k = 0;
  int l = 0;
  int addlines = 0;
  int naddlines  = 0;
  long int n_adds = 0;
  
  SEXP success;
  int * restrict i_success;

  double std;
  double * mtmpR;
  double * mtmpI;

  // success output
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // local pointer to the success output
  i_success = LOGICAL( success );

  // set the success output (will always be 1 at the moment..)
  *i_success = 1;


  



  // noise whitening (divide A and m with sqrt(var) )
  atmpR = acpyR;
  atmpI = acpyI;
  itmp = icpy;
  mtmpR = mcpyR;
  mtmpI = mcpyI;

  // Go through all theory matrix rows
  for( l = 0 ; l < nr ; ++l ){

    std = sqrt(*vcpy);
    
    // Go through all range gates
    for( i = 0 ; i < n ; ++i ){

      // divide only if this sample will be used
      if(*itmp){
	*atmpR = *atmpR / std;
	*atmpI = *atmpI / std;
      }
      
      // Increment the theory matrix counter
      ++atmpR;
      ++atmpI;
      ++itmp;
      
    }

    // divide the measurement with std
    *mtmpR = *mtmpR / std;
    *mtmpI = *mtmpI / std;
    
    // Increment the variance and measurement vector counters
    ++vcpy;
    ++mtmpR;
    ++mtmpI;

  }



  
  


  // Go through all theory matrix rows
  for( l = 0 ; l < nr ; ++l ){

    // Pointers to y-vector and Fisher information matrix
    ytmpR = yR;
    ytmpI = yI;
    qtmpR = qR;
    
    // Go through all range gates
    for( i = 0 ; i < n ; i+=8 ){
      
      // check if there are any non-zero data in the next 8 elements
      // naddlines is needed to avoid overflow at end of the vector
      addlines = 0;
      naddlines = 0;
      for ( k = 0 ; ( k < 8 ) & ( k < (n-i) ) ; ++k ){
	addlines += *icpy;
	++icpy;
	++naddlines;
      }
      
      // if there is something to add
      if (addlines){
	// add the information to real and imaginary parts of matrix Q
#pragma GCC ivdep
	for (k = 0 ; k<naddlines ; ++k ){
	  
	  // Add information, the imaginary part is always zero
	  *qtmpR += ( *acpyR * *acpyR + *acpyI * *acpyI ); // / *vcpy;
	  
	  // Add the corresponding measurement to the y-vector
	  *ytmpR += ( *mcpyR * *acpyR + *mcpyI * *acpyI ); // / *vcpy;
	  *ytmpI += ( *mcpyI * *acpyR - *mcpyR * *acpyI ); // / *vcpy;
	  
	  // Increment the information matrix and measurement vector counters
	  ++qtmpR;
	  ++acpyR;
	  ++acpyI;
	  ++ytmpR;
	  ++ytmpI;
	}
	// count added theory matrix elements and measurement rows
	n_adds += naddlines;
	
      }else{
	// move forward if only zeros were found
	qtmpR += naddlines;
	acpyR += naddlines;
	acpyI += naddlines;
	ytmpR += naddlines;
	ytmpI += naddlines;
      }
	
    }
    
    // Increment the variance and measurement vector counters
    ++mcpyR;
    ++mcpyI;
    //    ++vcpy;
    
  }

  // total number of floating point operations.
  //  *flop_count += 15.*((double)(n_adds));
  *flop_count += 12.*((double)(n_adds));
  
  UNPROTECT(1);

  return(success);

}
