// file:fishs_add_clutter.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"
/* 

   A special version of fisher solver for ground clutter
   estimation. Assumes unit variance and adds only one
   row at a time. 

   Arguments:
    Qvec  Upper triangular part of Fisher information matrix
    yvec  Modified measurement vector
    arow  One row of theory matrix
    meas  Measurement
    nx    Number of unknowns

*/

void fishs_add_clutter( const SEXP Qvec , const SEXP yvec , Rcomplex * arow , Rcomplex * meas , const int nx )
{
  Rcomplex *q = COMPLEX(Qvec);
  Rcomplex *y = COMPLEX(yvec);
  int n  = nx;
  int i  = 0;
  int j  = 0;

  Rcomplex * restrict qtmp;
  Rcomplex * restrict acpy = arow;
  Rcomplex * restrict atmp;
  Rcomplex * restrict ytmp;
  Rcomplex * restrict mcpy = meas;

  // Pointers to y-vector and Fisher information matrix
  ytmp = y;
  qtmp = q;

  // Go through all range gates
  for( i = 0 ; i < n ; ++i ){
    
    // Second pointer to the theory matrix
    atmp = acpy;
    
    // Go through all columns in the upper triangular part
    for( j = 0 ; j < ( n - i ) ; ++j ){

      // Add information
      qtmp->r += ( acpy->r * atmp->r + acpy->i * atmp->i );
      qtmp->i += ( acpy->r * atmp->i - acpy->i * atmp->r );

      // Increment the second theory matrix counter
      ++atmp;

      // Increment the information matrix counter
      ++qtmp;

    }

    // Add the corresponding measurement to the y-vector
    ytmp->r += ( mcpy->r * acpy->r + mcpy->i * acpy->i );
    ytmp->i += ( mcpy->i * acpy->r - mcpy->r * acpy->i );
    
    // Increment the y-vector counter
    ++ytmp;
    
    // Increment the theory matrix counter
    ++acpy; 
    
  }
  
}


