// file:clutter_meas.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/*

  Ground clutter suppression. This function adds clutter
  signal measurements to an inverse problem. The function
  clutter_subtract subtracts clutter contribution from a
  signal.

  Arguments:
   tcdata  Complex transmitter samples
   tidata  Transmitter sample indices
   rcdata  Complex receiver samples
   ridata  Receiver sample indices
   ndata   Data vector length
   rmin    Minimum range
   rmax    Maximum range
   Qvec    Upper triangular part of Fisher information matrix
   yvec    Modified measurement vector

  Returns:
   nrow    Number of measurement rows in the inverse problem

*/
SEXP clutter_meas( const SEXP tcdata , const SEXP tidata , const SEXP rcdata , const SEXP ridata , const SEXP ndata , const SEXP rmin , const SEXP rmax , SEXP Qvec , SEXP yvec )
{
  Rcomplex *tcd = COMPLEX( tcdata );
  int *tid = LOGICAL( tidata );
  Rcomplex * rcd = COMPLEX( rcdata );
  int *rid = LOGICAL( ridata );
  const int nd = *INTEGER( ndata );
  const int r0 = *INTEGER( rmin );
  const int r1 = *INTEGER( rmax );

  int i;
  int r;
  int isum;
  int nx;
  SEXP nrow;
  int nr;

  // Output
  PROTECT( nrow = allocVector( INTSXP , 1 ) );

  // Make sure that the data vectors contain non-zero
  // values only at points in which the logical vectors
  // are not set
  for( i = 0 ; i < nd ; ++i){
    if( tid[i]==0 ){
      tcd[i].r = 0.0;
      tcd[i].i = 0.0;
    }
    if( rid[i]==0 ){
      rcd[i].r = 0.0;
      rcd[i].i = 0.0;
    }
  }

  // Initialization
  nr = 0;
  nx = r1 - r0 + 1;
  r = 0;
  isum = 0;
  // Sum tx indices and set r
  for( i = 0 ; i <= r1 ; ++i ){
    // The largest range is corresponds to index 0,
    // after nx samples we will be below rmin.
    if( i < nx ) isum += tid[i];
    // Increment r
    ++r;
    // Set r to zero if a transmitter sample is meat
    if( tid[i] ) r = 0;
    // increment the rx data pointer
    ++rcd;
  }

  // Go through all data points
  for( i = r1 ; i < nd ; ++i ){
    // Set r = 0 if a transmitter sample is meat
    if( tid[i] ) r = 0;
    // Are we below rmax?
    if( r <= r1 ){
      // Are we above rmin?
      if( r >= r0 ){
	// Are the pulses within the clutter ranges?
	if( isum ){
	  // Is this receiver sample usable?
	  if( rid[i] ){
	    // Add a measurement
	    fishs_add_clutter( Qvec , yvec , tcd , rcd , nx );
	    // Increment measurement row counter
	    ++nr;
	  }
	}
      }
    }
    // Update counters if this was not the last sample
    if( i < nd ){
      isum -= tid[ i - r1 ];
      isum += tid[ i - r0 + 1 ];
      ++r;
      ++rcd;
      ++tcd;
    }
  }

  // Copy the number of rows to output
  *INTEGER( nrow ) = nr;

  UNPROTECT(1);

  // Return number of measured rows
  return( nrow );

}
