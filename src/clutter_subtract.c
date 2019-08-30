// file:clutter_subtract.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"
/*

  Ground clutter suppression. This function subtracts clutter
  signal from data.

  Arguments:
   tcdata  Complex transmitter samples
   tidata  Transmitter sample indices
   rcdata  Complex receiver samples
   ridata  Receiver sample indices
   ndata   Data vector length
   rmin    Minimum range
   rmax    Maximum range
   cldata  Measured clutter signal profile

  Returns:
   nrow    Number of points at which clutter
           signal was suppressed

*/

SEXP clutter_subtract( const SEXP tcdata , const SEXP tidata , const SEXP rcdata , const SEXP ridata , const SEXP ndata , const SEXP rmin , const SEXP rmax , const SEXP cldata )
{
  Rcomplex *tcd = COMPLEX( tcdata );
  int *tid = LOGICAL( tidata );
  Rcomplex * rcd = COMPLEX( rcdata );
  int *rid = LOGICAL( ridata );
  Rcomplex *cld = COMPLEX( cldata );
  const int nd = *INTEGER( ndata );
  const int r0 = *INTEGER( rmin );
  const int r1 = *INTEGER( rmax );

  int i;
  int j;
  int k;
  int r;
  int isum;
  int nx;
  SEXP nrow;
  int nr;
  Rcomplex clsum;
  Rcomplex * tcd2;
  Rcomplex * cld2;

  // Output
  PROTECT( nrow = allocVector( INTSXP , 1 ) );

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
  for( i = r1 ; i < ( nd - nx )  ; ++i ){
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
	    // Calculate clutter contribution and subtract it
	    clsum.r = 0.;
	    clsum.i = 0.;
	    tcd2 = tcd;
	    cld2 = cld;
	    for( j = 0 ; j < nx ; ++j ){
	      clsum.r += tcd2->r * cld2->r - tcd2->i * cld2->i;
	      clsum.i += tcd2->r * cld2->i + tcd2->i * cld2->r;
	      ++tcd2;
	      ++cld2;
	    }
	    rcd->r -= clsum.r;
	    rcd->i -= clsum.i;
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
