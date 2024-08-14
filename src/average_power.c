// file:average_power.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/*
  Average power vector for variance estsimation

  The algorithm proceeds as follows

  1. Locate falling edges of pulses  from idatatx
  2. locate the first falling edge at least maxrange samples
     from the beginning, give this pulse the pulse index 0
  3. Pick maxrange samples from idatatx from  immediately 
     *before* the first falling edge
  4. At all other falling edges, compare the maxrange points
     before the edge with the samples picked in (3)
  5. If the vectors compared in (4) are identical, also this
     pulse is given pulse index 0, repeat for all pulses
  6. If there pulses are left without an index, select the
     first of them and repeat steps (4) and (5) to give 
     these pulses the index 1. 
  7. Continue with indices 2, 3, ... 
     until all pulses have an index
  8. When all pulses have indices, calculate average 
     power profiles from pulses with identical indices


  Arguments:
   cdata    Complex receiver samples
   idatatx  Transmitter sample index vector
   idatarx  Receiver sample index vector
   ndata    Number of points in data vectors
   maxrange Maximum range for power profile estimation 
   nminave  Minimum number of samples to be averaged

  Returns:
   pdata    Average power vector. The first element contains the 
            ratio largest pulse index / number of pulses.

 */

SEXP average_power( SEXP cdata , SEXP idatatx , SEXP idatarx , SEXP ndata , SEXP maxrange , SEXP nminave )
{
  Rcomplex * cd = COMPLEX( cdata );
  int * idtx = LOGICAL( idatatx );
  int * idrx = LOGICAL( idatarx );
  int nd = *INTEGER( ndata );
  int maxr = *INTEGER( maxrange );
  int nmin = *INTEGER( nminave );

  SEXP pdata;
  double *pd;
  double *ptmp;
  int *pedges;
  int nedges;
  int *pinds;
  int *nsamp;
  int k, i, j;
  int pindcur;
  int pindmax;
  int p1;
  int sameamb;
  int r;
  int ippend;
  int ntot;
  double ptot;

  ntot = 0;
  ptot = .0;

  // Inspect the TX index vector
  // to make sure that 1 is exactly 1
  for( k = 0 ; k < nd ; ++k ) idtx[ k ] = idtx[ k ] ? 1 : 0 ;

  // Allocate the power vector
  PROTECT( pdata = allocVector( REALSXP , nd ) );

  // A pointer to the power vector
  pd = REAL( pdata );

  // Initialise to zero
  for( k = 0 ; k < nd ; ++k ) pd[ k ] = 0.;

  // Allocate a temporary vector for 
  // power profile calculation
  ptmp = R_Calloc( nd , double );

  // Initialise to zero
  for( k = 0 ; k < nd ; ++k ) ptmp[ k ] = 0.;

  // Allocate a vector for sample counter
  nsamp = R_Calloc( nd , int );

  // Initialise to zero
  for( k = 0 ; k < nd ; ++k ) nsamp[ k ] = 0;

  // Allocate a vector for pulse edge positions
  // (this could be shorter if needed)
  pedges = R_Calloc( nd , int );

  // Initialise to zero
  for( k = 0 ; k < nd ; ++k ) pedges[ k ] = 0;

  // Allocate a vector for pulse indices
  pinds = R_Calloc( nd , int );

  // Initialise to -1
  for( k = 0 ; k < nd ; ++k ) pinds[ k ] = -1;


  // Locate all falling edges of pulses
  nedges = 0;
  for( k = 0 ; k < ( nd - 1 ) ; ++k )
    {
      if( idtx[ k ] )
	{
	  if( !(idtx[ k + 1] ) )
	    {
	      pedges[ nedges++ ] = k;
	    }
	}
    }

  // The first falling pulse edge at least
  // maxr samples from the beginning
  p1 = nedges;
  for( k = 0 ; k < nedges ; ++k )
    {
      if( pedges[ k ] > maxr )
	{
	  p1 = k;
	  break;
	}
    }

  // Inspect the tx indices and give a unique index for
  // each unique 0-lag range-ambiguity function
  pindcur = 0;
  for( k = p1 ; k < nedges ; ++k )
    {
      // pinds < 0 for pulses that do not yet have an index
      if( pinds[ k ] < 0 )
	{
	  // Go through all the pulses
	  for( i = k ; i < nedges ; ++i )
	    {
	      // Compare only with pulses that 
	      // do not yet have an index
	      if( pinds[ i ] < 0 )
		{
		  // Inspect the points just before this pulse
		  sameamb = 1;
		  for( j = 0 ; j < maxr ; ++j )
		    {
		      if( (idtx[ pedges[ k ] - j ]) != (idtx[ pedges[ i ] - j ]) )
			{
			  sameamb = 0;
			  break;
			}
		    }
		  // If the ambiguities were identical,
		  // assign the pulse with the index pindcur
		  if( sameamb ) pinds[ i ] = pindcur;
		}
	    }
	  // Increment pindcur
	  ++pindcur;
	}
    }

  // There may be a pulse / pulses without an index
  // in the begin of data vector.
  // Give them an index if possible
  if( p1 > 0 )
    {
      for( i = p1 ; i < nedges ; ++i )
	{
	  sameamb = 1;
	  for( j = 0 ; j < pedges[ p1 - 1 ] ; ++j )
	    {
	      if( idtx[ pedges[ p1 - 1 ] - j ] != idtx[ pedges[ i ] - j ] )
		{
		  sameamb = 0;
		  break;
		}
	    }
	  if( sameamb )
	    {
	      pinds[ p1 - 1 ] = pinds[ i ];
	      break;
	    }
	}
      // Give a new index for the  pulse p1-1 if it did not match
      // with any of the exisiting ones. Pulses before p1-1 will 
      // not be used and they do not need an index.
      if( pinds[ p1 - 1 ] < 0 ) pinds[ p1 - 1 ] = pindcur;
    }

  // Store the largest pind
  pindmax = pindcur;

  // We have now an index for each pulse that needs one. Pulses
  // with equal indices have similar power profile range ambiguity
  // functions and their signal powers can be averaged.
  // Now we will walk through all different pulse indices,
  // calculate the correspondign power-profiles, and
  // store the results in appropriate places in the average
  // power vector

  // Start from the first falling edge, or
  // one point before if necessary
  if( p1 > 0 ) --p1 ;

  // Go through all pulses
  for( k = p1 ; k < nedges ; ++k )
    {

      // The indices will be set to -1 after processing,
      // an index >= indicates that the point has not 
      // yet been processed
      if( pinds[ k ] >= 0 )
	{

	  // Initialise the temporary power vector to zero
	  for( i = 0 ; i < nd ; ++i ) ptmp[ i ] = 0.;
	  
	  // Initialise the sample counter to zero
	  for( i = 0 ; i < nd ; ++i ) nsamp[ i ] = 0 ;
	  
	  // Check remaining pulses and try to find
	  // the same index
	  for( j = k ; j < nedges ; ++j )
	    {
	      // If a matching index is found, add power from the
	      // ipp to the temporary profile and increment sample
	      // counter accordingly
	      if( pinds[ j ] == pinds[ k ] )
		{

		  // Find distance to the next pulse end  (must not
		  // stop at pulse start in order to facilitate 
		  // bistatic operation)
		  if( ( j + 1 ) >= nedges )
		    {
		      ippend = nd - pedges[ j ];
		    }
		  else
		    {
		      ippend = pedges[ j + 1 ] - pedges[ j ];
		    }
		  for( i = 0 ; i < ippend ; ++i )
		    {
		      r = pedges[ j ] + i;
		      // This cuts off points that are too close to 
		      // the beginning of the data vector
		      if( r >= maxr )
			{
			  if( idrx[ r ] )
			    {
			      ptmp[ i ]  += cd[ r ].r * cd[ r ].r + cd[ r ].i * cd[ r ].i;
			      nsamp[ i ] += 1;
			      ptot +=  cd[ r ].r * cd[ r ].r + cd[ r ].i * cd[ r ].i;
			      ++ntot;
			    }
			}
		    }
		}
	    }

	  // Divide the summed powers by 
	  // the number of summed samples
	  for( i = 0 ; i < nd ; ++i )
	    {
	      if( nsamp[ i ] >= nmin ){
		ptmp[ i ] /= (double) nsamp[ i ];
	      }else{
		ptmp[ i ] = -1.;
	      }
	    }

	  // Go through the indices again and copy the power
	  //  values to appropriate places Set pinds to -1 at
	  //  points that have already been visited
	  pindcur = pinds[ k ];
	  for( j = k ;  j < nedges ; ++j )
	    {
	      if( pinds[ j ] == pindcur )
		{
		  if( ( j + 1 ) >= nedges )
		    {
		      ippend = nd - pedges[ j ];
		    }
		  else
		    {
		      ippend = pedges[ j + 1 ] - pedges[ j ];
		    }
		  for( i = 0 ; i < ippend ; ++i )
		    {
		      r = pedges[ j ] + i;
		      pd[ r ] = ptmp[ i ];
		    }
		  pinds[ j ] = -1;
		}
	    }

	}

    }

  // Put the grand average power to points that did not have 
  // enough averaged samples (they are set to -1 at this point)
  ptot /= (float)ntot;
  for( i = 0 ;  i < nd ; ++i ){
    if( pd[ i ] < 0.) pd[ i ] = ptot;
  }
  
  // Store the ratio pindmax / nedges to the first data point.
  // If the ratio is large the power estimation will not perform
  // very well.
  // The power value in this point cannot ever be needed in LPI.
  pd[0] = (float)pindmax / (float)nedges;

  // Free the temporary allocations
  Free(ptmp);
  Free(nsamp);
  Free(pinds);
  Free(pedges);
 
  UNPROTECT(1);
  return(pdata);

}
