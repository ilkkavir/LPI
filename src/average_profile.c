// file:average_profile.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/*
  Average lag-profile vector for speeding up
  the inversion process. Each average is 
  calculated over samples from  the same point in 
  the repeated code cycle

  The complicated structure is used because
  measuremnts may contain additional sync 
  times which need to be skipped.


  Arguments:
   cdata  Complex lagged product vector
   idata  Index vector for cdata
   ndata  Data vector length
   N_CODE Code cycle length

  Returns:
   success 1 if the processing was successful, 0 otherwise
 
*/

SEXP average_profile( SEXP cdata , SEXP idata , SEXP ndata , SEXP N_CODE)
{
  Rcomplex * cd = COMPLEX( cdata );
  int * id = LOGICAL( idata );
  int nd = *INTEGER( ndata );
  int ncode = *INTEGER( N_CODE );

  double *aver;
  double *avei;
  R_len_t *nave;
  R_len_t k;
  R_len_t ind1 , ind2, ipp_count;
  SEXP success;
  int *isuccess;

  // Allocate the return value and initialise it
  PROTECT(success = allocVector(LGLSXP,1));
  isuccess = LOGICAL(success);
  *isuccess = 1;

  // Allocate the average vectors,
  // real and imaginary parts separately
  aver = (double*) R_Calloc( nd , double );
  avei = (double*) R_Calloc( nd , double );

  // Initialise to zero
  for( k = 0 ; k < nd ; ++k ){
    aver[ k ] = 0.;
    avei[ k ] = 0.;
  }

  // Allocate vector for data sample counter
  nave = R_Calloc( nd , R_len_t );

  // Initialise to zero
  for( k = 0 ; k < nd ; ++k ) nave[ k ] = 0;

  // Start from begniing of the data vctor
  ind1 = 0;
  ind2 = 0;

  // Search for the start of the first pulse
  while( ( id[ind1] == 0 ) & ( ind1 < nd ) ) ++ind1;
  while( ( id[ind2] == 0 ) & ( ind2 < nd ) ) ++ind2;
  ipp_count = 0;

  // Repeat until end of data
  while( ind2 < nd ){

    // At this point we should be at pulse starts, loop until
    // we hit a point at which both pulses have ended.
    while( id[ind1] | id[ind2]){
      aver[ind1] += cd[ind2].r;
      avei[ind1] += cd[ind2].i;
      ++nave[ind1];
      ++ind1;
      ++ind2;
      if(ind2==nd) break;
    }

    if(ind2==nd) break;

    // Add power values until either of the indices
    // hits the next pulse
    while( (id[ind1]==0) & (id[ind2]==0)){
      aver[ind1] += cd[ind2].r; 
      avei[ind1] += cd[ind2].i;
      ++nave[ind1];
      ++ind1;
      ++ind2;
      if(ind2==nd) break;
    }

    if(ind2==nd) break;

    // Make sure that both indices point to a pulse start,
    // increment if necessary (This takes possible sync 
    // times into account)
    while( ( id[ind1] == 0 ) & ( ind1 < nd ) ) ++ind1;
    while( ( id[ind2] == 0 ) & ( ind2 < nd ) ) ++ind2;

    if(ind2==nd) break;

    // Increment the ipp counter
    ++ipp_count;
    if( ipp_count == ncode ){
      ipp_count = 0;
      ind1 = 0;
      while( id[ind1] == 0 ) ++ind1;
    }
  }

  // Divide the summed values with number of summed pulses
  for( k = 0 ; k < nd ; ++k ){
    if( nave[ k ] ){
      aver[k] /= (double)nave[k];
      avei[k] /= (double)nave[k];
    }
  }


  // Now there are averaged values available for one code
  // cycle, copy the valeus to make furhter analysis
  // simpler. Start from beginning of the data vector.
  ind1 = 0;
  ind2 = 0;

  // Search for the start of the first pulse
  while( ( id[ind1] == 0 ) & ( ind1 < nd ) ) ++ind1;
  while( ( id[ind2] == 0 ) & ( ind2 < nd ) ) ++ind2;
  ipp_count = 0;

  // Repeat until end of data
  while( ind2 < nd ){
    // At this point we should be at pulse starts,
    // loop until both pulses have ended
    while( id[ind1] | id[ind2]){
      cd[ind2].r = aver[ind1];
      cd[ind2].i = avei[ind1];
      ++ind1;
      ++ind2;
      if(ind2==nd) break;
    }

    if(ind2==nd) break;

    // Add power values until either of
    // the indices hits the next pulse
    while( (id[ind1]==0) & (id[ind2]==0)){
      cd[ind2].r = aver[ind1];
      cd[ind2].i = avei[ind1];
      ++ind1;
      ++ind2;
      if(ind2==nd) break;
    }

    if(ind2==nd) break;

    // Make sure that both indices point to a pulse start
    while( ( id[ind1] == 0 ) & ( ind1 < nd ) ) ++ind1;
    while( ( id[ind2] == 0 ) & ( ind2 < nd ) ) ++ind2;

    if(ind2==nd) break;

    // Increment the ipp counter
    ++ipp_count;
    if( ipp_count == ncode ){
      ipp_count = 0;
      ind1 = 0;
      while( id[ind1] == 0 ) ++ind1;
    }
  }

  // Free the temporary vectors
  Free(nave);
  Free(aver);
  Free(avei);

  UNPROTECT(1);

  return( success );

}
