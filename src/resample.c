// file:resample.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"
 
/*
  Resampling with linear interpolation. Reduces to a simple 
  boxcar filter when the filter length is an integer 
  multiple of the original sample interval. 

  Final sample rate must be smaller than or
  equal to the original one. 

  This function overwrites existing data vectors

  Arguments:
   cdata    Complex data samples
   idata    Index vector for cdata
   ndata    Data vector length
   nup      Upsamling factor
   nfilter  Filter length on upsampled data
            (final length is nfilter / nup)
   nfirst   Decimation start index
   nfirstfrac start point within the boxcar filter in upsampled units
   ipartial 0 if partial matched with filter
            should not be accepted in idata vector

  Returns:
   success  1 if resampling was successful, 0 otherwise

*/

SEXP resample( SEXP cdata , SEXP idata , SEXP ndata , SEXP nup , SEXP nfilter , SEXP nfirst , SEXP nfirstfrac, SEXP ipartial )
{

  Rcomplex * restrict cd = COMPLEX(cdata);
  int * restrict id = LOGICAL(idata);
  int nd = *INTEGER(ndata);
  const  int nu = *INTEGER(nup);
  const int nf = *INTEGER(nfilter);
  const int ns = *INTEGER(nfirst);
  const int nsf = *INTEGER(nfirstfrac);
  const int ipar = *LOGICAL(ipartial);
  uint64_t i, j, k, l, m, n;
  double frac;
  Rcomplex tmpsum;
  int tmpi[2];

  // For the return value
  SEXP success;
  int * restrict isuccess;

  // Allocate the return value and initialise it
  PROTECT(success = allocVector(LGLSXP,1));
  isuccess = LOGICAL(success);
  *isuccess = 1;

  /*
    i the current filter start point in upsampled data
    j the current point inside the (upsampled) boxcar filter
    k the current point within the original data vector
    l the current point within the resampled data vector
   */

  i = ns * nu ;  // Starting point in upsampled units
  //  j = nu-1;      // We are originally at the 
  //                 // beginning of the boxcar filter
  j = nsf+nu-1;  // increment with nu-1, we will use the full sample
                 // cd[k] in any case. The first resampled one will
                 // be wrong if nsf/=0, but we could not help this if 
                 // nsf < 0 in any case. 
  k = ns;        // Starting point in original sampling
  l = 0;         // Current point in the final filtered and
                 // decimated data vector, start filling 
                 // from beginning
  tmpsum.r = 0.; // Initialise the temp filter sum to zero
  tmpsum.i = 0.;
  tmpi[0] = 1;
  tmpi[1] = 0;


  while( ( ( i + nf ) / nu ) <= nd ){ // Current filter start + filter length <= data length
    while( j < nf ){         // One filter length of data
      tmpsum.r += cd[ k ].r; // Add the current point to the filter sum
      tmpsum.i += cd[ k ].i;
      tmpi[0] *= id[k];
      tmpi[1] += id[k];
      j += nu;               // Jump  to the next point that actually needs to be calculated
      ++k ;                  // Increment the sample counter of the original data vector
    }
    //    // Fraction of the k'th sample in the original data
    //    // vector that will go to l+1'th resampled point
    //    frac = ( (double)( j - nf + 1 ) ) / (double)nu;
    // not like this, it will create effectively two filters...
    // this should be better
    frac = 0.;
    if( ( j - nf + 1 ) ==  nu ) frac = 1.;

       //
       // the whole fraction thing could be removed, but the above lines will fix 
       // this for the time being. IV 2016-02-16.
       //
       // ... on the other hand, this will be rather easy to convert into upsamling, if that
       // would ever be needed?
       //





    // Now k could be beyond the data vector length,
    // check that it is not
    if( k < nd ){                    
      // Add the fraction that belongs to the k'th point
      tmpsum.r += ( 1. - frac )*cd[k].r;             
      tmpsum.i += ( 1. - frac )*cd[k].i;
      if( frac < .99999 ) tmpi[0] *= id[k];
      if( frac < .99999 ) tmpi[1] += id[k];
      // Now tmpsum is ready, copy its contents to
      // the l'th element of the data vector
      cd[l].r = tmpsum.r;                            
      cd[l].i = tmpsum.i;
      id[l] = ipar ? tmpi[1] : tmpi[0];
      // Put the remaining fraction of
      // k'th sample to the tmpsum
      tmpsum.r = frac*cd[k].r;        
      tmpsum.i = frac*cd[k].i;
      tmpi[0] = ( frac < .00001 ) ? 1 : id[k];
      tmpi[1] = ( frac < .00001 ) ? 0 : id[k];
      // One filter length backwards
      j -= nf;                      
      // The sample where we ended in the previous step was
      //  already added to tmpsum, jump to the next one             
      j += nu;                     
      // Move one filter length forwards 
      /*
      i += nf;                     
      ++k;
      */
      ++l;
    }

    // i and k must be incremented also at end of data to get us out of the loop
    i += nf;                     
    ++k;
  }

  // If we were exactly at end of data frac is unity, we will still get one more sample
  // k was incremented after hitting the end of data
  if( k == ( nd + 1 ) ){
    if( frac > .9999999 ){
      cd[l].r = tmpsum.r;
      cd[l].i = tmpsum.i;
      id[l] = ipar ? tmpi[1] : tmpi[0];
      ++l;
    }
  }

  *(INTEGER(ndata)) = l;

  // remove protection from the return value
  UNPROTECT(1);

  // return the variable success only, the data is now stored
  // in the R vectors 'cdata', 'idatar', and 'idatai'
  return(success);

}


/*
  Resampling with linear interpolation. Reduces to a simple 
  boxcar filter when the filter length is an integer 
  multiple of the original sample interval. 

  Final sample rate must be smaller than 
  or equal to the original one. 

  This function allocates new data vectors

  Arguments:
   cdata    Complex data samples
   idata    Index vector for cdata
   ndata    Data vector length
   nup      Upsamling factor
   nfilter  Filter length on upsampled data (final length
            is nfilter / nup)
   nfirst   Decimation start index
   ipartial 0 if partial matched with filter should not be
            accepted in idata vector

  Returns:
   ans      A list with components:
            cdata    Resampled complex data vector
            idata    Index vector for cdata
            ndata    Data vector length
            success  1 if resampling was successful,
                     0 otherwise

*/


SEXP resample_R( SEXP cdata , SEXP idata , SEXP ndata , SEXP nup , SEXP nfilter , SEXP nfirst , SEXP nfirstfrac, SEXP ipartial)
{
  SEXP ans;
  SEXP cdata_new;
  SEXP idata_new;
  SEXP ndata_new;
  SEXP s;
  SEXP names;
  char *cnames[4] = {"cdata","idata","ndata","success"};
  Rcomplex * restrict cnew;
  Rcomplex * restrict cold;
  int * restrict inew;
  int * restrict iold;
  uint64_t k;
  PROTECT_INDEX cpind=0;
  PROTECT_INDEX ipind=0;


  // Output list ans[[1]] = cdata , ans[[2]] = idata ,
  // ans[[3]] = ndata , ans[[4]] = success
  PROTECT( ans = allocVector( VECSXP , 4 ) );

  // Allocate the new complex vector
  PROTECT_WITH_INDEX( cdata_new = allocVector( CPLXSXP , *(INTEGER(ndata)) ) , &cpind );

  // Allocate the new logical vector
  PROTECT_WITH_INDEX( idata_new = allocVector( LGLSXP , *(INTEGER(ndata)) ) , &ipind );

  // Allocate the new ndata variable
  PROTECT( ndata_new = allocVector( INTSXP , 1 ) );

  // A pointer to the new cdata vector
  cnew = COMPLEX( cdata_new );

  // A pointer to the old cdata vector
  cold = COMPLEX( cdata );

  // A pointer to the new idata vector
  inew = LOGICAL( idata_new );

  // A pointer to the old idata vector
  iold = LOGICAL( idata );

  // Copy data from old cdata to new cdata
  for( k = 0 ; k < *(INTEGER(ndata)) ; ++k ){
    cnew[k].r = cold[k].r;
    cnew[k].i = cold[k].i;
  }

  // Copy data from old idata to new idata
  for( k = 0 ; k < *(INTEGER(ndata)) ; ++k ){
    inew[k] = iold[k];
  }

  // Use the same pointers to copy old ndata to new ndata
  inew = INTEGER( ndata_new);
  iold = INTEGER( ndata );
  *inew = *iold;

  // The  success logical
  PROTECT( s = allocVector( LGLSXP , 1 ) );

  // The actual resampling
  s = resample( cdata_new , idata_new , ndata_new , nup ,  nfilter ,  nfirst , nfirstfrac , ipartial );

  // Reallocate the vectors to match with the new data length
  SET_LENGTH( cdata_new , *INTEGER(ndata_new) );
  REPROTECT( cdata_new , cpind );
  SET_LENGTH( idata_new , *INTEGER(ndata_new) );
  REPROTECT( idata_new , ipind );

  // Collect the data into the return list
  SET_VECTOR_ELT( ans , 0 , cdata_new );
  SET_VECTOR_ELT( ans , 1 , idata_new );
  SET_VECTOR_ELT( ans , 2 , ndata_new );
  SET_VECTOR_ELT( ans , 3 , s );

  // Set the name attributes
  PROTECT( names = allocVector( STRSXP , 4 ));
  SET_STRING_ELT( names , 0 , mkChar( cnames[0] ) );
  SET_STRING_ELT( names , 1 , mkChar( cnames[1] ) );
  SET_STRING_ELT( names , 2 , mkChar( cnames[2] ) );
  SET_STRING_ELT( names , 3 , mkChar( cnames[3] ) );
  setAttrib( ans , R_NamesSymbol , names);

  UNPROTECT(6);

  return(ans);

}


