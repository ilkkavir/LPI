// file:theory_rows.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.


#include "LPI.h"

/*
  Make theory matrix rows and measurement vectors.

  This function overwrites existing data vectors
  
  Arguments:
   camb        Complex range ambiguity functions
   iamb        Index vector of range ambiguity functions
   cprod       Complex lagged product vector
   iprod       Index vector of lagged products
   rvar        Measurement variance vector
   ndata       Data vector length
   ncur        Current sample index
   nend        Last sample index to use (in this call)
   rlims       Range gate limits
   nranges     Number of range gates
   arows       Complex theory rows
   irows       Theory row indices
   mvec        Inversion measurement vector
   mvar        Inversion measurement variances
   nrows       Number of theory rows produced during 
               this call
   background  0 if additional background term is not used
   remoterx    0 if measurements TX times should not be used

  Returns:
   success     0 if no theory rows were produced _and_ end of
               data was reached, 1 otherwise
 */


SEXP theory_rows_r( SEXP camb , SEXP iamb , SEXP cprod , SEXP iprod , SEXP rvar , SEXP ndata , SEXP ncur , SEXP nend , SEXP rlims , SEXP nranges , SEXP arowsR , SEXP arowsI , SEXP irows , SEXP mvecR , SEXP mvecI , SEXP mvar , SEXP nrows , SEXP background , SEXP remoterx )
{
  const Rcomplex * restrict amb = COMPLEX(camb);
  const int * restrict amb_i = LOGICAL(iamb);
  const Rcomplex * restrict prod =  COMPLEX(cprod);
  const int * restrict prod_i = LOGICAL(iprod);
  const double * restrict var =  REAL(rvar);
  int n_cur = *INTEGER(ncur);
  int n_end = *INTEGER(nend);
  const int * restrict r_lims = INTEGER(rlims);
  const int n_ranges = *INTEGER(nranges);
  const int n_data = *INTEGER(ndata);
  const int bg = *LOGICAL(background);
  const int remrx = *LOGICAL(remoterx);
  double * restrict aR = REAL(arowsR);
  double * restrict aI = REAL(arowsI);
  int * restrict i_rows = LOGICAL(irows);
  double * restrict mR = REAL(mvecR);
  double * restrict mI = REAL(mvecI);
  double * restrict m_var = REAL(mvar);
  SEXP success;
  int * restrict i_success;
  int n_rows;
  R_len_t k;
  R_len_t n_start;
  R_len_t i;
  R_len_t j;
  R_len_t subi;
  R_len_t addi;
  R_len_t gati;
  int r_min;
  int r_lim;
  int r_max;
  int r_cur;

  
  // Check that n_end <= n_data
  n_end = ( n_data > n_end ? n_end : n_data );
  
  // Check that n_cur <= n_data
  n_cur = ( n_data > n_cur ? n_cur : n_data );

  // Success output
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // Local pointer to the success output
  i_success = LOGICAL( success );

  // Set the success output
  *i_success = 1;

  // The lowest range gate limnit - 1
  r_min = r_lims[0] - 2 ;

  // Samples with non-zero range ambiguity
  // function at heights below r_lim
  // will not be used in the theory matrix
  // Initialize r_min for monostatic reception
  r_lim = r_min;
  //  -1 (all samples accpected) for remote reception
  if( remrx ) r_lim = -1;

  // The highest range gate limit
  r_max = r_lims[n_ranges] + 1;

  // Make the first theory row.
  n_start = n_cur;
  // If we are too close to start of data
  // skip points as necessary
  if( n_start < r_lims[ n_ranges ] ) n_start = r_lims[ n_ranges ];

  // Make sure that we did not yet pass the end point
  if( n_start < n_end ){
    // Go through all range-gates
    for( i = 0 ; i <  n_ranges ; ++i ){
      // Initialize the theory matrix to zero
      aR[i] = .0;
      aI[i] = .0;
      i_rows[i] = 0;

      // Add contribution from all ranges
      // integrated to this gate
      for( j = r_lims[i] ; j < r_lims[ i + 1 ] ; ++j ){

        // In amb_i == 0 points there might be erroneous
	// values from previously calculated lags, 
        // it is thus extremely important to check
	// amb_i before addition / subtraction!
        if(amb_i[ n_start - j ]){
          aR[i] += amb[ n_start - j ].r;
          aI[i] += amb[ n_start - j ].i;
          i_rows[i] += amb_i[ n_start - j ];
        }
      }
    }

    // The last gate will be 1 or 0, depending on whether
    // the background ACF will be suppressed or not.
    aR[ n_ranges ] = ( bg == 0 ? 0.0 : 1.0);
    aI[ n_ranges ] = 0.0;
    i_rows[ n_ranges ]   = ( bg == 0 ? 0 : 1 );

  // If the first row could not be formed
  // set success to false and return
  }else{
    *i_success = 0;
  }

  // From this point on all possible theory rows will  be
  // formed but only those with indprod set are stored,
  // others are immediately overwritten

  // Number of stored rows
  n_rows = 0;

  // Range from the latest pulse
  r_cur = r_max;
  for( k = (n_start-r_max) ; k < n_start ; ++k ){
    if( k >= 0 ){
      if(amb_i[k]){
        r_cur = 0;
      }else{
        ++r_cur;
      }
    }
  }

  // Use all data points from n_start to n_end
  for( k = n_start ; k < n_end ; ++k ){

    // If this data point will be used (!=0 for clarity,
    // the prod_i vector may contains values larger than 1)
    if( (prod_i[k] != 0) & (r_cur > r_lim) & (r_cur < r_max)){

      // Copy data to the measurement vector
      mR[n_rows] = prod[k].r;
      mI[n_rows] = prod[k].i;
      m_var[n_rows]   = var[k];

      // Copy the current theory vectors to the next one.
      for( i = 0 ; i <  ( n_ranges + 1 ) ; ++i ){
       	i_rows[ ( n_rows + 1 ) * ( n_ranges + 1 ) + i ]   = i_rows[ n_rows * ( n_ranges + 1 ) + i ];
        // Set the theory rows exactly to zero at points
	// where the index vector is zero. This makes 
	// identification of blind ranges much easier.
        if(i_rows[ n_rows  * ( n_ranges + 1 ) + i ]==0){
          aR[ ( n_rows + 1 ) * ( n_ranges + 1 ) + i ] = 0.0;
          aI[ ( n_rows + 1 ) * ( n_ranges + 1 ) + i ] = 0.0;
          aR[ n_rows * ( n_ranges + 1 ) + i ] = 0.0;
          aI[ n_rows * ( n_ranges + 1 ) + i ] = 0.0;
        // Otherwise copy the theory matrix row
        }else{
          aR[ ( n_rows + 1 ) * ( n_ranges + 1 ) + i ] = aR[ n_rows * ( n_ranges + 1 ) + i ];
          aI[ ( n_rows + 1 ) * ( n_ranges + 1 ) + i ] = aI[ n_rows * ( n_ranges + 1 ) + i ];
        }
      }

      // Increment the theory row counter
      ++n_rows;

    }

    // Now form the next theory row using the previous
    // one and the range limit indices
    for( i = 0 ; i < n_ranges ; ++i ){
      // Index in the theory matrix 
      // (that is stored as a vector)
      gati = n_rows * ( n_ranges + 1 ) + i;
      // Index of the data point that 
      // will be added to this gate
      addi = k - r_lims[i] + 1;
      // Index of the data point that
      // will be subtracted from this gate
      subi = k - r_lims[i+1] + 1;

      // Do additions / subtractions only if the point
      // contains a non-zero ambiguity value
      if( amb_i[ addi ] ){
        aR[ gati ] += amb[ addi ].r;
        aI[ gati ] += amb[ addi ].i;
        i_rows[ gati ]   += amb_i[ addi ];
      }
      if( amb_i[ subi ] ){
        aR[ gati ] -= amb[ subi ].r;
        aI[ gati ] -= amb[ subi ].i;
        i_rows[ gati ]   -= amb_i[ subi ];
      }
	
    }

    // Count samples to exclude everything that contains
    // echoes from below the first gate
    if( amb_i[ k ] ){
      r_cur = 0;
    }else{
      ++r_cur;
    }

  }

  // Write the row count to the output variable
  *( INTEGER( nrows ) ) = n_rows;

  // Update the current position in the data vector
  *( INTEGER( ncur ) ) = n_end;

  
  UNPROTECT(1);

  return(success);

}
