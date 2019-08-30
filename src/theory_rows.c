// file:theory_rows.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.


#include "LPI.h"
/*
  Make theory matrix rows and measurement vectors.

  This function allocates new data vectors.
  
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
   fitsize     0 if the vectors should not be reallocated to
               match the final data size.
   background  0 if additional background term is not used
   remoterx    0 if measurements TX times should not be used

  
  Returns:
   ans       A list with elements
              arows   Theory matrix rows
              irows   Theory row indices
              m       Inversion measurement vector
              var     Measurement variances
              nrows   Number of theory rows produced
              success Logical, set if all processing
                      was successful

 */

SEXP theory_rows_alloc( SEXP camb , SEXP iamb , SEXP cprod , SEXP iprod , SEXP rvar , SEXP ndata , SEXP ncur , SEXP nend  , SEXP rlims , SEXP nranges , SEXP fitsize , SEXP background , SEXP remoterx )
{
  const int n_cur = *INTEGER(ncur);
  const int n_end = *INTEGER(nend);
  const int n_ranges = *INTEGER(nranges);
  const int fit_size = *LOGICAL(fitsize);
  SEXP ans;
  SEXP arows;
  SEXP irows;
  SEXP mvec;
  SEXP mvar;
  SEXP success;
  SEXP nrows;
  SEXP names;
  int n_rows;
  const char * c_names[6] =  {"arows","irows","m","var","nrows","success"};
  PROTECT_INDEX arind =  0;
  PROTECT_INDEX irind = 0;
  PROTECT_INDEX mind = 0;
  PROTECT_INDEX vind = 0;


  // Output list
  PROTECT( ans = allocVector( VECSXP , 5 ) );

  // A vector for the theory matrix rows
  PROTECT_WITH_INDEX( arows = allocVector( CPLXSXP , ( ( n_end - n_cur + 1 ) * ( n_ranges + 1) ) ) , & arind );

  // A vector for the theory matrix indices
  PROTECT_WITH_INDEX( irows = allocVector( LGLSXP , ( ( n_end - n_cur + 1 ) * ( n_ranges + 1) ) ) , & irind );

  // A vector for the measurements
  PROTECT_WITH_INDEX( mvec = allocVector( CPLXSXP , ( n_end - n_cur + 1 ) ) , & mind );

  // A vector for the measurement errors
  PROTECT_WITH_INDEX( mvar = allocVector( REALSXP , ( n_end - n_cur + 1 ) ) , & vind );

  // Number of rows for the R output
  PROTECT( nrows = allocVector( INTSXP , 1 ) );

  // Success output
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // Call the theory_rows function to actually make the rows
  success = theory_rows( camb , iamb , cprod , iprod , rvar , ndata , ncur , nend , rlims ,\
                         nranges , arows , irows , mvec , mvar , nrows , background , remoterx );

  // Read the row count
  n_rows = *(INTEGER(nrows));

  // Reallocate the vectors to match with the data lengths
  if(fit_size){
    SET_LENGTH( arows , ( n_rows * ( n_ranges + 1 ) ) );
    REPROTECT( arows , arind );
    SET_LENGTH( irows , ( n_rows * ( n_ranges + 1 ) ) );
    REPROTECT( irows , irind );
    SET_LENGTH( mvec , n_rows );
    REPROTECT( mvec , mind );
    SET_LENGTH( mvar , n_rows );
    REPROTECT( mvar , vind );
  }

  // Collect the data into the return list
  SET_VECTOR_ELT( ans , 0 , arows   );
  SET_VECTOR_ELT( ans , 1 , irows   );
  SET_VECTOR_ELT( ans , 2 , mvec    );
  SET_VECTOR_ELT( ans , 3 , mvar    );
  SET_VECTOR_ELT( ans , 4 , nrows   );
  SET_VECTOR_ELT( ans , 5 , success );

  // Set the names attributes
  PROTECT( names = allocVector( STRSXP , 5 ) );
  SET_STRING_ELT( names , 0 , mkChar( c_names[0] ) );
  SET_STRING_ELT( names , 1 , mkChar( c_names[1] ) );
  SET_STRING_ELT( names , 2 , mkChar( c_names[2] ) );
  SET_STRING_ELT( names , 3 , mkChar( c_names[3] ) );
  SET_STRING_ELT( names , 4 , mkChar( c_names[4] ) );
  SET_STRING_ELT( names , 5 , mkChar( c_names[5] ) );
  setAttrib( ans , R_NamesSymbol , names);


  UNPROTECT(7);

  return(ans);

}





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


SEXP theory_rows( SEXP camb , SEXP iamb , SEXP cprod , SEXP iprod , SEXP rvar , SEXP ndata , SEXP ncur , SEXP nend , SEXP rlims , SEXP nranges , SEXP arows , SEXP irows , SEXP mvec , SEXP mvar , SEXP nrows , SEXP background , SEXP remoterx )
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
  Rcomplex * restrict a_rows = COMPLEX(arows);
  int * restrict i_rows = LOGICAL(irows);
  Rcomplex * restrict m_vec = COMPLEX(mvec);
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
      a_rows[i].r = .0;
      a_rows[i].i = .0;
      i_rows[i] = 0;

      // Add contribution from all ranges
      // integrated to this gate
      for( j = r_lims[i] ; j < r_lims[ i + 1 ] ; ++j ){

        // In amb_i == 0 points there might be erroneous
	// values from previously calculated lags, 
        // it is thus extremely important to check
	// amb_i before addition / subtraction!
        if(amb_i[ n_start - j ]){
          a_rows[i].r += amb[ n_start - j ].r;
          a_rows[i].i += amb[ n_start - j ].i;
          i_rows[i] += amb_i[ n_start - j ];
        }
      }
    }

    // The last gate will be 1 or 0, depending on whether
    // the background ACF will be suppressed or not.
    a_rows[ n_ranges ].r = ( bg == 0 ? 0.0 : 1.0);
    a_rows[ n_ranges ].i = 0.0;
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
  for( k = (n_start-r_min) ; k < n_start ; ++k ){
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
      m_vec[n_rows].r = prod[k].r;
      m_vec[n_rows].i = prod[k].i;
      m_var[n_rows]   = var[k];

      // Copy the current theory vectors to the next one.
      for( i = 0 ; i <  ( n_ranges + 1 ) ; ++i ){
       	i_rows[ ( n_rows + 1 ) * ( n_ranges + 1 ) + i ]   = i_rows[ n_rows * ( n_ranges + 1 ) + i ];
        // Set the theory rows exactly to zero at points
	// where the index vector is zero. This makes 
	// identification of blind ranges much easier.
        if(i_rows[ n_rows  * ( n_ranges + 1 ) + i ]==0){
          a_rows[ ( n_rows + 1 ) * ( n_ranges + 1 ) + i ].r = 0.0;
          a_rows[ ( n_rows + 1 ) * ( n_ranges + 1 ) + i ].i = 0.0;
          a_rows[ n_rows * ( n_ranges + 1 ) + i ].r = 0.0;
          a_rows[ n_rows * ( n_ranges + 1 ) + i ].i = 0.0;
        // Otherwise copy the theory matrix row
        }else{
          a_rows[ ( n_rows + 1 ) * ( n_ranges + 1 ) + i ].r = a_rows[ n_rows * ( n_ranges + 1 ) + i ].r;
          a_rows[ ( n_rows + 1 ) * ( n_ranges + 1 ) + i ].i = a_rows[ n_rows * ( n_ranges + 1 ) + i ].i;
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
        a_rows[ gati ].r += amb[ addi ].r;
        a_rows[ gati ].i += amb[ addi ].i;
        i_rows[ gati ]   += amb_i[ addi ];
      }
      if( amb_i[ subi ] ){
        a_rows[ gati ].r -= amb[ subi ].r;
        a_rows[ gati ].i -= amb[ subi ].i;
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
