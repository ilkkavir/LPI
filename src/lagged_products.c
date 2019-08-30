// file:lagged_products.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"
/*
  Calculate lagged products of a signal 
  and its complex conjugate. 

  This function allocates new data vectors

  Arguments:
   cdata1  ndata1 vector of complex signal samples
   cdata2  ndata2 vector of complex signal samples
   idata1  ndata1 integer vector of usable 
           RX sample positions
   idata2  ndata2 integer vector of usable
           RX sample positions
   ndata1  Number of samples in cdata1 and idata1
   ndata2  Number of samples in cdata2 and idata2
   lag     Lag
  
  Returns:
   ans       A list with elements
              cdata   Complex vector of lagged products
              idata   Index vector for cdata
              ndata   Data vector length
              success Logical, set if all processing
                      was successful
*/

SEXP lagged_products_alloc( SEXP cdata1 , SEXP cdata2 , SEXP idata1 , SEXP idata2 , SEXP ndata1 , SEXP ndata2 , SEXP lag)
{
  Rcomplex *cd1 = COMPLEX(cdata1);
  Rcomplex *cd2 = COMPLEX(cdata2);
  int *id1 = LOGICAL(idata1);
  int *id2 = LOGICAL(idata2);
  int *nd1 = INTEGER(ndata1);
  int *nd2 = INTEGER(ndata2);
  int *l = INTEGER(lag);

  SEXP ans;
  SEXP lcdata;
  Rcomplex *lcd;
  SEXP lidata;
  int *lid;
  SEXP success;
  int *isuccess;
  SEXP ndata;
  int *nd;
  SEXP names;
  char *cnames[4] = {"cdata","idata","ndata","success"};
  int k=0;

  // Allocate the return value list
  PROTECT( ans  = allocVector( VECSXP , 4 ) );

  // Allocate the ndata output
  PROTECT( ndata = allocVector( INTSXP , 1 ) );

  // A local pointer to ndata
  nd = INTEGER( ndata );

  // Output data length will be minimum of the two
  //  input data lengths, minus the time-lag
  *nd = *nd1 - *l;
  if( *nd1 > *nd2 ) *nd = *nd2 - *l;

  // Allocate the lagged product vector
  PROTECT( lcdata = allocVector( CPLXSXP , *nd ) );

  // A local pointer to the lagged product vector
  lcd = COMPLEX( lcdata );

  // Allocate an index vector for the lagged products
  PROTECT( lidata = allocVector( LGLSXP , *nd ) );

  // A local pointer to the lagged product vector
  lid = LOGICAL( lidata );

  // Allocate the success return value
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // A local pointer to the success value
  isuccess = LOGICAL( success );
  *isuccess = 1;

  // The actual lagged product calculation
  for( k = 0 ; k < *nd ; ++k ){

    // Calculate the index vector point
    lid[k] = (id1[k] * id2[k+ *l]);

    // Calculate the actual data product only if the index vector is set
    if(lid[k]){
      lcd[k].r = cd1[k].r * cd2[k+ *l].r + cd1[k].i * cd2[k+ *l].i;
      lcd[k].i = -cd1[k].r * cd2[k+ *l].i + cd1[k].i * cd2[k+ *l].r;
    }
  }


  // Collect the return values under the list "ans"
  SET_VECTOR_ELT( ans , 0 , lcdata );
  SET_VECTOR_ELT( ans , 1 , lidata );
  SET_VECTOR_ELT( ans , 2 , ndata );
  SET_VECTOR_ELT( ans , 3 , success );

  // Set the name attributes
  PROTECT( names = allocVector( STRSXP , 4 ) );
  SET_STRING_ELT( names , 0 , mkChar( cnames[0] ) );
  SET_STRING_ELT( names , 1 , mkChar( cnames[1] ) );
  SET_STRING_ELT( names , 2 , mkChar( cnames[2] ) );
  SET_STRING_ELT( names , 3 , mkChar( cnames[3] ) );
  setAttrib( ans , R_NamesSymbol , names );

  UNPROTECT(6);

  return(ans);

}




/*
  Calculate lagged products of a signal
  and its complex conjugate. 

  This function overwrites existing data vectors

  Arguments:
   cdata1  ndata1 vector of complex signal samples
   cdata2  ndata2 vector of complex signal samples
   idata1  ndata1 integer vector of usable
           RX sample positions
   idata2  ndata2 integer vector of usable
           RX sample positions
   cdatap  complex vector for the lagged products
   idatap  integer vector for the lagged product indices
   ndata1  Number of samples in cdata1 and idata1
   ndata2  Number of samples in cdata2 and idata2
   lag     Lag

  Returns:
   success 1 if processing was succesful, 0 otherwise
  
*/

SEXP lagged_products( SEXP cdata1 , SEXP cdata2 , SEXP idata1 , SEXP idata2 , SEXP cdatap ,\
		      SEXP idatap , SEXP ndata1 , SEXP ndata2 , SEXP lag )
{
  Rcomplex *cd1      =  COMPLEX(cdata1);
  Rcomplex *cd2      =  COMPLEX(cdata2);
  int      *id1      =  LOGICAL(idata1);
  int      *id2      =  LOGICAL(idata2);
  Rcomplex *cdp      =  COMPLEX(cdatap);
  int      *idp      =  LOGICAL(idatap);
  int       nd1      = *INTEGER(ndata1);
  int       nd2      = *INTEGER(ndata2);
  int       l        = *INTEGER(lag)   ;
  SEXP      success                    ;
  int      *isuccess                   ;
  int       k        =  0              ;
  int       npr                        ;

  // Output data length will be minimum of the
  //  two input data lengths, minus the time-lag
  npr = nd1 - l;
  if( nd1 > nd2 ) npr = nd2 - l;

  // Allocate the success return value
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // A local pointer to the success value
  isuccess = LOGICAL( success );
  *isuccess = 1;

  // The actual lagged product calculation
  for( k = 0 ; k < npr ; ++k ){

    // The logical vector
    idp[k] = (id1[k] * id2[k+ l]);

    // Multiply the actual data points only
    //  if the logical vector is set
    if(idp[k]){
      cdp[k].r = cd1[k].r * cd2[k+ l].r + cd1[k].i * cd2[k+ l].i;
      cdp[k].i = cd1[k].r * cd2[k+ l].i - cd1[k].i * cd2[k+ l].r;
    }
  }

  // Set the logical vector to false at
  // points where it cannot be calculated
  for( k = 0 ; k < l ; ++k ){
    idp[npr+k] = 0;
  }

  UNPROTECT(1);

  return(success);

}



/*
  Real-valued lagged products for variance estimation.

  No Index vectors, because they are carried with
  the complex vectors.

  This function overwrites existing data vectors

  Arguments:
   rdata1  ndata1 vector of real signal samples
   rdata2  ndata2 vector of real signal samples
   prdata  real vector for the lagged products
   ndata1  Number of samples in rdata1
   ndata2  Number of samples in rdata2
   lag     Lag

  Returns:
   success 1 if processing was successful, 0 otherwise

*/
SEXP lagged_products_r( SEXP rdata1 , SEXP rdata2 , SEXP prdata , SEXP ndata1 ,\
			SEXP ndata2 , SEXP lag )
{
  double *rd1      =  REAL(rdata1)   ;
  double *rd2      =  REAL(rdata2)   ;
  double *prd      =  REAL(prdata)   ;
  int     nd1      = *INTEGER(ndata1);
  int     nd2      = *INTEGER(ndata2);
  int     l        = *INTEGER(lag)   ;
  SEXP    success                    ;
  int    *isuccess                   ;
  int     k        =  0              ;
  int     npr                        ;

  // Output data length will be minimum of the two input
  // data lengths, minus the time-lag
  npr = nd1 - l;
  if( nd1 > nd2 ) npr = nd2 - l;

  // Allocate the success return value
  PROTECT( success = allocVector( LGLSXP , 1 ) );

  // A local pointer to the success value
  isuccess = LOGICAL( success );
  *isuccess = 1;

  // The actual lagged product calculation
  for( k = 0 ; k < npr ; ++k ){
    prd[k] =  rd1[k] * rd2[k+ l];
  }

  UNPROTECT(1);

  return(success);

}
