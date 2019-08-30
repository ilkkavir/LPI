// file:read_gdf_data.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

#include "LPI.h"

/*
  Read IQ data from .gdf files

  This function allocates new data vectors.

  Arguments:
   ndata      Total number of data points to read
   nfiles     Number of data files
   filepaths  nfiles data file paths
   istart     nfiles start indices
   iend       nfiles end indices
   bigendian  logical, 0 if files are little-endian

  Returns:
   ans       A list with elements
              cdata   Complex data samples
             idatar  Lowest bits from real part (PPS)
             idatai  Lowest bits from imaginary part (TX)
             ndata   Data vector length
             success Logical, set if all requested data was read

 */

SEXP read_gdf_data_R( SEXP ndata , SEXP nfiles , SEXP filepaths, SEXP istart , SEXP iend , SEXP bigendian)
{

  SEXP ans;
  SEXP cdata;
  SEXP idatar;
  SEXP idatai;
  SEXP n;
  SEXP s;
  SEXP names;
  char *cnames[5] = {"cdata","idatar","idatai","ndata","success"};

  // Output list ans[[1]] = cdata , ans[[2]] = idatar ,
  // ans[[3]] = idatai , ans[[4]] = ndata , ans[[5]] = success
  PROTECT( ans = allocVector( VECSXP , 5 ) ); 
  // The cdata vector
  PROTECT( cdata = allocVector( CPLXSXP , *(INTEGER(ndata)) ) );
  // The idatar vector
  PROTECT( idatar = allocVector( LGLSXP , *(INTEGER(ndata)) ) );
  // The idatai vector
  PROTECT( idatai = allocVector( LGLSXP , *(INTEGER(ndata)) ) );
  // The ndata vector
  PROTECT( n = allocVector( INTSXP , 1 ) );
  // the read success logical
  PROTECT( s = allocVector( LGLSXP , 1 ) );

  // The actual reading
  s = read_gdf_data( cdata , idatar , idatai , ndata , nfiles, filepaths, istart , iend , bigendian );

  // Collect the data vectors into the list
  SET_VECTOR_ELT( ans , 0 , cdata  );
  SET_VECTOR_ELT( ans , 1 , idatar );
  SET_VECTOR_ELT( ans , 2 , idatai );
  SET_VECTOR_ELT( ans , 3 , ndata  );
  SET_VECTOR_ELT( ans , 4 , s      );

  // Set the name attributes
  PROTECT( names = allocVector( STRSXP , 5 ));
  SET_STRING_ELT( names , 0 , mkChar( cnames[0] ) );
  SET_STRING_ELT( names , 1 , mkChar( cnames[1] ) );
  SET_STRING_ELT( names , 2 , mkChar( cnames[2] ) );
  SET_STRING_ELT( names , 3 , mkChar( cnames[3] ) );
  SET_STRING_ELT( names , 4 , mkChar( cnames[4] ) );
  setAttrib( ans , R_NamesSymbol , names);

  UNPROTECT(7);

  return(ans);

}

/*
  Read IQ data from .gdf files

  This function reads the data to pre-allocated vectors

  Arguments:
   cdata      ndata complex vector for data samples
   idatar     ndata integer vector for lowest bits
              in real part
   idatai     ndata integer vector for lowest bits
              in imaginary part
   ndata      Total number of data points read from files
   nfiles     Number of data files
   filepaths  nfiles data file paths
   istart     nfiles start indices
   iend       nfiles end indices
   bigendian  logical, 0 if files are little-endian

  Returns:
   success    1 if all data was read, 0 otherwise

 */
SEXP read_gdf_data( SEXP cdata , SEXP idatar , SEXP idatai , SEXP ndata , SEXP nfiles , SEXP filepaths , SEXP istart , SEXP iend , SEXP bigendian )
{
  // Pointers to the R variables
  Rcomplex *cd = COMPLEX(cdata);
  int *idr = LOGICAL(idatar);
  int *idi = LOGICAL(idatai);
  int *nd = INTEGER(ndata);
  int *nf = INTEGER(nfiles);
  const char *fpath;
  int *is = INTEGER(istart);
  int *ie = INTEGER(iend);
  int be = *LOGICAL(bigendian);
  // For the return value
  SEXP success;
  int *isuccess;
  // Counters and other temporary variables
  register uint64_t k, n, kd;
  register int16_t ir;
  register int16_t ii;
  uint8_t *rblock;
  FILE *fp;
  size_t rr;

  // Allocate the return value and initialise it
  PROTECT(success = allocVector(LGLSXP,1));
  isuccess = LOGICAL(success);
  *isuccess = 1;

  // Data point counter
  kd = 0;

  // Allocate memory for the data
  rblock = (uint8_t*) malloc(*nd*4);

  for( k=0 ; k<(*nf) ; ++k ) {

    // Select the data file from input list
    fpath = CHAR(STRING_ELT(filepaths,k));

    // Open the data file for reading
    fp = fopen( fpath , "r" );

    // If the fopen failed set success
    // to false and break the read loop
    if(fp == 0){
      *isuccess = 0;
      break;
    }

    // Seek for the starting point
    fseek( fp , (is[k]*4) , SEEK_SET);

    // Read data from file number k
    rr = fread( &(rblock[(kd*4)]) , 1 , ((ie[k]-is[k]+1) * 4) , fp );

    // Success will be set only if every single sample
    // is successfully read
    *isuccess = *isuccess && (rr==(ie[k]-is[k]+1)*4);

    // Close the file
    fclose(fp);

    for( n=is[k] ; n<=ie[k] ; ++n ){

      // Conversion from 8-bit to 16-bit integers
      if( be ){ // big endian
	ir = ( rblock[kd*4] << 8 ) + rblock[kd*4+1];
	ii = ( rblock[kd*4+2] << 8 ) + rblock[kd*4+3];
      }else{ // little endian
	ir = ( rblock[kd*4+1] << 8 ) + rblock[kd*4];
	ii = ( rblock[kd*4+3] << 8 ) + rblock[kd*4+2];
      }
      // Put the lowest bits to the index vectors
      idr[kd] = ir & 0x0001;
      idi[kd] = ii & 0x0001;

      // Remove the lowest bits from the complex data vector
      ir &= 0xfffe;
      ii &= 0xfffe;

      // Conversion to double floats 
      // (which are used in Rcomplex data type)
      cd[kd].r = (double)(ir);
      cd[kd].i = (double)(ii);

      // Increment the sample counter
      ++kd;

    }

  }
  // Free memory
  free(rblock);

  // Remove protection from the return value
  UNPROTECT(1);

  // Return the variable success only, the data is now stored
  // in the R vectors 'cdata', 'idatar', and 'idatai'
  return(success);
  
}

