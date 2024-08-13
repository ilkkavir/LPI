// file:LPI.h
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

// Data types and function prototypes

#include <R.h>
#include <math.h>
#include <stdint.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Complex.h>
#include <R_ext/Constants.h>

//static const double pi=3.1415926535;
#define AMB_N_INTERP  5


// gdf file input
SEXP read_gdf_data_R( SEXP ndata , SEXP nfiles , SEXP filepaths , SEXP istart , SEXP iend , SEXP bigendian);
SEXP read_gdf_data( SEXP cata , SEXP idatar , SEXP idatai , SEXP ndata , SEXP nfiles, SEXP filepaths , SEXP istart , SEXP iend , SEXP bigendian);

// Frequency mixing
SEXP mix_frequency_R( SEXP cdata , SEXP ndata , SEXP frequency);
SEXP mix_frequency( SEXP cdata , SEXP ndata , SEXP frequency);

// Index adjustments
SEXP index_adjust_R( SEXP idata , SEXP ndata , SEXP shifts );
SEXP index_adjust( SEXP idata , SEXP ndata , SEXP shifts );

// Lagged products
SEXP lagged_products_alloc( SEXP cdata1 , SEXP cdata2 , SEXP idata1 , SEXP idata2 , SEXP ndata1 , SEXP ndata2 , SEXP lag);
SEXP lagged_products( SEXP cdata1 , SEXP cdata2 , SEXP idata1 , SEXP idata2 , SEXP cdatap , SEXP idatap , SEXP ndata1 , SEXP ndata2 , SEXP lag );
SEXP lagged_products_r( SEXP rdata1 , SEXP rdata2 , SEXP prdata , SEXP ndata1 , SEXP ndata2 , SEXP lag );

// Theory matrix construction
SEXP theory_rows_alloc( SEXP camb , SEXP iamb , SEXP cprod , SEXP iprod , SEXP rvar , SEXP ndata , SEXP ncur , SEXP nend , SEXP rlims , SEXP nranges , SEXP fitsize , SEXP background, SEXP remoterx ); 
SEXP theory_rows( SEXP camb , SEXP iamb , SEXP cprod , SEXP iprod , SEXP rvar , SEXP ndata , SEXP ncur , SEXP nend , SEXP rlims , SEXP nranges , SEXP arows , SEXP irows , SEXP mvec , SEXP mvar , SEXP nrows , SEXP background, SEXP remoterx );

// Inverse problem solvers
SEXP fishs_add( SEXP Qvec , SEXP yvec , const SEXP arows , const SEXP irows , const SEXP meas , const SEXP var , const SEXP nx , const SEXP nrow );
SEXP deco_add( SEXP Qvec , SEXP yvec , const SEXP arows , const SEXP irows , const SEXP meas , const SEXP var , const SEXP nx , const SEXP nrow );
SEXP dummy_add( SEXP msum , SEXP vsum , SEXP rmin , SEXP rmax , SEXP mdata , SEXP mambig , SEXP iamb , SEXP iprod , SEXP edata , SEXP ndata );

// All data preparations collected together
SEXP prepare_data( SEXP cdata , SEXP idata , SEXP ndata , SEXP frequency, SEXP shifts , SEXP nup , SEXP nfilter , SEXP nfirst , SEXP nfirstfrac , SEXP ipartial );

// Average signal power in points withe identical IPPs and pulse lengths
SEXP average_power( SEXP cdata , SEXP idatatx , SEXP idatarx , SEXP ndata , SEXP maxrange , SEXP nminave);

// Average lag profile
SEXP average_profile( SEXP cdata , SEXP idata , SEXP ndata , SEXP N_CODE);

// Resampling
SEXP resample( SEXP cdata , SEXP idata , SEXP ndata , SEXP nup , SEXP nfilter , SEXP nfirst , SEXP nfirstfrac , SEXP ipartial);
SEXP resample_R( SEXP cdata , SEXP idata , SEXP ndata , SEXP nup , SEXP nfilter , SEXP nfirst , SEXP nfirstfrac , SEXP ipartial);

// Range ambiguity function calculation with optional interpolation
SEXP range_ambiguity( SEXP cdata1 ,SEXP cdata2 , SEXP idata1 , SEXP idata2 , SEXP cdatap , SEXP idatap , SEXP ndata1 ,  SEXP ndata2 ,  SEXP lag );

// Ground clutter suppression
SEXP clutter_meas( const SEXP tcdata , const SEXP tidata , const SEXP rcdata , const SEXP ridata , const SEXP ndata , const SEXP rmin ,  const SEXP rmax , SEXP Qvec , SEXP yvec );
SEXP clutter_subtract( const SEXP tcdata , const SEXP tidata , SEXP rcdata , const SEXP ridata , const SEXP ndata , const SEXP rmin , const SEXP rmax , const SEXP cldata );
void fishs_add_clutter( SEXP Qvec , SEXP yvec , Rcomplex * arow , Rcomplex * meas , const int nx );
		      
