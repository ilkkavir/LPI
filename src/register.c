// file:register.c
// (c) 2010- University of Oulu, Finland
// Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
// Licensed under FreeBSD license.

// R registration of C functions

#include "LPI.h"
static const R_CallMethodDef callMethods[20] = {
  { "read_gdf_data_R"       , (DL_FUNC) & read_gdf_data_R       , 6 } , 
  { "mix_frequency_R"       , (DL_FUNC) & mix_frequency_R       , 3 } , 
  { "index_adjust_R"        , (DL_FUNC) & index_adjust_R        , 3 } , 
  { "lagged_products_alloc" , (DL_FUNC) & lagged_products_alloc , 7 } ,
  { "lagged_products"       , (DL_FUNC) & lagged_products       , 9 } ,
  { "lagged_products_r"     , (DL_FUNC) & lagged_products_r     , 6 } ,
  { "fishs_add"             , (DL_FUNC) & fishs_add             , 7 } ,
  { "theory_rows_alloc"     , (DL_FUNC) & theory_rows_alloc     , 13} ,
  { "theory_rows"           , (DL_FUNC) & theory_rows           , 17} ,
  { "prepare_data"          , (DL_FUNC) & prepare_data          , 10} ,
  { "average_power"         , (DL_FUNC) & average_power         , 6 } ,
  { "deco_add"              , (DL_FUNC) & deco_add              , 7 } ,
  { "average_profile"       , (DL_FUNC) & average_profile       , 4 } ,
  { "dummy_add"             , (DL_FUNC) & dummy_add             , 10} ,
  { "resample"              , (DL_FUNC) & resample              , 8 } ,
  { "resample_R"            , (DL_FUNC) & resample_R            , 8 } ,
  { "range_ambiguity"       , (DL_FUNC) & range_ambiguity       , 9 } ,
  { "clutter_meas"          , (DL_FUNC) & clutter_meas          , 9 } ,
  { "clutter_subtract"      , (DL_FUNC) & clutter_subtract      , 8 } ,
  { NULL , NULL , 0 }
};

void R_init_LPI(DllInfo *info)
{
  R_registerRoutines( info , NULL , callMethods , NULL , NULL );
}


