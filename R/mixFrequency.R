mixFrequency <- function( cdata , ndata , frequency ){
#
# Complex frequency mixing. 
#
# INPUT:
#  cdata     a complex vector of input data
#  ndata     number of samples in the data vector (or number of samples to use from the beginning)
#  frequency mixing frequency, the data vector will be multiplied with a complex sinusoid exp(1i*2*pi*frequency*k),
#            where k is the index in cdata vector, starting from 0
#
# OUPUT:
#  cdata     the cdata vector after frequency mixig
#  success   a logical value TRUE if the frequency mixing was successful, otherwise FALSE
#
#
#
  
  storage.mode(ndata) <- "integer"
  
  return( .Call( "mix_frequency_R" , cdata , ndata , frequency ) )

}
