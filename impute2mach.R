#we need to use the local libs
library('GenABEL')
library('DatABEL')

#read commandline args
args <- commandArgs(trailing=TRUE)

#conversion of impute files to mach format files
#Args:
#args[[1]]= chr
#args[[2]]= geno file
#args[[3]]= info file
#args[[4]]= sample file
#args[[5]]= output file
#Every file name is intended to be used with ABSOLUTE PATH

chr <- args[[1]]
geno_file <- args[[2]]
info_file <- args[[3]]
sample_file <- args[[4]]
output_file <- args[[5]]

impute2mach(geno_file, info_file, sample_file, output_file, maketextdosefile=TRUE,order_info_snp_pos_freq1_info_qual_type=c(2,3,6:9))

#q(save="yes")
