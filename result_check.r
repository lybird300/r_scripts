#Script to check results from analyses

args <- commandArgs(trailing=TRUE)
#args[[1]] = result table
#args[[2]] = gzipped flag
#args[[3]] = trait

#read the table
result_table <- args[[1]]
gz_flag <- args[[2]]
trait <- args[[3]]

result_table
gz_flag
trait

#load pval check function
source("/nfs/users/nfs_m/mc14/Work/r_scripts/pval_check.r")
# load plot formatter
source("/nfs/users/nfs_m/mc14/Work/r_scripts/qqman.R")

#create the corrected result and plots
corrected_result <- pvalcheck(result_table,gz_flag,trait)

#print the new table
gz <- gzfile(paste(result_table,"corrected.txt.gz",sep="_"), open="w"))
	write.table(corrected_result$pvals, gz, row.names=F, quote=F, sep=" ",col.names=T) 
close(gz)
	
