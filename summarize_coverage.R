#Summarize coverage for a single bam file

#read commandline args
args <- commandArgs(trailing=TRUE)

#conversion of impute files to mach format files
#Args:
#args[[1]]= input file
#Every file name is intended to be used with ABSOLUTE PATH

input_file <- args[[1]]

cov_info <- read.table(input_file, header=FALSE,sep="\t")
colnames(cov_info) <- c("CHROM","POS","COVERAGE")
cov_info$CHROM <- as.factor(cov_info$CHROM)
cov_info$POS <- as.factor(cov_info$POS)


write.table(summary(cov_info),file=paste(input_file,'summary.txt',sep='_'), sep="\t", row.names=FALSE, col.names=TRUE, quote=F)



#q(save="yes")

