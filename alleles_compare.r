#09/19/2012

#R script for extract matching alleles
#
rm(list=ls())

#read commandline args
args <- commandArgs(trailing=TRUE)
 
#Args list:
#args[[1]]: first map file
#args[[2]]: second map file
#args[[3]]: chr

first_filepath <- args[[1]]
second_filepath <- args[[2]]
chr <- args[[3]]

#Upload chr data:
first_file <- read.table(first_filepath,skip=1,header=F,sep='\t', stringsAsFactors=F, comment.char="")
second_file <- read.table(second_filepath,header=F,sep='\t', stringsAsFactors=F, comment.char="")

#take only the first 5 columns
first_file <- first_file[,1:5]
second_file <- second_file[,1:5]

#set colnames
colnames(first_file) <- c('CHROM','POS','ID','REF','ALT')
colnames(second_file) <- c('CHROM','POS','ID','REF','ALT')

#Merge the fileset and output if the alt and ref allele match
merged_file <- merge(first_file,second_file,by=c('CHROM','POS'),sort=T)

colnames(merged_file) <- c('CHROM','POS','F_ID','F_REF','F_ALT','S_ID','S_REF','S_ALT')

#now check if the alleles match
matching <- merged_file[which(merged_file$F_REF == merged_file$S_REF & merged_file$F_ALT == merged_file$S_ALT),]
matching <- cbind(matching,MATCH='TRUE')

#non matching lines
non_matching <- merged_file[-which(merged_file$F_REF == merged_file$S_REF & merged_file$F_ALT == merged_file$S_ALT),]
non_matching <- cbind(non_matching,MATCH='FALSE')

#concat files
merged <- rbind(matching,non_matching)

merged <- merged[order(merged$POS),]

#now write the output of this merging...and summary
write.table(summary(merged),file=paste("merged_summary_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(merged,file=paste("merged_chr",chr,".map",sep=""),sep="\t",col.names=T,quote=F,row.names=F)





