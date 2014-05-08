#Script to calculate genotype discordance
#We have all files in plink format
#
rm(list=ls())
require(gplots)

#read commandline args
args <- commandArgs(trailing=TRUE)

#set args
#ref_alt_file <- '370K/110_unrel_overlapping_biallelic_SNP_CR_HWE_MISS.frq'
#ref_map_file <- '370K/110_unrel_overlapping_biallelic_SNP_CR_HWE_MISS.map'
ref_alt_file <- args[[1]]
ref_map_file <- args[[2]]

#Upload a file used to establish what is the Ref and what the Alt Allele.
#first we need to remove all spaces and convert to tabs from the frq file..
#sed -i 's/ \+/        /g' biallelic_overl.SNP.unfilt.geno.seq.VB.frq

#now we read the map file to join position and rsid and minor allele column
ref_allele_tbl <- read.table(ref_alt_file, sep="\t",header=T, stringsAsFactors=F, comment.char="")
ref_map_tbl <- read.table(ref_map_file, sep="\t",header=F, stringsAsFactors=F, comment.char="")
#now we can remove useless columns from ref_allele
ref_allele_tbl$A1 <- NULL
ref_allele_tbl$NCHROBS <- NULL
ref_map_tbl$V3 <- NULL
gc()
gc()

colnames(ref_map_tbl) <- c("CHR","rsID","POS")
#now merge this file to create a set that has to be used for set the reference allele in files from SEQ or 2.5M imputation
merge_tbl <- merge(ref_allele_tbl,ref_map_tbl, by.x="SNP",by.y="rsID",sort=F)
ref_tbl <- cbind(merge_tbl[,c(2)],merge_tbl[,c(6)],merge_tbl[,c(1)],as.character(merge_tbl[,c(3)]),merge_tbl[,c(4)])
colnames(ref_tbl) <- c("CHR","POS","SNP","MAJOR","MAF")
#write the ref table
write.table(ref_tbl,file="ref_discordance_table.txt",sep="\t",col.names=T,quote=F,row.names=F)


